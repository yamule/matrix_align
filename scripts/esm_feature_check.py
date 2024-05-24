import torch
import esm
import sys,re,os,gzip
import copy;
import numpy as np;
import gc;
import random;

# タンパク質配列を ESM にかける際に、長いタンパク質についても分割して処理し、切断部周辺はいくつかスキップして重複部分は平均 Representation として
# 出力するスクリプト




import argparse;

def check_bool(v):
    v = v.lower();
    if v == "true" or v == "1":
        return True;
    if v == "false" or v == "0":
        return False;
    raise Exception("true or false or 1 or 0 are expected.");
        

parser = argparse.ArgumentParser();
parser.add_argument("--infile",help='Multi-FASTA フォーマットのファイル',required= True) ;
parser.add_argument("--outfile",required= True) ;
parser.add_argument("--random_seed",required=False,default=None,type=int) ;
parser.add_argument("--max_length",required=False,default=1000,help='これより長い配列は前後どちらかが最初の段階で削られる。',type=int);
parser.add_argument("--min_fragment_length",required=False,default=20,help='前後を削ってこれより短い場合は検証されない。',type=int);
parser.add_argument("--num_samples",required=False,default=10000,help='fasta ファイルから何本使用するか。',type=int);
parser.add_argument("--device",required=True);
parser.add_argument("--model_path",required=True);
parser.add_argument("--normalize",required=False,default=False,type=check_bool);
parser.add_argument("--batch_size",required=False,default=20,type=int);

args = parser.parse_args();

print(args);
model_path = args.model_path;
max_length = args.max_length; # esm に渡す文字列の最大長
batch_size = args.batch_size; # model に与える配列数
num_samples = args.num_samples;
min_fragment_length = args.min_fragment_length;
normalize = args.normalize;

if args.random_seed is not None:
    random.seed(123);
    
infile = args.infile;
outfile = args.outfile;

ddev = args.device;

# Load ESM-2 model
model, alphabet = esm.pretrained.load_model_and_alphabet(model_path)
if ddev == "cuda":
    model = model.eval().cuda();
else:
    model = model.eval();

batch_converter = alphabet.get_batch_converter()

def loadFasta(filename):
    if filename.endswith("gz"):
        fin = gzip.open(filename,"rt");
    else:
        fin = open(filename,"r");
    ret = [];
    cdict = dict();
    cdict["seq"] = "";
    ret.append(cdict);

    for ll in fin:
        mat = re.search("[\s]*>",ll);
        if(not mat == None):
            cdict = dict();
            ret.append(cdict);
            nmat = re.search("[\s]*>[\s]*([^\s]+)",ll);
            if(not nmat == None):
                cdict["name"] = nmat.group(1);
                cdict["desc"] = "";
            dmat = re.search("[\s]*>[\s]*([^\s]+)[\s]+([^\s][^\r\n]*)",ll);
            if(not dmat == None):
                cdict["desc"] = dmat.group(2);
            cdict["seq"] = "";
        else:
            cdict["seq"] += re.sub("[\s]","",ll).upper();
            
    if(len(ret[0]["seq"]) == 0):
        ret.pop(0);
    fin.close();
    return ret;


# sys.stderr.write(str(alphabet.to_dict())+"\n");
i_to_a = {};
for kk,vv in alphabet.to_dict().items():
    i_to_a[vv] = kk; # 重複は気にしない


fass_ = loadFasta(infile);
fass = [];
for ff in list(fass_):
    sseq = re.sub(r"[^A-Z]","",ff["seq"].upper());
    if len(sseq) < min_fragment_length+1: # 1 残基は削るので
        continue;
    sseq = re.sub(r"[JUZBOX]","X",sseq);
    fass.append(ff);

samples = [];
used = {};
if num_samples > len(fass):
    sys.stderr.write("Number of sequences are lower than num_samples. {} vs {}\n".format(len(fass),num_samples));
    samples.extend(fass);
else:
    while len(samples) < num_samples:
        ppin = random.randrange(0,len(fass));
        if ppin in used:
            continue;
        samples.append(fass[ppin]);
        used[ppin] = 100;
del used;
del fass;
with open(outfile+".samples","wt") as fout:
    for ss in samples:
        fout.write(">"+ss["name"]+" "+ss["desc"]+"\n");
        fout.write(ss["seq"]+"\n");

name_used = {};
for ss in list(samples):
    nname = ss["name"];
    counter = 0;
    while nname in name_used:
        counter += 1;
        nname = ss["name"]+".".str(counter);
    name_used[nname] = 100;
    ss["name"] = nname;
if normalize:
    normalizerfile = "tmp."+str(random.random())+".dat";
    normalizeout = open(normalizerfile,"wt");
value_diff_sum = None;
while len(samples) > 0:
    data = [];
    print("remained",len(samples));
    while len(data) < batch_size:
        if len(samples) == 0:
            break;
        ff = samples.pop();
        sseq = ff["seq"];
        if len(sseq) > max_length:
            if random.random() > 0.5:
                data.append([
                    ff["name"],sseq[:max_length]
                ]);
            else:
                data.append([
                    ff["name"],sseq[-max_length:]
                ]);
        else:
            data.append([
                ff["name"],sseq
            ]);
    shifter = [
        (0,0), (1,0), (25, 0), (50, 0), (100,  0), (200,  0),
               (1,1), (25,25), (50,50), (100,100), (200,200),
               (0,1), ( 0,25), ( 0,50), (  0,100), (  0,200),
    ];

    #shifter = [
    #    (0,0), (1,0),
    #           (1,1),
    #           (0,1)
    #];

    base_value = {};
    for sii,sss in enumerate(list(shifter)):
        data_fragment = [];
        for dd in list(data):
            sseq = dd[1];
            if len(sseq) < sss[0]+sss[1]+min_fragment_length:
                continue;
            if sss[1] != 0:
                data_fragment.append(
                    [dd[0],dd[1][sss[0]:-sss[1]]]
                );
            else:
                data_fragment.append(
                    [dd[0],dd[1][sss[0]:]]
                );
        if len(data_fragment) < 1:
            continue;
        batch_labels, batch_strs, batch_tokens = batch_converter(data_fragment)
        batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
        data_fragment.clear();
        batch_tokens = torch.utils._pytree.tree_map(lambda x:x.to(ddev),batch_tokens);
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=False);
        token_representations = results["representations"][33];
        del results;
        # Generate per-sequence representations via averaging
        # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
        lines = [];
        for i, tokens_len in enumerate(batch_lens):
            ress = token_representations[i, 1 : tokens_len - 1].to("cpu").numpy().astype(np.float64);
            if sii == 0:
                assert sss[0] == 0 and sss[1] == 0;
                base_value[batch_labels[i]] = ress;
                for ii in range(ress.shape[0]):
                    buffs = [];
                    for jj in range(ress.shape[1]):
                        buffs.append("{:.7}".format(ress[ii,jj]));
                    lines.append("\t".join(buffs));
            else:
                for ii in range(ress.shape[0]):
                    diff = np.abs(base_value[batch_labels[i]][sss[0]+ii] - ress[ii]);
                    if value_diff_sum is None:
                        value_diff_sum = np.zeros_like(diff,dtype=np.float64);
                    value_diff_sum += diff;
        normalizeout.write("\n".join(lines));
        normalizeout.write("\n");
if normalize:
    normalizeout.close();
    vsum = np.zeros_like(value_diff_sum,dtype=np.float64);
    vcount = 0;
    with open(normalizerfile)as fin:
        for ll in fin:
            ll = re.sub(r"[\r\n]","",ll);
            if len(ll) == 0:
                continue;
            ptt = re.split(r"\t",ll);
            assert len(ptt) == vsum.shape[0];
            vcount += 1;
            for ii in range(len(ptt)):
                vsum[ii] += float(ptt[ii]);
    vave = vsum/float(vcount);
    vvar = np.zeros_like(value_diff_sum,dtype=np.float64);
    with open(normalizerfile)as fin:
        for ll in fin:
            ll = re.sub(r"[\r\n]","",ll);
            if len(ll) == 0:
                continue;
            ptt = re.split(r"\t",ll);
            assert len(ptt) == vsum.shape[0];
            for ii in range(len(ptt)):
                bf = float(ptt[ii]);
                vvar[ii] += (bf-vave[ii])*(bf-vave[ii]);
    vvar = vvar/float(vcount);
    vstd = np.sqrt(vvar);
    #os.remove(normalizerfile);
    with open(outfile+".stats","wt") as fout:
        for ii in range(vvar.shape[0]):
            fout.write("index:\t{}\tsum:\t{}\tvar:\t{}\tcount:{}\n".format(ii,vsum[ii],vvar[ii],vcount));

    res = [];
    for ii in range(value_diff_sum.shape[0]):
        if vstd[ii] == 0.0:
            sys.stderr.write("STD  in column "+str(ii+1)+" is zero???\n");
            res.append([ii,value_diff_sum[ii]]);
        else:
            res.append([ii,value_diff_sum[ii]/vstd[ii]]);
            
    res = list(sorted(res,key=lambda x:x[1]));
    with open(outfile,"wt") as fout:
        for rr in list(res):
            fout.write(str(rr[0])+"\t"+str(rr[1])+"\n");

else:
    res = [];
    for ii in range(value_diff_sum.shape[0]):
        res.append([ii,value_diff_sum[ii]]);
    res = list(sorted(res,key=lambda x:x[1]));
    with open(outfile,"wt") as fout:
        for rr in list(res):
            fout.write(str(rr[0])+"\t"+str(rr[1])+"\n");

