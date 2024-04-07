import torch
import esm
import sys,re,os,gzip
import copy;
import numpy as np;

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
parser.add_argument("--crop_length",required= True,help='断片化後の長さ',type=int);
parser.add_argument("--model_path",required= True);
parser.add_argument("--shift_length",required= True,help='次の断片を作る際の移動量',type=int);
parser.add_argument("--cut_length",required= True,help='断片化された際の境界の残基についてはいくつか削除する',type=int);
parser.add_argument("--device",required=True);
parser.add_argument("--batch_size",required=False,default=20,type=int);
parser.add_argument("--round",required=False,default=7,help='小数点以下で丸める際の桁数',type=int);

args = parser.parse_args();

print(args);
model_path = args.model_path;
crop_length = args.crop_length; # esm に渡す文字列の最大長
shift_length = args.shift_length; # 複数に分割される際の開始点の移動量
batch_size = args.batch_size; # model に与える配列数
cut_length = args.cut_length; # 複数に分割された際に分割点に近い部分のデータをどれくらい捨てるか

infile = args.infile;
outfile = args.outfile;
rounder = args.round;
ddev = args.device;
assert outfile.endswith(".gz"), "Output file must ends with .gz ."

# Load ESM-2 model
model, alphabet = esm.pretrained.load_model_and_alphabet(model_path)
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

# crop_length に分割し、(開始位置, 文字列) のタプルで返す
def crop_sequence(aastr,crop_length,shift_length):
    aastr = re.sub(r"[\s]","",aastr.upper());
    assert re.search(r"[^A-Z]",aastr) is None, "Unexpected letter "+aastr;
    assert "J" not in aastr;

    if len(aastr) <= crop_length:
        return [(0,aastr),];
    ret = [];
    offset = 0;
    seqlen = len(aastr);

    while True:
        enn = offset+crop_length;
        if enn > seqlen:
            enn = seqlen;
        offset = enn - crop_length;
        ret.append((offset,copy.deepcopy(aastr[offset:enn])));
        if enn == seqlen:
            break;
        offset += shift_length;
    return ret;

# (offset,representations) というリストを与えて所定のパラメータを基準に結合する
def merge_representations(replist_,crop_length,cut_length):
    if len(replist_) == 1:
        return replist_[0][1];
    replist = list(sorted(replist_,key=lambda x:x[0])) # 0 番目の要素にはオフセット位置が入っており、これを基準にソートする
    repnum = len(replist);
    full_length = replist[-1][0] + crop_length;
    counter = np.zeros((full_length,));
    values = np.zeros((full_length,replist[0][1].shape[-1]));
    for ii in range(repnum):
        st = replist[ii][0];
        en = st+crop_length;
        pst = 0;
        pen = crop_length;
        
        if ii != 0:
            st += cut_length;
            pst += cut_length;
        if ii != repnum-1:
            en -= cut_length;
            pen -= cut_length;
        
        values[st:en] += replist[ii][1][pst:pen];
        counter[st:en] += 1;
    
    return values/counter[:,None];

# sys.stderr.write(str(alphabet.to_dict())+"\n");
i_to_a = {};
for kk,vv in alphabet.to_dict().items():
    i_to_a[vv] = kk; # 重複は気にしない

"""
デバッグ用
crop_length = 600; # esm に渡す文字列の最大長
shift_length = 150; # 複数に分割される際の開始点の移動量
cut_length = 75; # 複数に分割された際に分割点に近い部分のデータをどれくらい捨てるか
batch_size = 20; # model に与える配列数


crop_length = 25; # esm に渡す文字列の最大長
shift_length = 12; # 複数に分割される際の開始点の移動量
cut_length = 4; # 複数に分割された際に分割点に近い部分のデータをどれくらい捨てるか
outdir = "nogit/";
rounder = 7;
"""

fass = loadFasta(infile);

# Prepare data (first 2 sequences from ESMStructuralSplitDataset superfamily / 4)
#fass = [ # デバッグ用
#    {"name":"protein1", "seq":"MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"},
#]

full_seq = {};
seq_desc = {};
for ff in fass:
    if ff["name"] in full_seq:
        sys.stderr.write("Duplicated name found. "+ff["name"]+"\n");
        raise Exception();
    full_seq[ff["name"]] = ff["seq"];
    seq_desc[ff["name"]] = ff["desc"];
buff = {};
seq_fragment = {};
remained = [];
data = [];
if outfile.endswith(".gz"):
    fout = gzip.open(outfile,"wt");
else:
    fout = open(outfile,"wt");
format_string = None;
replen = None;
while len(fass) > 0 or len(remained) > 0:
    while len(data) < batch_size:
        if len(remained) == 0:
            break;
        data.append(remained.pop());

    while len(data) < batch_size:
        if len(fass) == 0:
            break;
        ff = fass.pop();
        cropped = crop_sequence(ff["seq"],crop_length,shift_length);
        seq_fragment[ff["name"]] = {};
        for cc in cropped:
            offset,seqq = cc;
            fragmentname = ff["name"]+"##"+str(offset);
            seq_fragment[ff["name"]][fragmentname] = [offset,None]
            if len(data) < batch_size:
                data.append([
                    fragmentname,seqq
                ]);
            else:
                remained.append([
                    fragmentname,seqq
                ]);
    

    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
    data.clear();
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
    token_representations = results["representations"][33]

    # Generate per-sequence representations via averaging
    # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
    for i, tokens_len in enumerate(batch_lens):
        mat = re.search(r"(.+)##([^#]+)$",batch_labels[i]);
        assert mat is not None;
        sname = mat.group(1);
        seq_fragment[sname][batch_labels[i]][1] = copy.deepcopy(token_representations[i, 1 : tokens_len - 1].to("cpu").numpy());

    completed = [];
    for seqname in list(seq_fragment.keys()):
        flag = True;
        for jj in list(seq_fragment[seqname].keys()):
            if seq_fragment[seqname][jj][1] is None:
                flag = False;
                break;
        if flag:
            plist = [];
            for jj in list(seq_fragment[seqname].keys()):
                plist.append(seq_fragment[seqname][jj]);
            completed.append(seqname);
            repres = merge_representations(plist,crop_length,cut_length);
            print("processed:",seqname,"fragments:",len(plist),flush=True);
            #print("\t".join([str(x) for x in repres]));
            replen_ = len(repres[0]);
            if replen is None:
                replen = replen_;
                format_string = ("\t{:."+str(rounder)+"}")*replen;
            else:
                assert replen == replen_;
            linebuff = [];
            linebuff.append(">"+seqname+" "+seq_desc[seqname]+"\n");
            
            for jj in range(len(repres)):
                linebuff.append(full_seq[seqname][jj]);
                linebuff.append(format_string.format(*repres[jj]))
                #for mm in range(replen):
                #    linebuff.append(format_string.format(repres[jj][mm]));
                linebuff.append("\n");
            
            fout.write("".join(linebuff));
            del linebuff;
    for cc in completed:
        seq_fragment.pop(cc);

fout.close();

assert len(seq_fragment) == 0;
assert len(data) == 0;
assert len(remained) == 0;
assert len(fass) == 0;

