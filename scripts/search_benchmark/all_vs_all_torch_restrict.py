import re,os,sys,gzip;
import numpy as np;
import argparse;
import subprocess;
import torch;
import math;
import shutil;
import gc;
import random;


# 一部のみ使用して global per channel normalization をかける
# 検証用

EPSILON = 1.0e-7;

def check_bool(v):
    v = v.lower();
    if v == "true" or v == "1":
        return True;
    if v == "false" or v == "0":
        return False;
    raise Exception("true or false or 1 or 0 are expected.");
        
parser = argparse.ArgumentParser();
parser.add_argument("--targetdir",required= True) ;
parser.add_argument("--outdir",required=True) ;
parser.add_argument("--device",required=False,default="cuda") ;
parser.add_argument("--num_samples",required=False,default=None,type=int) ;
parser.add_argument("--sample_superfamily",required=False,default=None) ;
parser.add_argument("--pigz",required=False,default=True,type=check_bool) ;
parser.add_argument("--batch_size",required=False,default=20,type=int) ;
parser.add_argument("--random_seed",required=False,default=123,type=int) ;


# 使用しない
# parser.add_argument("--unbiased_global_stats",required=False,default=False,type=check_bool) ;
# parser.add_argument("--check_fragment",required=False,default=True,type=check_bool) ; 

args = parser.parse_args();
print(args);

random.seed(args.random_seed);

# どちらか一方のみが None
assert (args.sample_superfamily is None) != (args.num_samples is None);

targetdir = args.targetdir;
outdir = args.outdir;
use_pigz = args.pigz;
batch_size=args.batch_size;
ddev = torch.device(args.device);
use_unbiased_global_stats = False;
check_fragment = False;
global_per_channel_normalization = True;
num_samples = args.num_samples;
sample_superfamily = args.sample_superfamily;

EPSILON_TENSOR=torch.tensor(EPSILON, dtype=torch.float32, device=ddev)

if use_pigz and shutil.which("pigz") is None:
    raise Exception("pigz is not found.");

if not torch.cuda.is_available() and ddev.type == 'cuda':
    raise Exception("CUDA is not available.");

if not os.path.exists(outdir):
    os.mkdir(outdir);
else:
    raise Exception("Please remove "+outdir);

def check_first_line(infile):
    ffin = gzip.open(infile,"rt");
    l = re.sub(r"[\r\n]","",ffin.readline());
    ffin.close();
    assert l[0] == ">";

    mat = re.search(r">([^\s]+)[\s]+([^\s].+)",ll);
    if mat:    
        return {"name":mat.group(0),"desc":mat.group(1)};
    else:
        raise Exception("????"+infile+"\n"+l+"\n");
    

def load_mat(infile):
    ret = [];
    # {name, desc, seq, value}
    with gzip.open(infile,"rt") as fin:
        current = None;
        for ll in fin:
            ll = re.sub(r"[\s]+$","",ll);
            if len(ll) == 0:
                continue;
            if ll.startswith("#"):
                continue;
            if ll.startswith("//"):
                if current is not None:
                    ret.append(current);
                current = None;
                continue;
            if ll.startswith(">"):
                if current is not None:
                    ret.append(current);
                    current = None;
                mat = re.search(r">([^\s]+)[\s]+([^\s].+)",ll);
                current = {};
                if mat:
                    current["name"] = mat.group(1);
                    current["desc"] = mat.group(2);
                else:
                    mat = re.search(r">([^\s]+)",ll);
                    if mat:
                        current["name"] = mat.group(1);
                        current["desc"] = "";
                    else:
                        raise Exception("Unexpected line "+ll);
                current["seq"] = [];
                current["value"] = [];
                continue;
            ptt = re.split(r"[\s]+",ll);
            assert len(ptt[0]) == 1, "Unexpected line "+ll 
            current["seq"].append(ptt[0]);
            current["value"].append(
                torch.tensor([float(xx) for xx in ptt[1:]],dtype=torch.float32,device=ddev)
            );

        if current is not None:
            ret.append(current);
            current = None;
    return ret;

allfiles_ = list(sorted(os.listdir(targetdir)));
allfiles = [];
for aa in allfiles_:
    if aa.endswith("mat.gz"):
        allfiles.append(
            os.path.join(targetdir,aa)
        );
    else:
        sys.stderr.write(aa+" was skipped.\n");
stats = None;
if use_unbiased_global_stats:
    statsfile = os.path.join(outdir,"stats.unbiased.dat");
else:
    statsfile = os.path.join(outdir,"stats.dat");

vsiz = None;
if global_per_channel_normalization:
    ssum = [];
    mmean = [];
    vvar = [];
    headseq = None;
    valcount = 0;
    if num_samples is not None:
        random.shuffle(allfiles);
    usednames = {};
    saout = open(statsfile+".samples","wt");
    for aa in list(allfiles):
        if sample_superfamily is not None:
            chk = check_first_line(aa);
            familyname = re.split(r"[\s]+",chk["desc"])[0];
            pcc = re.split(r"\.",familyname)
            sfname = pcc[0]+"."+pcc[1]+"."+pcc[2];
            if sfname != sample_superfamily:
                continue;

        c = load_mat(aa);
        if len(c) != 1:
            raise Exception(aa+" only one entry per file is expected.");
        if vsiz is None:
            vsiz = len(c[0]["value"][0]);
            headseq = c[0];
            ssum = torch.zeros((vsiz,),dtype=torch.float32,device=ddev)
        else:
            assert len(c[0]["value"][0]) == vsiz, "Inconsistent value sizes detected."+aa+"\n";
        for cc in c:
            saout.write(cc["name"]+" "+cc["desc"]+"\n");
            assert cc["name"] not in usednames;
            usednames[cc["name"]] = 100;

            for jj in range(len(cc["value"])):
                valcount += 1;
                assert vsiz == len(cc["value"][jj]),cc["name"]+" position "+str(jj)+" has different value length with "+headseq["name"] +"\n"+str(len(cc["value"][jj]))+" vs "+str(vsiz);
                ssum += cc["value"][jj];

            if num_samples is not None and len(usednames) >= num_samples:
                break;
        if num_samples is not None and len(usednames) >= num_samples:
            break;
    
    assert len(usednames) > 10; # 適当

    mmean = ssum/float(valcount);
    vvar = torch.zeros_like(mmean);
    for aa in list(allfiles):
        if sample_superfamily is not None:
            chk = check_first_line(aa);
            familyname = re.split(r"[\s]+",chk["desc"])[0];
            pcc = re.split(r"\.",familyname)
            sfname = pcc[0]+"."+pcc[1]+"."+pcc[2];
            if sfname != sample_superfamily:
                continue;
        
        c = load_mat(aa);
        if len(c) != 1:
            raise Exception(aa+" only one entry per file is expected.");
        for cc in c:
            if cc["name"] not in usednames:
                continue;
            assert usednames[cc["name"]] == 100;

            usednames[cc["name"]] = 50;
            for jj in range(len(cc["value"])):
                vc = cc["value"][jj];
                vvar += (mmean -vc)*(mmean-vc);

    for kk in list(usednames.keys()):
        assert usednames[kk] == 50;

    if use_unbiased_global_stats:
        vvar /= float(valcount-1);
    else:
        vvar /= float(valcount);

    sstd = torch.sqrt(vvar);
    saout.close();
    with open(statsfile,"wt") as fout:
        for ii in range(vsiz):
            fout.write(
                "index:\t{}\tvar:\t{:.7f}\tmean:\t{:.7f}\tcount:\t{}\n".format(ii,float(vvar[ii]),float(mmean[ii]),valcount)
            );

    for ii in range(sstd.shape[0]):
        if sstd[ii] == 0: # stdev が 0 の場合、あとの処理で mean が引かれて全部 0 になるはず。
            sstd[ii] = 1;
    stats = {
        "std":sstd,
        "var":vvar,
        "mean":mmean
    };

allfiles = list(sorted(allfiles));

allvalues_average = [];
allvalues_max = [];
allvalues_min = [];

allvalues_median = [];
allvalues_05 = [];
allvalues_95 = [];

numfiles = len(allfiles);
name_desc = [];

basename_to_value = {};
name_to_index = {};
for ii in range(numfiles):
    aa = allfiles[ii];
    c = load_mat(aa);
    if vsiz is None:
        vsiz = len(c[0]["value"][0]);

    for cc in list(c):
        
        basename = cc["name"];
        fragmentindex = 0;
        familyname = re.split(r"[\s]+",cc["desc"])[0]; # 最初のカラムに Family ID が入っている想定
        if check_fragment:
            mat = re.search(r"^(.+)%([0-9]+)$",cc["name"]);
            if mat:
                basename = mat.group(1);
                fragmentindex = int(mat.group(2));

        if basename not in name_to_index:
            basename_to_value[basename] = [];
            nameindex = len(name_desc);
            name_to_index[basename] = nameindex;
            name_desc.append(
                (basename,familyname) 
            );
        
        assert name_desc[name_to_index[basename]][1] == familyname,cc["name"]+" "+cc["desc"]+"\n"+str(name_desc[name_to_index[basename]])+" ???";
        values_all = [];
        for spos in range(len(cc["value"])):
            vv = cc["value"][spos]
            assert len(vv) == vsiz, "Inconsistent value sizes detected."+aa+"\n";
            if stats is not None:
                assert (stats["std"] != 0.0).all(); # 0 の場合は 1 が入っており、mean を引いて 0 になるはず
                xvalue = (vv-stats["mean"])/stats["std"];
            else:
                xvalue = vv;
            values_all.append(xvalue);
        values_all = torch.stack(values_all,dim=0);
        values_all = torch.permute(values_all,(1,0));

        ave = [];
        mmax = [];
        mmin = [];
        v05 = [];
        v95 = [];
        mmed = [];
        qspan = torch.tensor([0.0,0.05,0.5,0.95,1.0],dtype=torch.float32,device=ddev);
        for vii in range(vsiz):
            ave.append(
                float(values_all[vii].mean())
            );

            qres = torch.quantile(input=values_all[vii], q=qspan);
            mmin.append(float(qres[0]));
            v05.append(float(qres[1]));
            mmed.append(float(qres[2]));
            v95.append(float(qres[3]));
            mmax.append(float(qres[4]));
            del qres;
        del values_all;
        basename_to_value[basename].append([[nameindex,fragmentindex,False]
        ,{"av":ave,"ma":mmax,"mi":mmin,"me":mmed,"5":v05,"95":v95}]);
    del c;
    gc.collect();

tagkeys = ["av","ma","mi","me","5","95"];
allvalues_source = {};
for tt in list(tagkeys):
    allvalues_source[tt] = [];

globalid_to_basedata = [];
for kk in list(basename_to_value.keys()):
    vlist = list(sorted(basename_to_value[kk],key=lambda x:x[0][1]));
    baseid = name_to_index[kk];
    tmpp = []
    for vv in list(vlist):
        for tt in list(tagkeys):
            allvalues_source[tt].append(vv[1][tt]);
        tmpp.append(vv[0])
    globalid_to_basedata.extend(tmpp);
    globalid_to_basedata[-1][-1] = True;# 最後のフラグメントは True

del basename_to_value;
gc.collect();

def dot_product(a,b):
    assert len(a.shape) == 2;
    assert len(b.shape) == 2;
    return (a*b).sum(dim=-1);

def cos_sim(a,b):
    assert len(a.shape) == 2;
    assert len(b.shape) == 2;
    anorm = torch.sqrt((a*a).sum(dim=-1, keepdim=True));
    bnorm = torch.sqrt((b*b).sum(dim=-1, keepdim=True));

    if (anorm < EPSILON).any() or (bnorm < EPSILON).any():
        sys.stderr.write("Warning: Extremely low values were found in cos_sim.\n"+str(a)+"\n"+str(b)+"\n");

    amask = anorm > 0.0;
    bmask = bnorm > 0.0;
    anorm = torch.where(amask, anorm, EPSILON_TENSOR);
    bnorm = torch.where(bmask, bnorm,EPSILON_TENSOR);
    ret =  torch.where(torch.squeeze(amask*bmask,dim=-1),((a/anorm)*(b/bnorm)).sum(dim=-1), torch.tensor(0.0, dtype=torch.float32, device=ddev));
    return ret;

def euc_dist(a,b):
    assert len(a.shape) == 2;
    assert len(b.shape) == 2;
    c = a-b;
    return torch.sqrt((c*c).sum(dim=-1));

def euc_dist_norm(a,b):
    assert len(a.shape) == 2;
    assert len(b.shape) == 2;
    anorm = torch.sqrt((a*a).sum(dim=-1, keepdim=True));
    bnorm = torch.sqrt((b*b).sum(dim=-1, keepdim=True));
    
    if (anorm < EPSILON).any() or (bnorm < EPSILON).any():
        sys.stderr.write("Warning: Extremely low values were found in euc_dist_norm.\n"+str(a)+"\n"+str(b)+"\n");

    anorm = torch.where(anorm == 0.0, EPSILON_TENSOR, anorm);
    bnorm = torch.where(bnorm == 0.0, EPSILON_TENSOR, bnorm);
    c = (a/anorm)-(b/bnorm);
    return torch.sqrt((c*c).sum(dim=-1));

def correl(a, b):
    assert len(a.shape) == 2;
    assert len(b.shape) == 2;
    mean_a = a.mean(dim=1, keepdim=True);
    mean_b = b.mean(dim=1, keepdim=True);
    a_centered = a - mean_a;
    b_centered = b - mean_b;
    numerator = (a_centered * b_centered).sum(dim=1);
    denominator = torch.sqrt((a_centered ** 2).sum(dim=1) * (b_centered ** 2).sum(dim=1));
    
    if (denominator < EPSILON).any():
        sys.stderr.write("Warning: Extremely low values were found in correl.\n"+str(a)+"\n"+str(b)+"\n");

    zero_denominator = denominator == 0;
    denominator = torch.where(zero_denominator, EPSILON_TENSOR, denominator);
    return torch.where(zero_denominator, torch.tensor(0.0, dtype=torch.float32, device=ddev),  numerator / denominator);

fragment_counter = len(globalid_to_basedata);
for (stag,ttag) in [
    ("average","av")
    ,("max","ma")
    ,("min","mi")
    ,("v05","5")
    ,("v95","95")
    ,("median","me")
    ]:
    print("calc",stag,flush=True);
    allvalues = torch.tensor(allvalues_source[ttag],dtype=torch.float32,device=ddev);
    del allvalues_source[ttag];
    gc.collect();

    funcs = [
        ("dot",dot_product,True),("cos_sim",cos_sim,True),("euc_dist",euc_dist,False),("euc_dist_norm",euc_dist_norm,False),("correl",correl,True)
    ];
    res = {};
    for ff in list(funcs):
        res[ff[0]] = [];
    prev_index = -1;
    for fragmentindex in range(fragment_counter):
        currenttargetindex = globalid_to_basedata[fragmentindex][0];
        if currenttargetindex != prev_index:
            # 初期化されていない場合エラーを発生させて終了する
            assert len(res[funcs[0][0]]) == 0, "Error in code.";
        prev_index = currenttargetindex;

        query_name = name_desc[currenttargetindex][0];
        query_desc = name_desc[currenttargetindex][1];
        arr_i = allvalues[fragmentindex];
        
        num_batches = math.ceil(fragment_counter/batch_size);
        arr_i_expanded = arr_i.unsqueeze(0).repeat(batch_size, 1);

        with torch.no_grad():
            for jj in range(num_batches):
                end_index = min((jj+1)*batch_size, fragment_counter);
                arr_j = allvalues[jj*batch_size:end_index];
                current_siz = arr_j.shape[0];
                for tag,func, _  in list(funcs):
                    batch_res = func(arr_i_expanded[:current_siz],arr_j).detach().cpu().tolist();
                    for kkk in range(current_siz):
                        globalindex = jj*batch_size+kkk;
                        if globalid_to_basedata[globalindex][0] == currenttargetindex:
                            continue;
                        res[tag].append((globalid_to_basedata[globalindex][0],float(batch_res[kkk])));
                    del batch_res;

        if globalid_to_basedata[fragmentindex][-1]:
            for tag,func,reverser in list(funcs):
                outname = os.path.join(outdir,"res_"+str(currenttargetindex)+"."+stag+"."+tag+".dat");

                assert not os.path.exists(outname),outname+" already exists!";

                with open(outname,"wt") as fout:
                    fout.write(">"+query_name+" "+query_desc+"\n");
                    t = list(sorted(res[tag],key=lambda x:x[1],reverse=reverser));
                    processed = {};
                    for tt in list(t):
                        sequenceindex = tt[0];
                        if sequenceindex in processed:
                            continue;
                        assert sequenceindex != currenttargetindex;
                        processed[sequenceindex] = 100;
                        fout.write("{}\t{}\t{:.7f}".format(name_desc[sequenceindex][0],name_desc[sequenceindex][1],tt[1])+"\n");
                if use_pigz:
                    subprocess.run(["pigz",outname],check=True);
            # 前の配列の全フラグメントが処理されたので初期化
            res = {};
            for ff in list(funcs):
                res[ff[0]] = [];
    del allvalues;
    gc.collect();
    assert globalid_to_basedata[fragment_counter-1][-1],"???";