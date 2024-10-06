import re,os,sys,gzip;
import numpy as np;
import argparse;
import subprocess;
import torch;
import math;
import shutil;

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
parser.add_argument("--global_per_channel_normalization",required=True,type=check_bool) ;
parser.add_argument("--pigz",required=False,default=True,type=check_bool) ;
parser.add_argument("--batch_size",required=False,default=20,type=int) ;
parser.add_argument("--unbiased_global_stats",required=False,default=False,type=check_bool) ;

args = parser.parse_args();

targetdir = args.targetdir;
global_per_channel_normalization = args.global_per_channel_normalization;
outdir = args.outdir;
use_pigz = args.pigz;
use_unbiased_global_stats = args.unbiased_global_stats;
batch_size=args.batch_size;
ddev = torch.device(args.device);

EPSILON_TENSOR=torch.tensor(EPSILON, dtype=torch.float32, device=ddev)

if use_pigz and shutil.which("pigz") is None:
    raise Exception("pigz is not found.");

if not torch.cuda.is_available() and ddev.type == 'cuda':
    raise Exception("CUDA is not available.");

if not os.path.exists(outdir):
    os.mkdir(outdir);

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
            ptt = re.split(r"[\s]",ll);
            assert len(ptt[0]) == 1, "Unexpected line "+ll 
            current["seq"].append(ptt[0]);
            current["value"].append(
                [float(xx) for xx in ptt[1:]]
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
if global_per_channel_normalization:
    ssum = [];
    mmean = [];
    vvar = [];
    vsiz = None;
    headseq = None;
    valcount = 0;
    for aa in list(allfiles):
        c = load_mat(aa);
        if len(c) == 0:
            raise Exception(aa+" does not have data.");
        if vsiz is None:
            vsiz = len(c[0]["value"][0]);
            headseq = c[0];
            for _ in range(vsiz):
                ssum.append(0);
                mmean.append(0);
                vvar.append(0);
        else:
            assert len(c[0]["value"][0]) == vsiz, "Inconsistent value sizes detected."+aa+"\n";
        for cc in c:
            for jj in range(len(cc["value"])):
                valcount += 1;
                assert vsiz == len(cc["value"][jj]),cc["name"]+" position "+str(jj)+" has different value length with "+headseq["name"] +"\n"+str(len(cc["value"][jj]))+" vs "+str(vsiz);
                for ii in range(vsiz):
                    ssum[ii] += cc["value"][jj][ii];

    for ii in range(vsiz):
        mmean[ii] = ssum[ii]/float(valcount);

    for aa in list(allfiles):
        c = load_mat(aa);
        for cc in c:
            for jj in range(len(cc["value"])):
                for ii in range(vsiz):
                    vvar[ii] += (mmean[ii]-cc["value"][jj][ii])*(mmean[ii]-cc["value"][jj][ii]);
    sstd = [];
    for ii in range(vsiz):
        if use_unbiased_global_stats:
            vvar[ii] /= valcount-1;
        else:
            vvar[ii] /= valcount;
        sstd.append(
            np.sqrt(vvar[ii])
        );
    stats = {
        "std":sstd,
        "var":vvar,
        "mean":mmean
    };
    with open(statsfile,"wt") as fout:
        for ii in range(vsiz):
            fout.write(
                "index:\t{}\tvar:\t{:.7f}\tmean:\t{:.7f}\tcount:\t{}\n".format(ii,vvar[ii],mmean[ii],valcount)
            );

allvalues_average = [];
allvalues_max = [];
allvalues_min = [];
numfiles = len(allfiles);
name_desc = [];

sample_counter = 0;
for ii in range(numfiles):
    aa = allfiles[ii];
    c = load_mat(aa);
    if vsiz is None:
        vsiz = len(c[0]["value"][0]);
    else:
        assert len(c[0]["value"][0]) == vsiz, "Inconsistent value sizes detected."+aa+"\n";

    for cc in list(c):
        name_desc.append(
            (cc["name"],re.split(r"[\s]+",cc["desc"])[0]) # 最初のカラムに Family ID が入っている想定
        );
        sample_counter += 1;
        ssum = [0.0 for jj in range(vsiz)];
        mmax = [None for jj in range(vsiz)];
        mmin = [None for jj in range(vsiz)];
        for vv in list(cc["value"]):
            for vii in range(vsiz):
                
                if stats is not None:
                    if stats["std"][vii] == 0:
                        xvalue = vv[vii]-stats["mean"][vii];
                    else:
                        xvalue = (vv[vii]-stats["mean"][vii])/stats["std"][vii];
                else:
                    xvalue = vv[vii];
                
                ssum[vii] += xvalue;
                if mmax[vii] is None:
                    mmax[vii] = xvalue;
                    mmin[vii] = xvalue;
                else:
                    mmax[vii] = max([mmax[vii],xvalue]);
                    mmin[vii] = min([mmin[vii],xvalue]);

        for vii in range(vsiz):
            ssum[vii] /= len(cc["value"]);
        allvalues_average.append(ssum);
        allvalues_max.append(mmax);
        allvalues_min.append(mmin);


for (stag,val) in [("average",allvalues_average),("max",allvalues_max),("min",allvalues_min)]:
    houtname = os.path.join(outdir,stag+".dat");
    with open(houtname,"wt") as fout:
        for ii in range(sample_counter):
            fout.write("\t".join(name_desc[ii])+"\t"+"\t".join(["{:.7f}".format(float(x)) for x in val[ii]])+"\n");

allvalues_average = torch.tensor(allvalues_average,dtype=torch.float32,device=ddev);
allvalues_max = torch.tensor(allvalues_max,dtype=torch.float32,device=ddev);
allvalues_min = torch.tensor(allvalues_min,dtype=torch.float32,device=ddev);

def dot_product(a,b):
    assert len(a.shape) == 2;
    assert len(b.shape) == 2;
    return (a*b).sum(dim=-1);

def cos_sim(a,b):
    assert len(a.shape) == 2;
    assert len(b.shape) == 2;
    anorm = torch.sqrt((a*a).sum(dim=-1, keepdim=True));
    bnorm = torch.sqrt((b*b).sum(dim=-1, keepdim=True));
    amask = anorm > 0.0;
    bmask = bnorm > 0.0;
    anorm = torch.where(amask, anorm, EPSILON_TENSOR);
    bnorm = torch.where(bmask, bnorm,EPSILON_TENSOR);
    return torch.where(amask*bmask,((a/anorm)*(b/bnorm)).sum(dim=-1), torch.tensor(0.0, dtype=torch.float32, device=ddev));

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
    zero_denominator = denominator == 0;
    denominator = torch.where(zero_denominator, EPSILON_TENSOR, denominator);
    return torch.where(zero_denominator, torch.tensor(0.0, dtype=torch.float32, device=ddev),  numerator / denominator);

    
for (stag,allvalues) in [("average",allvalues_average),("max",allvalues_max),("min",allvalues_min)]:
    for ii in range(sample_counter):
        query_name = name_desc[ii][0];
        query_desc = name_desc[ii][1];
        arr_i = allvalues[ii];
        
        funcs = [
            ("dot",dot_product,True),("cos_sim",cos_sim,True),("euc_dist",euc_dist,False),("euc_dist_norm",euc_dist_norm,False),("correl",correl,True)
        ];
        res = {};
        for ff in list(funcs):
            res[ff[0]] = [];

        num_batches = math.ceil(sample_counter/batch_size);
        arr_i_expanded = arr_i.unsqueeze(0).repeat(batch_size, 1);

        for jj in range(num_batches):
            end_index = min((jj+1)*batch_size, sample_counter);
            arr_j = allvalues[jj*batch_size:end_index];
            current_siz = arr_j.shape[0];
            for tag,func in list(funcs):
                with torch.no_grad():
                    batch_res = func(arr_i_expanded[:current_siz],arr_j).detach().cpu().tolist();
                for kkk in range(current_siz):
                    globalindex = jj*batch_size+kkk;
                    if globalindex == ii:
                        continue;
                    res[tag].append((globalindex,batch_res[kkk]));

        for tag,func,reverser in list(funcs):
            outname = os.path.join(outdir,"res_"+str(ii)+"."+stag+"."+tag+".dat");
            with open(outname,"wt") as fout:
                fout.write(">"+query_name+" "+query_desc+"\n");
                t = list(sorted(res[tag],key=lambda x:x[1],reverse=reverser));
                for tt in list(t):
                    fout.write("\t".join(name_desc[tt[0]])+"\t"+"{:.7f}".format(tt[1])+"\n");
            if use_pigz:
                subprocess.run(["pigz",outname],check=True);