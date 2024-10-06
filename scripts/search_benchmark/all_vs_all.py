import re,os,sys,gzip;
import numpy as np;
import argparse;
import subprocess;

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
parser.add_argument("--global_per_channel_normalization",required=True,type=check_bool) ;
parser.add_argument("--pigz",required=False,default=True,type=check_bool) ;
parser.add_argument("--unbiased_global_stats",required=False,default=False,type=check_bool) ;


args = parser.parse_args();

targetdir = args.targetdir;
global_per_channel_normalization = args.global_per_channel_normalization;
outdir = args.outdir;
use_pigz = args.pigz;
use_unbiased_global_stats = args.unbiased_global_stats;

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
                "index:\t{}\tvar:\t{:.7f}\tmean:\t{:.7f}\n".format(ii,vvar[ii],mmean[ii])
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
            (cc["name"],re.split(r"[\s]+",cc["desc"])[0])
        );
        sample_counter += 1;
        ssum = [0.0 for jj in range(vsiz)];
        mmax = [None for jj in range(vsiz)];
        mmin = [None for jj in range(vsiz)];
        for vv in list(cc["value"]):
            for vii in range(vsiz):
                
                if stats is not None:
                    if stats["std"][vii] == 0:
                        xvalue = 0.0;
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
        allvalues_average.append(np.array(ssum));
        allvalues_max.append(np.array(mmax));
        allvalues_min.append(np.array(mmin));

for (stag,val) in [("average",allvalues_average),("max",allvalues_max),("min",allvalues_min)]:
    houtname = os.path.join(outdir,stag+".dat");
    with open(houtname,"wt") as fout:
        for ii in range(sample_counter):
            fout.write("\t".join(name_desc[ii])+"\t"+"\t".join(["{:.7f}".format(x) for x in val[ii].tolist()])+"\n");


def dot_product(a,b):
    assert len(a.shape) == 1;
    assert len(b.shape) == 1;
    return (a*b).sum();

def cos_sim(a,b):
    assert len(a.shape) == 1;
    assert len(b.shape) == 1;
    anorm = np.sqrt((a*a).sum())+EPSILON;
    bnorm = np.sqrt((b*b).sum())+EPSILON;
    return ((a/anorm)*(b/bnorm)).sum();

def euc_dist(a,b):
    assert len(a.shape) == 1;
    assert len(b.shape) == 1;
    c = a-b;
    return np.sqrt((c*c).sum());

def euc_dist_norm(a,b):
    assert len(a.shape) == 1;
    assert len(b.shape) == 1;
    anorm = np.sqrt((a*a).sum())+EPSILON;
    bnorm = np.sqrt((b*b).sum())+EPSILON;
    c = (a/anorm)-(b/bnorm);
    return np.sqrt((c*c).sum());

def correl(a,b):
    assert len(a.shape) == 1;
    assert len(b.shape) == 1;

    if np.all(a == a[0]) or np.all(b == b[0]):
        return 0.0;

    return np.corrcoef(a,b)[0,1];

for (stag,allvalues) in [("average",allvalues_average),("max",allvalues_max),("min",allvalues_min)]:
    for ii in range(sample_counter):
        query_name = name_desc[ii][0];
        query_desc = name_desc[ii][1];
        arr_i = allvalues[ii];
        
        funcs = [
            ("dot",dot_product),("cos_sim",cos_sim),("euc_dist",euc_dist),("euc_dist_norm",euc_dist_norm),("correl",correl)
        ];
        res = {};
        for ff in list(funcs):
            res[ff[0]] = [];
        for jj in range(sample_counter):
            if ii == jj:
                continue;
            arr_j = allvalues[jj];
            for tag,func in list(funcs):
                res[tag].append((jj,func(arr_i,arr_j)));
        for tag,func in list(funcs):
            outname = os.path.join(outdir,"res_"+str(ii)+"."+stag+"."+tag+".dat");
            with open(outname,"wt") as fout:
                fout.write(">"+query_name+" "+query_desc+"\n");
                t = list(sorted(res[tag],key=lambda x:x[1]));
                for tt in list(t):
                    fout.write("\t".join(name_desc[tt[0]])+"\t"+"{:.7f}".format(tt[1])+"\n");
            if use_pigz:
                subprocess.run(["pigz",outname],check=True);