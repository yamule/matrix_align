import re,os,sys,gzip;
import numpy as np;
import argparse;
import h5py;

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
parser.add_argument("--rowwise_normalization",required=True,type=check_bool) ;

args = parser.parse_args();

targetdir = args.targetdir;
rowwise_normalization = args.rowwise_normalization;
outdir = args.outdir;
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
                continue;
            if ll.startswith(">"):
                if current is not None:
                    ret.append(current);
                    current = None;
                mat = re.search(r">([^\s]+)[\s]+([^s].+)",ll);
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
stats = None;
statsfile = os.path.join(outdir,"stats.dat");
if rowwise_normalization:
    ssum = [];
    mmean = [];
    vvar = [];
    vsiz = 0;
    headseq = None;
    valcount = 0;
    for aa in list(allfiles):
        c = load_mat(aa);
        if len(c) == 0:
            raise Exception(aa+" does not have data.");
        if vsiz == 0:
            vsiz = len(c[0]["value"][0]);
            headseq = c[0];
            for _ in range(vsiz):
                ssum.append(0);
                mmean.append(0);
                vvar.append(0);
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

vsiz = 0;
allvalues = [];
numfiles = len(allfiles);
name_desc = [];
nsample_current = numfiles;
sample_counter = 0;
for ii in range(numfiles):
    aa = allfiles[ii];
    c = load_mat(aa);
    if vsiz == 0:
        vsiz = len(c[0]["value"][0]);
        
    for cc in list(c):
        name_desc.append(
            (cc["name"],re.split(r"[\s]+",cc["desc"])[0])
        );
        sample_counter += 1;
        ssum = [0.0 for jj in range(vsiz)];
        for vv in list(cc["value"]):
            for vii in range(vsiz):
                if stats is not None:
                    if stats["std"][vii] == 0:
                        ssum[vii] += vv[vii];
                    else:
                        ssum[vii] += (vv[vii]-stats["mean"][vii])/stats["std"][vii];
                else:
                    ssum[vii] += vv[vii];
        for vii in range(vsiz):
            ssum[vii] /= len(cc["value"]);
        allvalues.append(np.array(ssum));

houtname = os.path.join(outdir,"averaged.dat");
with open(houtname,"wt") as fout:
    for ii in range(sample_counter):
        fout.write("\t".join(name_desc[ii])+"\t"+"\t".join(["{:.7f}".format(x) for x in allvalues[ii].tolist()])+"\n");

def dot_product(a,b):
    assert len(a.shape) == 1;
    assert len(b.shape) == 1;
    return (a*b).sum();

def cos_sim(a,b):
    assert len(a.shape) == 1;
    assert len(b.shape) == 1;
    anorm = np.sqrt((a*a).sum());
    bnorm = np.sqrt((b*b).sum());
    if anorm == 0 or bnorm == 0:
        return 0.0;
    return ((a/anorm)*(b/bnorm)).sum();

def euc_dist(a,b):
    assert len(a.shape) == 1;
    assert len(b.shape) == 1;
    c = a-b;
    return np.sqrt((c*c).sum());

def euc_dist_norm(a,b):
    assert len(a.shape) == 1;
    assert len(b.shape) == 1;
    anorm = np.sqrt((a*a).sum());
    bnorm = np.sqrt((b*b).sum());
    c = (a/anorm)-(b/bnorm);
    return np.sqrt((c*c).sum());

def correl(a,b):
    assert len(a.shape) == 1;
    assert len(b.shape) == 1;
    return np.corrcoef(a,b)[0,1];

all_results = [];
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
        outname = os.path.join(outdir,"res_"+str(ii)+"."+tag+".dat");
        with open(outname,"wt") as fout:
            fout.write(">"+query_name+" "+query_desc+"\n");
            t = list(sorted(res[tag],key=lambda x:x[1]));
            for tt in list(t):
                fout.write("\t".join(name_desc[tt[0]])+"\t"+"{:.7f}".format(tt[1])+"\n");
