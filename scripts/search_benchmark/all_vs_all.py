import re,os,sys,gzip;
import numpy as np;
import argparse;

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
parser.add_argument("--rowwise_normalization",required=True) ;
parser.add_argument("--function_type",help="correl, dist, dot, cos",required=True) ;

args = parser.parse_args();

targetdir = args.targetdir;
function_type = args.function_type;
rowwise_normalization = args.rowwise_normalization;
outdir = args.outdir;
if os.path.exists(outdir):
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
        mmean.append(ssum[ii]/float(valcount));
    
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
                "index:\t{}\tvar:\t{}\tmean:\t{}\n".format(ii,vvar[ii],mmean[ii])
            );
