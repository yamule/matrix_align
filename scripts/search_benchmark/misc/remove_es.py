import re,os,sys,subprocess,gzip;
import argparse;

def check_bool(v):
    v = v.lower();
    if v == "true" or v == "1":
        return True;
    if v == "false" or v == "0":
        return False;
    raise Exception("true or false or 1 or 0 are expected.");      

parser = argparse.ArgumentParser();
parser.add_argument("--indir",required= True) ;
parser.add_argument("--outdir",required= True) ;

args = parser.parse_args();
indir = args.indir;
outdir = args.outdir;

os.mkdir(outdir);

allfiles = os.listdir(indir);
for ff in allfiles:
    if not ff.endswith("dat.gz"):
        continue;
    inname = os.path.join(indir,ff);
    outname = os.path.join(outdir,ff);
    if os.path.exists(outname):
        raise Exception();
    lines = [];
    eflag = False;
    with gzip.open(inname,"rt") as fin:
        for ll in fin:
            if ll.startswith(">e") or ll.startswith(">fullex_e") or ll.startswith(">recov_e") :
                eflag = True;
                break;

            if ll.startswith("e") or ll.startswith("fullex_e") or ll.startswith("recov_e") :
                continue;
            lines.append(ll);
    if eflag:
        continue;
    with gzip.open(outname,"wt") as fout:
        fout.write("".join(lines));
