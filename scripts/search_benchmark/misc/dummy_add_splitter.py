import re,os,sys,gzip;
import argparse;
import random;
def check_bool(v):
    v = v.lower();
    if v == "true" or v == "1":
        return True;
    if v == "false" or v == "0":
        return False;
    raise Exception("true or false or 1 or 0 are expected.");
        
# dummy_add.py で作成したファイルを分割する
parser = argparse.ArgumentParser();
parser.add_argument("--indir",required= True) ;
parser.add_argument("--outdir",required=True) ;
parser.add_argument("--split_length", required=False, default=100,type=int) ;

random.seed(123);

max_ex_length = 600;

args = parser.parse_args();
indir = args.indir;
outdir = args.outdir;
split_length = args.split_length;

os.mkdir(outdir);
allfiles = os.listdir(indir);
for ff in list(sorted(allfiles)):
    if not ff.endswith("mat.gz"):
        continue;
    infile = indir+"/"+ff;
    outbase = outdir+"/"+re.sub(r"\.mat\.gz$","",ff);
    with gzip.open(infile,"rt") as fin:
        lines = fin.readlines();
    head = lines.pop(0);
    if lines[-1].startswith("//"):
        lines.pop();
    
    for ll in list(lines):
        assert re.search("^[A-Z]+[\s]+",ll);
    slen = len(lines); # 全長
    mat = re.search(r"N=([0-9]+)",head);
    nlen = int(mat.group(1)); # N 末に付加された配列の長さ
    mat = re.search(r"C=([0-9]+)",head);
    clen = int(mat.group(1)); # C 末に付加された配列の長さ
    
    origlen = slen - nlen - clen;

    ptt = re.split(r"[\s]+",head[1:]); # > 記号があるので 1 飛ばす
    
    seqname = ptt[0];
    familyname = ptt[1];

    if origlen < split_length*0.2:
        st = max([0,nlen-random.randint(0,int(split_length*0.5))]);
        lastpos = min([slen-clen+random.randint(0,int(split_length*0.5)),slen]);
    else:
        st = max([0,nlen-random.randint(0,int(split_length*0.4))]);
        lastpos = min([slen-clen+random.randint(0,int(split_length*0.4)),slen]);

    splitcount = 0;
    while True:
        en = st+split_length;

        if en > lastpos:
            en = lastpos;
            st = max([0,en-split_length]);

        outname = outbase+"_"+str(splitcount)+".mat.gz";
        assert not os.path.exists(outname);
        outseqname = seqname+"%"+str(splitcount);
        with gzip.open(outname,"wt") as fout:
            fout.write(">"+outseqname+" "+familyname+"\n");
            for ll in list(lines[st:en]):
                fout.write(ll);
        splitcount += 1;
        st = st+split_length//2;
        if en == lastpos:
            break;
