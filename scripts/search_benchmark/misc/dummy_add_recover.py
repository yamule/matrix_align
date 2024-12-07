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
        
# dummy_add.py で作成したファイルについて、dummy 部分を削除したものを保存する
parser = argparse.ArgumentParser();
parser.add_argument("--indir",required= True) ;
parser.add_argument("--outdir",required=True) ;


args = parser.parse_args();
indir = args.indir;
outdir = args.outdir;

os.mkdir(outdir);
allfiles = os.listdir(indir);
for ff in list(sorted(allfiles)):
    if not ff.endswith("mat.gz"):
        continue;
    infile = indir+"/"+ff;
    outfile = outdir+"/"+re.sub(r"\.mat.gz","",ff)+"_recov.mat.gz";
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

    st = 0+nlen; # ヘッダは削除したので開始は 0 行目から
    en = st + origlen;

    ptt = re.split(r"[\s]+",head[1:]); # > 記号があるので 1 飛ばす
    
    seqname = re.sub(r"fullex_","recov_",ptt[0]);
    familyname = ptt[1];

    assert not os.path.exists(outfile);

    with gzip.open(outfile,"wt") as fout:
        
        fout.write(">"+seqname+" "+familyname+"\n");
        for ll in list(lines[st:en]):
            fout.write(ll);
    