import re,sys,os,gzip;


# XXXXX と X が連続で 5 つある配列を除外する。
# JUZBOX はすべて X に変更する。
# 同一アミノ酸配列がある場合はまとめ、標準出力に出力する。
# SCOPe を処理する場合、同一アミノ酸配列が別ファミリーに所属している場合除外する。

def loadFasta(filename):
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
            cdict["seq"] += re.sub("[\s]","",ll);
            
    if(len(ret[0]["seq"]) == 0):
        ret.pop(0);
    fin.close();
    return ret;
infile = sys.argv[1];
outfile = sys.argv[2];
fass = loadFasta(infile);

printed = {};
familyname = {};
sameseq = {};
errorseq = {};
for ff in fass:

    ff["seq"] = re.sub(r"[JUZBOX]","X",ff["seq"].upper());

    if "scope" in infile:
        mmat = re.search(r"^[\s]*([^\s]+)",ff["desc"]);
        family = mmat.group(1);
    if "pdb_seqres" in infile:
        if "mol:protein" not in ff["desc"]:
            continue;
    if "XXXXX" in ff["seq"]:
        continue;
    if ff["seq"] in printed:
        if "scope" in infile:
            if family != familyname[ff["seq"]]:
                errorseq[ff["seq"]] = 100;
                sys.stderr.write("Different family id "+family+" vs "+familyname[ff["seq"]]+", "+printed[ff["seq"]]["name"]+" vs "+ff["name"]+", "+ff["seq"]+"\n");
        sameseq[ff["seq"]].append(ff["name"]);
        continue;
    else:
        printed[ff["seq"]] = ff;
        sameseq[ff["seq"]] = [];
    sameseq[ff["seq"]].append(ff["name"]);
    if "scope" in infile:
        mmat = re.search(r"^[\s]*([^\s]+)",ff["desc"]);
        familyname[ff["seq"]] = family;
        
with open(outfile,"wt") as fout:
    for kk in list(sorted(printed.keys())):
        ff = printed[kk];
        if ff["seq"] in errorseq:
            continue;
        fout.write(">"+ff["name"]+" "+ff["desc"]+"\n");
        fout.write(ff["seq"]+"\n");
        print(">"+ff["name"]);
        print(",".join(sameseq[ff["seq"]]));
    

