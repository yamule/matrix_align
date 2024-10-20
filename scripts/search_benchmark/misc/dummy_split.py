import re,os,sys;
import argparse;
import random;
def check_bool(v):
    v = v.lower();
    if v == "true" or v == "1":
        return True;
    if v == "false" or v == "0":
        return False;
    raise Exception("true or false or 1 or 0 are expected.");
        
# REPRESENTATION の POOLING が断片化された場合でも有効か検証するための、断片化後の配列を作成するためのスクリプト
# infile 内の配列の前後に split_length*0.4 未満のアミノ酸残基を付加し、split_length に断片化し、outfile に保存する
# 付加される配列は dummy_source から選ばれ、N 末側に付加される場合は C 末側を、C 末側に付加される場合は N 末側が付加される
# 断片化は split_length の長さを N 末から split_length//2 シフトしながら切り出し、最後の断片の長さが足りない場合は 最終残基インデクス-split_length+1  からとる
parser = argparse.ArgumentParser();
parser.add_argument("--infile",help='Multi-FASTA フォーマットのファイル',required= True) ;
parser.add_argument("--dummy_source",help='Multi-FASTA フォーマットのファイル',required= True) ; # 最長 split_length*0.4 まで追加される
parser.add_argument("--outfile",required=True) ;
parser.add_argument("--split_length",required=True,type=int) ;# 断片化後の配列の長さ

random.seed(123);

args = parser.parse_args();
infile = args.infile;
dummy_source = args.dummy_source;
outfile = args.outfile;
split_length = args.split_length;

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

inseq = loadFasta(infile);
dummyseq = loadFasta(dummy_source);

with open(outfile,"wt") as fout:
    for ff in list(inseq):
        basename = ff["name"];
        baseseq = ff["seq"];
        if baseseq.startswith("M"): # N 末がメチオニンである場合削る
            baseseq = baseseq[1:];
        if len(baseseq) < split_length*0.2:
            sys.stderr.write("Warning:"+basename+" is really short and longer dummy sequence will be added.\n");
        exseq = ""; # split_length より短い場合は split_length になるようにする
        while len(exseq) < split_length:
            if len(baseseq) < split_length*0.2:
                dlen_n = random.randint(0,round(split_length*0.5));
                dlen_c = random.randint(0,round(split_length*0.5));
            else:
                dlen_n = random.randint(0,round(split_length*0.4));
                dlen_c = random.randint(0,round(split_length*0.4));

            nseq = "";
            cseq = "";
            while len(nseq) < dlen_n:
                nseq = dummyseq[random.randrange(0,len(dummyseq))]["seq"];
            while len(cseq) < dlen_c:
                cseq = dummyseq[random.randrange(0,len(dummyseq))]["seq"];
            
            
            nseq = nseq[:-1]; # 開始メチオニンのように末端を覚えていると良くないので末端残基は削る
            cseq = cseq[1:];

            if len(nseq) > dlen_n:
                nseq = nseq[-dlen_n:];
            if len(cseq) > dlen_c:
                cseq = cseq[:dlen_c];
            
            exseq = nseq+baseseq+cseq
        

        seqlen = len(exseq);
        st = 0;
        fragcount = 0;
        while True:        
            en = st+split_length;
            if en > seqlen:
                en = seqlen;
                st = max([0,en-split_length]);

            fout.write(">ex_"+basename+"%"+str(fragcount)+" "+ff["desc"]+"\n");
            fout.write(exseq[st:en]+"\n");
            fragcount += 1;
            st = st+split_length//2;
            if en == seqlen:
                break;