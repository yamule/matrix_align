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
        
# REPRESENTATION の POOLING が断片化された場合でも有効か検証するためのスクリプトだが、このスクリプト内では断片化は行わない。
# infile 内の配列の前後に max 600 アミノ酸残基の完全なタンパク質を付加し、outfile に保存する。
# desc に N=X C=X として付加された断片の長さを加える。
# 付加される配列は dummy_source から選ばれ、N 末側に付加される場合は C 末側を、C 末側に付加される場合は N 末側が付加される。
parser = argparse.ArgumentParser();
parser.add_argument("--infile",help='Multi-FASTA フォーマットのファイル',required= True) ;
parser.add_argument("--dummy_source",help='Multi-FASTA フォーマットのファイル',required= True) ; # 最長 split_length*0.4 まで追加される
parser.add_argument("--outfile",required=True) ;

random.seed(123);

max_ex_length = 600;

args = parser.parse_args();
infile = args.infile;
dummy_source = args.dummy_source;
outfile = args.outfile;

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
        methionine_flag = False;
        if baseseq.startswith("M"): # N 末がメチオニンである場合削る
            baseseq = baseseq[1:];
            methionine_flag = True;
        exseq = ""; # split_length より短い場合は split_length になるようにする
        

        nseq = "";
        cseq = "";
        counter = 0;

        cseq = dummyseq[random.randrange(0,len(dummyseq))]["seq"];
        nseq = dummyseq[random.randrange(0,len(dummyseq))]["seq"];
        while len(nseq) > max_ex_length:
            nseq = dummyseq[random.randrange(0,len(dummyseq))]["seq"];
            counter += 1;
            if counter > 10000:
                nseq = "";
                break;
        counter = 0;
        while len(cseq) > max_ex_length:
            cseq = dummyseq[random.randrange(0,len(dummyseq))]["seq"];
            counter += 1;
            if counter > 10000:
                cseq = "";
                break;
        
        if len(nseq) > 0:
            nseq = nseq[:-1]; # 開始メチオニンのように末端を覚えていると良くないので末端残基は削る
        if len(cseq) > 0:
            cseq = cseq[1:];

        exseq = nseq+baseseq+cseq
        if len(nseq) == 0 and methionine_flag:
            exseq = "M"+exseq;

        fout.write(">fullex_"+basename+" "+ff["desc"]+" N={} C={}".format(len(nseq),len(cseq))+"\n");
        fout.write(exseq+"\n");
