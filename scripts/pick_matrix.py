import re,os,sys,gzip;
import argparse;

def check_bool(v):
    v = v.lower();
    if v == "true" or v == "1":
        return True;
    if v == "false" or v == "0":
        return False;
    raise Exception("true or false or 1 or 0 are expected.");
        

parser = argparse.ArgumentParser();
parser.add_argument("--infile",help='Matrix 出力ファイル',required= True) ;
parser.add_argument("--targets",help='抽出したい配列の名前',required= True) ;
args = parser.parse_args();

if args.infile.endswith(".gz"):
    fin = gzip.open(args.infile,"rt");
else:
    fin = open(args.infile,"rt");

targett = set();
for pp in re.split(",",args.targets):
    pp = re.sub(r"[\s]","",pp);
    if len(pp) > 0:
        targett.add(">"+pp);
flag = False;
for ll in fin:
    if ll.startswith(">"):
        pt = re.split(r"[\s]+",ll);
        if pt[0] in targett:
            flag = True;
        else:
            flag = False;
    if flag:
        print(ll,end="");

fin.close();