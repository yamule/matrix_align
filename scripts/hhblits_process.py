import argparse;
import os,re,sys;
import subprocess;

parser = argparse.ArgumentParser()
parser.add_argument('--fasta_file')
parser.add_argument('--outdir')
parser.add_argument('--hhblits_path')
parser.add_argument('--hhdb_path')
parser.add_argument('--round',required=False,type=int,default=5)

args = parser.parse_args()
outdir = args.outdir;
fasta_file = args.fasta_file;
hhblits_path = args.hhblits_path;
hhdb_path = args.hhdb_path;
float_round = args.round;

hhoptions = [
    "-cpu","12",
    "-d",hhdb_path,
    "-n","3"
];

def loadFasta(filename):
    if filename.endswith("gz"):
        fin = gzip.open(filename,"rt");
    else:
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
            cdict["seq"] += re.sub("[\s]","",ll).upper();
            
    if(len(ret[0]["seq"]) == 0):
        ret.pop(0);
    fin.close();
    return ret;

if not os.path.exists(outdir):
    os.mkdir(outdir);

fastas = loadFasta(fasta_file);
for ff in fastas:
    pname = re.sub(r"[^\-.A-Za-z0-9]","_",ff["name"]);
    fasout = outdir+"/"+pname+".fas";
    hhmout = outdir+"/"+pname+".hhm";
    matout = outdir+"/"+pname+".mat";

    if os.path.exists(matout):
        sys.stderr.write(matout+" already exists. skipped.\n");
        sys.stderr.flush();
        continue;


    with open(fasout,"wt") as fout:
        fout.write(">"+ff["name"]+" "+ff["desc"]+"\n");
        fout.write(ff["seq"]+"\n");
    resout = outdir+"/"+pname+".hhr";
    a3mout = outdir+"/"+pname+".a3m";
    proc = subprocess.run(
        [hhblits_path,"-i",fasout,"-o",resout,"-oa3m",a3mout,"-ohhm",hhmout,*hhoptions],stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8"
    );
    os.system("gzip "+a3mout);
    sys.stderr.write(proc.stderr);
    sys.stdout.flush();
    sys.stderr.flush();

    matlines = [];
    with open(hhmout,"rt") as fin:
        lines = fin.readlines();
        sline = 0;
        for ii in range(len(lines)):
            if lines[ii].startswith("NULL"):
                sline = ii+1;
                if lines[ii+1].startswith("HMM"):
                    sline = ii+4;
                break;
        for ii in range(sline,len(lines)):
            if lines[ii].startswith("//"):
                if ii != len(lines)-1:
                    sys.stderr.write("Does "+hhmout+" have multiple hmm? multiple hmm is not supported.\n");
                break;
            if (ii-sline)%3 == 0:
                ptt = re.split(r"[\s]+",lines[ii]);
                pres = [ptt[0]];
                for jj in range(2,22):
                    if ptt[jj] == "*":
                        pres.append("0");
                    else:
                        pres.append(
                            str(
                                round(2.0**(float(ptt[jj])*-1/1000.0),float_round)
                            )
                        );
                matlines.append("\t".join(pres));
    with open(matout,"wt") as fout:
        fout.write(">"+ff["name"]+" "+ff["desc"]+"\n");
        fout.write("\n".join(matlines));
        fout.write("\n");
        