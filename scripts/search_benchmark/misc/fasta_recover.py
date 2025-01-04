import re,os,sys;
import gzip;

# ディレクトリに含まれる mat file から FASTA の配列を回復する

indir = sys.argv[1];
outfile = sys.argv[2];
allfiles = list(sorted(os.listdir(indir)));
fout = open(outfile,"wt");
for aa in list(allfiles):
    if not aa.endswith("mat.gz"):
        sys.stderr.write(aa+" was skipped.\n");
        continue;
    with gzip.open(indir+"/"+aa,"rt") as fin:
        for ll in fin:
            if ll.startswith(">"):
                fout.write(ll);
            else:
                fout.write(ll[0]);
        fout.write("\n");
