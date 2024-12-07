import os,sys,re;
import misc.smithwaterman;
import gzip;
import subprocess;

sw = misc.smithwaterman.SmithWaterman();

# PLM の比較では Family メンバより Superfamily メンバのほうがスコアが高くなっている例について
# ALIGNMENT の IDENTITY で調べた場合は FAMILY のほうが高いか確認するためのスクリプト

import argparse;

tmalign = "/home/ubuntu8/apps/TMalign/TMalign"
pdbstyle_dir = "/home/ubuntu8/data/data0/GMalign/matrix_align/nogit/large_data/id20/pdbstyle/";

def check_bool(v):
    v = v.lower();
    if v == "true" or v == "1":
        return True;
    if v == "false" or v == "0":
        return False;
    raise Exception("true or false or 1 or 0 are expected.");
def check_identity(res):
    stpos = -1;
    enpos = -1;
    for ii in range(len(res[0])):
        if stpos < 0:
            if res[0][ii] != "-" and res[1][ii] != "-":
                stpos = ii;
        if enpos < 0:
            ppos = len(res[0]) -1 -ii;
            if res[0][ppos] != "-" and res[1][ppos] != "-":
                enpos = ppos;

    shortlen = min([
        len(re.sub(r"[^A-Z]","",res[0].upper())),
        len(re.sub(r"[^A-Z]","",res[1].upper()))
    ])

    matchcount = 0;
    for ii in range(stpos,enpos+1):
        if res[0][ii] == res[1][ii]:
            matchcount += 1;

    return matchcount/max(
        [shortlen,enpos-stpos+1]
    );

def check_identity_tmalign(file1,file2):
    if not os.path.exists(file1) or not os.path.exists(file2):
        return -1;
    proc = subprocess.run([tmalign,file1,file2],stdout=subprocess.PIPE,encoding="utf-8");
    li = re.split(r"[\r\n]+",proc.stdout);
    spos = None;
    for ii in range(len(li)):
        if "denotes residue pairs of d" in li[ii]:
            spos = (ii+1,ii+3);
    assert spos is not None;
    return check_identity([li[spos[0]],li[spos[1]]])

def check_identity_sw(seq1,seq2):
    # SW のアラインメント長いか、短い方の配列の長さか
    # どちらか長い方で割る。
    # （参考：BLAST はアラインメント長。）
    res = sw.align(seq1["seq"],seq2["seq"]);
    return check_identity(res);

parser = argparse.ArgumentParser();
parser.add_argument("--fasta_file",help='全長の AA が入っている FASTA',required= True) ;
parser.add_argument("--target_dir",help='all vs all の結果が入っているディレクトリ',required= True) ;
parser.add_argument("--out_file",help='',required= True) ;

args = parser.parse_args();

fasta_file = args.fasta_file;
target_dir = args.target_dir;
out_file = args.out_file;



fass_ = sw.loadFasta(fasta_file);
name_to_fasta = {};
for ff in list(fass_):
    assert ff["name"] not in name_to_fasta;
    name_to_fasta[ff["name"]] = ff;

allfiles = os.listdir(target_dir);
for aa in list(allfiles):
    if not aa.endswith("gz"):
        continue;
    lower_is_better = "euc_dist" in aa;# euclidean distance は小さい値のほうが良い
    if not lower_is_better:
        assert "euc" not in aa; # ファイル名変えていないかチェック
    with gzip.open(target_dir+"/"+aa,"rt") as fin:
        lines = fin.readlines();
        head = lines.pop(0);
        mat = re.search(r">([^\s]+)[\s]+([^\s]+)",head);
        if mat:
            query = mat.group(1);
            qfamily = mat.group(2);
        else:
            raise Exception();
        fptt = re.split(r"\.",qfamily);
        qsuperfamily = ".".join(fptt[0:3]);

        hits = [];
        for ll in list(lines):
            ptt = re.split(r"[\s]+",re.sub(r"[\s]+$","",ll));
            hits.append(
                (ptt[0],ptt[1],float(ptt[-1])) # 名前、ファミリーID、スコア が入っているはず
            );
        
        if lower_is_better:
            hits_sorted = list(sorted(hits,key=lambda x:x[-1]));
        else:
            hits_sorted = list(reversed(sorted(hits,key=lambda x:x[-1])));

        fhit = None;
        supfhit = None;
        for hh in list(hits_sorted):
            if hh[0] == query:
                continue;
            if hh[1] == qfamily:
                if fhit is None:
                    fhit = hh;
                continue;
            fptt = re.split(r"\.",hh[1]);
            hsuperfamily = ".".join(fptt[0:3]);
            if hsuperfamily == qsuperfamily:
                if supfhit is None:
                    supfhit = hh;
            if fhit is not None and supfhit is not None:
                break;
        if fhit is not None and supfhit is not None:
            checkercode = "OK";
            if (not lower_is_better and fhit[-1] < supfhit[-1]) or (lower_is_better and fhit[-1] > supfhit[-1]):
                checkercode = "NG"
            
            fid = check_identity_sw(
                name_to_fasta[query],name_to_fasta[fhit[0]]
            );
            supfid = check_identity_sw(
                name_to_fasta[query],name_to_fasta[supfhit[0]]
            );

            fid_tm = check_identity_tmalign(
                pdbstyle_dir+"/"+query+".ent",pdbstyle_dir+"/"+fhit[0]+".ent"
            );

            supfid_tm = check_identity_tmalign(
                pdbstyle_dir+"/"+query+".ent",pdbstyle_dir+"/"+supfhit[0]+".ent"
            );

            li = "\t".join(str(x) for x in [checkercode,aa,query,qfamily,fhit[0],fhit[1],fhit[-1]
            ,fid,fid_tm,supfhit[0],supfhit[1],supfhit[-1],supfid,supfid_tm]);
            print(li);

