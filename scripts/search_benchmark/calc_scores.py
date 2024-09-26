import re,gzip,os,sys;
targetdir = sys.argv[1];

is_cath = False;
if "cath" in targetdir:
    is_cath = True;

def get_groupids(gcode):
    groups = re.split(r"\.",gcode);
    assert len(groups) > 3,"??? invalid format. "+gcode+"\n";
    g1 = ".".join(groups);
    g2 = ".".join(groups[:-1]);
    g3 = ".".join(groups[:-2]);
    return (g1,g2,g3);

targetdir 
allfiles = {};
group_1 = {}; 
group_2 = {}; 
group_3 = {};
file_to_group = {};
for aa in list(sorted(os.listdir(targetdir))):

    if not aa.endswith("dat.gz"):
        continue;
    ptt = re.split(r"\.",aa); # <スコアの種類>.dat.gz というファイル名規則になっている想定
    stype = ptt[-3];

    if len(file_to_group) > 1000:
        break;
    if stype not in allfiles:
        allfiles[stype] = [];
    filename = os.path.join(targetdir,aa);
    allfiles[stype].append(filename);
    
    with gzip.open(filename,"rt") as fin:
        ll = fin.readline();
        ptt = re.split(r"[\s]+",ll[1:]);
        name = ptt[0];
        g1,g2,g3 = get_groupids(ptt[1]);
        file_to_group[filename] = g1;
        if g1 not in group_1:
            group_1[g1] = set();
        if g2 not in group_2:
            group_2[g2] = set();
        if g3 not in group_3:
            group_3[g3] = set();

        group_1[g1].add(name);
        group_2[g2].add(name);
        group_3[g3].add(name);

for score_type in list(sorted(allfiles.keys())):
    g1_score_fp1 = 0;
    g2_score_fp1 = 0;
    g3_score_fp1 = 0;

    g1_score_fp10 = 0;
    g2_score_fp10 = 0;
    g3_score_fp10 = 0;

    g1_count = 0;
    g2_count = 0;
    g3_count = 0;

    g1_num = 0;
    g2_num = 0;
    g3_num = 0;

    for kk,vv in list(group_1.items()):
        if len(vv) > 1:
            g1_num += 1;
    for kk,vv in list(group_2.items()):
        if len(vv) > 1:
            g2_num += 1;
    for kk,vv in list(group_3.items()):
        if len(vv) > 1:
            g3_num += 1;

    higher_is_better = "higher_is_better";
    for ff in list(allfiles[score_type]):

        g1,g2,g3 = get_groupids(file_to_group[ff]);

        if (is_cath and len(group_2[g2]) == 1) or len(group_3[g3]) == 1:
            continue;
        with gzip.open(ff,"rt") as fin:
            alllines = fin.readlines();
            ptt = re.split(r"[\s]+",alllines[0][1:]);
            name = ptt[0];
            assert file_to_group[ff] == ptt[1];

            g1_fp1 = -1;
            g1_fp10 = -1;
            g2_fp1 = -1;
            g2_fp10 = -1;
            g3_fp1 = -1;
            g3_fp10 = -1;
            
            tpcount_g1 = 0;
            tpcount_g2 = 0;
            tpcount_g3 = 0;

            fpcount_g1 = 0;
            fpcount_g2 = 0;
            fpcount_g3 = 0;
            targetscores_ = [];
            for pii in range(1,len(alllines)):
                ptt = re.split(r"[\s]+",alllines[pii]); # name[\t]groupid[\t]score となっている想定。
                targetscores_.append((ptt[0],ptt[1],float(ptt[2])));
            if "euc" in score_type:
                targetscores = list(sorted(targetscores_,key=lambda x: x[-1]));
                higher_is_better = "lower_is_better";
            else:
                targetscores = list(sorted(targetscores_,key=lambda x: x[-1],reverse = True));
            del targetscores_;

            for ptt in list(targetscores):
                if ptt[0] == name:
                    continue;
                ############# 上位グループが同じでも FP とみなす。
                bg1,bg2,bg3 = get_groupids(ptt[1]);
                
                if bg1 == g1:
                    tpcount_g1 += 1;
                else:
                    fpcount_g1 += 1;
                    if fpcount_g1 == 1:
                        g1_fp1 = tpcount_g1;
                    if fpcount_g1 == 10:
                        g1_fp10 = tpcount_g1;

                if bg2 == g2:
                    tpcount_g2 += 1;
                else:
                    fpcount_g2 += 1;
                    if fpcount_g2 == 1:
                        g2_fp1 = tpcount_g2;
                    if fpcount_g2 == 10:
                        g2_fp10 = tpcount_g2;
        
                if bg3 == g3:
                    tpcount_g3 += 1;
                else:
                    fpcount_g3 += 1;
                    if fpcount_g3 == 1:
                        g3_fp1 = tpcount_g3;
                    if fpcount_g3 == 10:
                        g3_fp10 = tpcount_g3;

            #assert len(group_1[g1]) == tpcount_g1+1,g1+"\t"+str(len(group_1[g1]))+" vs "+str(tpcount_g1+1);
            #assert len(group_2[g2]) == tpcount_g2+1;
            #assert len(group_3[g3]) == tpcount_g3+1;
            if tpcount_g1 > 0:
                g1_score_fp1 += g1_fp1/float(tpcount_g1)/float(len(group_1[g1]))/g1_num;
                g1_score_fp10 += g1_fp10/float(tpcount_g1)/float(len(group_1[g1]))/g1_num;
            if tpcount_g2 > 0:
                g2_score_fp1 += g2_fp1/float(tpcount_g2)/float(len(group_2[g2]))/g2_num;
                g2_score_fp10 += g2_fp10/float(tpcount_g2)/float(len(group_2[g2]))/g2_num;
            if tpcount_g3 > 0:
                g3_score_fp1 += g3_fp1/float(tpcount_g3)/float(len(group_3[g3]))/g3_num;
                g3_score_fp10 += g3_fp10/float(tpcount_g3)/float(len(group_3[g3]))/g3_num;
    print("=========");
    def strline(*argg):
        return "\t".join([str(x) for x in argg]);


    print(strline(targetdir,score_type,higher_is_better,"g1","fp1","count:",len(allfiles[score_type]),"score:",g1_score_fp1));
    print(strline(targetdir,score_type,higher_is_better,"g1","fp10","count:",len(allfiles[score_type]),"score:",g1_score_fp10));
    print(strline(targetdir,score_type,higher_is_better,"g2","fp1","count:",len(allfiles[score_type]),"score:",g2_score_fp1));
    print(strline(targetdir,score_type,higher_is_better,"g2","fp10","count:",len(allfiles[score_type]),"score:",g2_score_fp10));
    print(strline(targetdir,score_type,higher_is_better,"g3","fp1","count:",len(allfiles[score_type]),"score:",g3_score_fp1));
    print(strline(targetdir,score_type,higher_is_better,"g3","fp10","count:",len(allfiles[score_type]),"score:",g3_score_fp10));



"""
一個しかないものについては削除する
残ったものについてスコアごとに読み込む
ソートする
全 TP を数える
FP 一個取るまで、10 個取るまでの TP 割合を計算する
同一 Family, など属性が同じ物の数で平均する
全体で平均する
"""
