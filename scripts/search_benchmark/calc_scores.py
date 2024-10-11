import re,gzip,os,sys;
import copy;
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
group_1_orig = {}; 
group_2_orig = {}; 
group_3_orig = {};
file_to_group = {};
for aa in list(sorted(os.listdir(targetdir))):

    if not aa.endswith("dat.gz"):
        continue;
    ptt = re.split(r"\.",aa); # .<max, average, min>.<スコアの種類>.dat.gz というファイル名規則になっている想定
    assert len(ptt) == 5;
    stype = ptt[-4]+"."+ptt[-3];
    
    #print(aa,stype);
    #raise Exception();

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
        if g1 not in group_1_orig:
            group_1_orig[g1] = 0;
        if g2 not in group_2_orig:
            group_2_orig[g2] = 0;
        if g3 not in group_3_orig:
            group_3_orig[g3] = 0;
    # if len(allfiles[stype]) > 1000: # デバッグ用
    #    break;

for score_type in list(sorted(allfiles.keys())):
    
    g1_count = copy.deepcopy(group_1_orig);
    g2_count = copy.deepcopy(group_2_orig);
    g3_count = copy.deepcopy(group_3_orig);

    g1_score_fp1 = copy.deepcopy(group_1_orig);
    g2_score_fp1 = copy.deepcopy(group_2_orig);
    g3_score_fp1 = copy.deepcopy(group_3_orig);

    g1_score_fp10 = copy.deepcopy(group_1_orig);
    g2_score_fp10 = copy.deepcopy(group_2_orig);
    g3_score_fp10 = copy.deepcopy(group_3_orig);


    g1_score_hit1 = copy.deepcopy(group_1_orig);
    g2_score_hit1 = copy.deepcopy(group_2_orig);
    g3_score_hit1 = copy.deepcopy(group_3_orig);

    g1_score_hit10 = copy.deepcopy(group_1_orig);
    g2_score_hit10 = copy.deepcopy(group_2_orig);
    g3_score_hit10 = copy.deepcopy(group_3_orig);

    g1_score_roc = copy.deepcopy(group_1_orig);
    g2_score_roc = copy.deepcopy(group_2_orig);
    g3_score_roc = copy.deepcopy(group_3_orig);

    g1_score_pr_roc = copy.deepcopy(group_1_orig);
    g2_score_pr_roc = copy.deepcopy(group_2_orig);
    g3_score_pr_roc = copy.deepcopy(group_3_orig);


    higher_is_better = "higher_is_better";
    for ff in list(allfiles[score_type]):

        g1,g2,g3 = get_groupids(file_to_group[ff]);

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
                    continue;
                else:
                    fpcount_g1 += 1;
                    if fpcount_g1 == 1:
                        g1_fp1 = tpcount_g1;
                    if fpcount_g1 == 10:
                        g1_fp10 = tpcount_g1;

                if bg2 == g2:
                    tpcount_g2 += 1;
                    continue;
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


            # 全 TP, FP を数えたので ROC-AUC を計算する
            roc_g1 = 0;
            roc_g2 = 0;
            roc_g3 = 0;

            tp_current_g1 = 0;
            fp_current_g1 = 0;
            tp_current_g2 = 0;
            fp_current_g2 = 0;
            tp_current_g3 = 0;
            fp_current_g3 = 0;

            tp_prev_g1 = 0;
            fp_prev_g1 = 0;
            tp_prev_g2 = 0;
            fp_prev_g2 = 0;
            tp_prev_g3 = 0;
            fp_prev_g3 = 0;

            for ptt in list(targetscores):
                if ptt[0] == name:
                    continue;
                ############# 上位グループが同じでも FP とみなす。
                bg1,bg2,bg3 = get_groupids(ptt[1]);

                if bg1 == g1:
                    tp_current_g1 += 1;
                    continue;
                else:
                    fp_current_g1 += 1;
                    roc_g1 += (tp_prev_g1+tp_current_g1)/float(fpcount_g1)/tpcount_g1/2.0;
                    fp_prev_g1 = fp_current_g1;
                    tp_prev_g1 = tp_current_g1;
        
                if bg2 == g2:
                    tp_current_g2 += 1;
                    continue;
                else:
                    fp_current_g2 += 1;
                    roc_g2 += (tp_prev_g2+tp_current_g2)/float(fpcount_g2)/tpcount_g2/2.0;
                    fp_prev_g2 = fp_current_g2;
                    tp_prev_g2 = tp_current_g2;
        
                if bg3 == g3:
                    tp_current_g3 += 1;
                else:
                    fp_current_g3 += 1;
                    roc_g3 += (tp_prev_g3+tp_current_g3)/float(fpcount_g3)/tpcount_g3/2.0;
                    fp_prev_g3 = fp_current_g3;
                    tp_prev_g3 = tp_current_g3;

            # pr_roc を計算する
            pr_roc_g1 = 0;
            pr_roc_g2 = 0;
            pr_roc_g3 = 0;

            tp_current_g1 = 0;
            fp_current_g1 = 0;
            tp_current_g2 = 0;
            fp_current_g2 = 0;
            tp_current_g3 = 0;
            fp_current_g3 = 0;

            tp_prev_g1 = 0;
            fp_prev_g1 = 0;
            tp_prev_g2 = 0;
            fp_prev_g2 = 0;
            tp_prev_g3 = 0;
            fp_prev_g3 = 0;

            for ptt in list(targetscores):
                if ptt[0] == name:
                    continue;
                ############# 上位グループが同じでも FP とみなす。
                bg1,bg2,bg3 = get_groupids(ptt[1]);

                if bg1 == g1:
                    tp_current_g1 += 1;
                    pprev = tp_prev_g1/(tp_prev_g1+fp_prev_g1);
                    pcurrent = tp_current_g1/float(tp_current_g1+fp_current_g1);
                    rp_roc_g1 += (pprev+pcurrent)/float(tpcount_g1)/2.0;
                    tp_prev_g1 = tp_current_g1;
                    fp_prev_g1 = fp_current_g1;
                    continue;
                else:
                    fp_current_g1 += 1;
                    
        
                if bg2 == g2:
                    tp_current_g2 += 1;
                    pprev = tp_prev_g2/(tp_prev_g2+fp_prev_g2);
                    pcurrent = tp_current_g2/float(tp_current_g2+fp_current_g2);
                    rp_roc_g2 += (pprev+pcurrent)/float(tpcount_g2)/2.0;
                    tp_prev_g2 = tp_current_g2;
                    fp_prev_g2 = fp_current_g2;
                    continue;
                else:
                    fp_current_g2 += 1;
        
                if bg3 == g3:
                    tp_current_g3 += 1;
                    pprev = tp_prev_g3/(tp_prev_g3+fp_prev_g3);
                    pcurrent = tp_current_g3/(tp_current_g3+fp_current_g3);
                    rp_roc_g3 += (pprev+pcurrent)/float(tpcount_g3)/2.0;
                    tp_prev_g3 = tp_current_g3;
                    fp_prev_g3 = fp_current_g3;
                else:
                    fp_current_g3 += 1;
                    
            if tpcount_g1 > 0:
                g1_count[g1] += 1;
                g1_score_fp1[g1] += g1_fp1/float(tpcount_g1);
                g1_score_fp10[g1] += g1_fp10/float(tpcount_g1);

                g1_score_roc[g1] += roc_g1;
                g1_score_pr_roc[g1] += pr_roc_g1;
                
                if g1_fp1 > 0:
                    g1_score_hit1[g1] += 1.0;
                if g1_fp10 > 0:
                    g1_score_hit10[g1] += 1.0;
                    
            if tpcount_g2 > 0:
                g2_count[g2] += 1;
                g2_score_fp1[g2] += g2_fp1/float(tpcount_g2);
                g2_score_fp10[g2] += g2_fp10/float(tpcount_g2);
                
                g2_score_roc[g2] += roc_g2;
                g2_score_pr_roc[g2] += pr_roc_g2;
                
                if g2_fp1 > 0:
                    g2_score_hit1[g2] += 1.0;
                if g2_fp10 > 0:
                    g2_score_hit10[g2] += 1.0;
                    

            if tpcount_g3 > 0:
                g3_count[g3] += 1;
                g3_score_fp1[g3] += g3_fp1/float(tpcount_g3);
                g3_score_fp10[g3] += g3_fp10/float(tpcount_g3);

                g3_score_roc[g3] += roc_g3;
                g3_score_pr_roc[g3] += pr_roc_g3;
                
                if g3_fp1 > 0:
                    g3_score_hit1[g3] += 1.0;
                if g3_fp10 > 0:
                    g3_score_hit10[g3] += 1.0;
            
    print("=========");
    num_queries  = {};
    score_sum_fp1 = {};
    score_sum_fp10 = {};

    # Top1 or 10 に TP が含まれていたかどうか https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae119/7735315
    score_sum_hit1 = {}; 
    score_sum_hit10 = {}; 

    score_sum_roc = {}; 
    score_sum_pr_roc = {}; 

    for (tag,count,score_fp1,score_fp10,score_hit1,score_hit10,score_roc,score_pr_roc) in [
        ("g1",g1_count,g1_score_fp1,g1_score_fp10,g1_score_hit1,g1_score_hit10,g1_score_roc,g1_score_pr_roc),
        ("g2",g2_count,g2_score_fp1,g2_score_fp10,g2_score_hit1,g2_score_hit10,g2_score_roc,g2_score_pr_roc),
        ("g3",g3_count,g3_score_fp1,g3_score_fp10,g3_score_hit1,g3_score_hit10,g3_score_roc,g3_score_pr_roc),
    ]:
        counted_class = 0;
        num_queries[tag]  = 0;
        for gg in list(count.keys()):
            if count[gg] > 0:# 同一クラスに存在するメンバ数で割って平均化する
                counted_class += 1; 
                num_queries[tag] += count[gg];
                score_fp1[gg] /= count[gg];
                score_fp10[gg] /= count[gg];
                score_hit1[gg] /= count[gg];
                score_hit10[gg] /= count[gg];
                score_roc[gg] /= count[gg];
                score_pr_roc[gg] /= count[gg];
                
        score_sum_fp1[tag] = 0;
        score_sum_fp10[tag] = 0;
        score_sum_hit1[tag] = 0;
        score_sum_hit10[tag] = 0;
        score_sum_roc[tag] = 0;
        score_sum_pr_roc[tag] = 0;
        for gg in list(count.keys()):
            if count[gg] > 0: # 全クラス数で割って更に平均化する
                score_fp1[gg] /= counted_class;
                score_fp10[gg] /= counted_class;
                score_hit1[gg] /= counted_class;
                score_hit10[gg] /= counted_class;
                score_roc[gg] /= counted_class;
                score_pr_roc[gg] /= counted_class;
                
                score_sum_fp1[tag] += score_fp1[gg];
                score_sum_fp10[tag] += score_fp10[gg];
                score_sum_hit1[tag] += score_hit1[gg];
                score_sum_hit10[tag] += score_hit10[gg];
                score_sum_roc[tag] += score_hit1[gg];
                score_sum_pr_roc[tag] += score_hit10[gg];


    def strline(*argg):
        return "\t".join([str(x) for x in argg]);


    print(strline(targetdir,score_type,higher_is_better,"g1","untilfp1","count:",num_queries["g1"],"score:",score_sum_fp1["g1"]));
    print(strline(targetdir,score_type,higher_is_better,"g2","untilfp1","count:",num_queries["g2"],"score:",score_sum_fp1["g2"]));
    print(strline(targetdir,score_type,higher_is_better,"g3","untilfp1","count:",num_queries["g3"],"score:",score_sum_fp1["g3"]));
    
    print(strline(targetdir,score_type,higher_is_better,"g1","hit_at_1","count:",num_queries["g1"],"score:",score_sum_hit1["g1"]));
    print(strline(targetdir,score_type,higher_is_better,"g1","hit_at_10","count:",num_queries["g1"],"score:",score_sum_hit10["g1"]));
    print(strline(targetdir,score_type,higher_is_better,"g2","hit_at_1","count:",num_queries["g2"],"score:",score_sum_hit1["g2"]));
    print(strline(targetdir,score_type,higher_is_better,"g2","hit_at_10","count:",num_queries["g2"],"score:",score_sum_hit10["g2"]));
    print(strline(targetdir,score_type,higher_is_better,"g3","hit_at_1","count:",num_queries["g3"],"score:",score_sum_hit1["g3"]));
    print(strline(targetdir,score_type,higher_is_better,"g3","hit_at_10","count:",num_queries["g3"],"score:",score_sum_hit10["g3"]));

    print(strline(targetdir,score_type,higher_is_better,"g1","roc","count:",num_queries["g1"],"score:",score_sum_roc["g1"]));
    print(strline(targetdir,score_type,higher_is_better,"g2","roc","count:",num_queries["g2"],"score:",score_sum_roc["g2"]));
    print(strline(targetdir,score_type,higher_is_better,"g3","roc","count:",num_queries["g3"],"score:",score_sum_roc["g3"]));
    
    print(strline(targetdir,score_type,higher_is_better,"g1","pr_roc","count:",num_queries["g1"],"score:",score_sum_pr_roc["g1"]));
    print(strline(targetdir,score_type,higher_is_better,"g2","pr_roc","count:",num_queries["g2"],"score:",score_sum_pr_roc["g2"]));
    print(strline(targetdir,score_type,higher_is_better,"g3","pr_roc","count:",num_queries["g3"],"score:",score_sum_pr_roc["g3"]));
    
    

"""
一個しかないものについては削除する
残ったものについてスコアごとに読み込む
ソートする
全 TP を数える
FP 一個取るまで、10 個取るまでの TP 割合を計算する
同一 Family, など属性が同じ物の数で平均する
全体で平均する
"""
