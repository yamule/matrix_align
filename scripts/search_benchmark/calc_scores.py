import re,gzip,os,sys;
import copy;
import math;

targetdir = sys.argv[1];

is_cath = False;
if "cath" in targetdir:
    is_cath = True;

# {score:score,label:label}
def calc_roc_auc(targets,reverse):
    ssorted = list(sorted(targets,key=lambda x:x["score"],reverse=reverse));
    prev_score = ssorted[0]["score"];
    samplenum = len(ssorted);
    positives = 0;
    for ss in list(ssorted):
        if ss["label"] == 1:
            positives += 1;
    negatives = samplenum - positives;
    tpcount = 0;
    fpcount = 0;
    ii = 0;
    points = [];
    points.append((0.0,0.0));
    if positives == 0 or negatives == 0:
        return float("nan");
    while True:
        while prev_score == ssorted[ii]["score"]:
            if ssorted[ii]["label"] == 1:
                tpcount += 1;
            else:
                fpcount += 1;
            ii+=1;
            if samplenum == ii:
                break;

        tpr = tpcount/float(positives);
        fpr = fpcount/float(negatives);
        points.append((tpr,fpr));
        if samplenum == ii:
            break;
        prev_score = ssorted[ii]["score"];
    points.append((1.0,1.0));
    ret = 0.0;
    for ii in range(len(points)-1):
        current_tpr = points[ii][0];
        current_fpr = points[ii][1];
        next_tpr = points[ii+1][0];
        next_fpr = points[ii+1][1];
        ret += (current_tpr+next_tpr)*(next_fpr-current_fpr)/2.0;
    return ret;

# {score:score,label:label}
def calc_pr_roc_auc(targets,reverse):
    ssorted = list(sorted(targets,key=lambda x:x["score"],reverse=reverse));
    prev_score = ssorted[0]["score"];
    samplenum = len(ssorted);
    positives = 0;
    for ss in list(ssorted):
        if ss["label"] == 1:
            positives += 1;
    negatives = samplenum - positives;
    if positives == 0 or negatives == 0:
        return float("nan");
    tpcount = 0;
    fpcount = 0;
    ii = 0;
    points = [];
    while True:
        while prev_score == ssorted[ii]["score"]:
            if ssorted[ii]["label"] == 1:
                tpcount += 1;
            else:
                fpcount += 1;
            ii+=1;
            if samplenum == ii:
                break;

        precision = tpcount/float(tpcount+fpcount);
        recall = tpcount/float(positives);
        if len(points) == 0:
            points.append((precision,0.0))
        points.append((precision,recall));
        if samplenum == ii:
            break;
        prev_score = ssorted[ii]["score"];
    ret = 0.0;
    for ii in range(len(points)-1):
        current_p = points[ii][0];
        current_r = points[ii][1];
        next_p = points[ii+1][0];
        next_r = points[ii+1][1];
        ret += (current_p+next_p)*(next_r-current_r)/2.0;
        
        
    return ret;


# {score:score,label:label}
def calc_average_precision_score(targets,reverse):
    ssorted = list(sorted(targets,key=lambda x:x["score"],reverse=reverse));
    prev_score = ssorted[0]["score"];
    samplenum = len(ssorted);
    positives = 0;
    for ss in list(ssorted):
        if ss["label"] == 1:
            positives += 1;
    negatives = samplenum - positives;
    if positives == 0 or negatives == 0:
        return float("nan");
    tpcount = 0;
    fpcount = 0;
    ii = 0;
    points = [];
    while True:
        while prev_score == ssorted[ii]["score"]:
            if ssorted[ii]["label"] == 1:
                tpcount += 1;
            else:
                fpcount += 1;
            ii+=1;
            if samplenum == ii:
                break;
        precision = tpcount/float(tpcount+fpcount);
        recall = tpcount/float(positives);
        points.append((precision,recall));
        if samplenum == ii:
            break;
        prev_score = ssorted[ii]["score"];
    ret = points[0][0]*points[0][1];
    for ii in range(len(points)-1):
        current_p = points[ii][0];
        current_r = points[ii][1];
        next_p = points[ii+1][0];
        next_r = points[ii+1][1];
        ret += next_p*(next_r-current_r);
        
    return ret;


def tester():
    from sklearn.metrics import roc_auc_score
    from sklearn.metrics import average_precision_score
    import random;
    for _ in range(100):
        chatgpt_sc = [];
        chatgpt_lab = [];
        mine = [];
        samevalue = {};
        for ii in range(1000):
            scc = round(random.random(),2);
            lab = 1 if random.random() > 0.5 else 0;
            mine.append({
                "score":scc,"label":lab
            });
            chatgpt_lab.append(lab);
            chatgpt_sc.append(scc);
        
        assert abs(calc_average_precision_score(mine,reverse=True)
         - average_precision_score(chatgpt_lab, chatgpt_sc)) < 0.0001;
        assert abs(calc_roc_auc(mine,reverse=True)
         - roc_auc_score(chatgpt_lab, chatgpt_sc)) < 0.0001;
    print("OK");

if sys.argv[1] == "test":
    tester();
    exit(0);



def get_groupids(gcode):
    groups = re.split(r"\.",gcode);
    assert len(groups) > 3,"??? invalid format. "+gcode+"\n";
    g1 = ".".join(groups);
    g2 = ".".join(groups[:-1]);
    g3 = ".".join(groups[:-2]);
    return (g1,g2,g3);

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

    g1_score_ave_prec = copy.deepcopy(group_1_orig);
    g2_score_ave_prec = copy.deepcopy(group_2_orig);
    g3_score_ave_prec = copy.deepcopy(group_3_orig);


    higher_is_better = False if "euc" in score_type else True;
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

            targetscores = list(sorted(targetscores_,key=lambda x: x[-1],reverse = higher_is_better));
            del targetscores_;

            allscores_g1 = [];
            allscores_g2 = [];
            allscores_g3 = [];
            for ptt in list(targetscores):
                if ptt[0] == name:
                    continue;
                ############# 上位グループが同じでも FP とみなす。
                bg1,bg2,bg3 = get_groupids(ptt[1]);
                
                if bg1 == g1:
                    tpcount_g1 += 1;
                    allscores_g1.append({"score":ptt[-1],"label":1});
                    continue;
                else:
                    fpcount_g1 += 1;
                    allscores_g1.append({"score":ptt[-1],"label":0});
                    if fpcount_g1 == 1:
                        g1_fp1 = tpcount_g1;
                    if fpcount_g1 == 10:
                        g1_fp10 = tpcount_g1;

                if bg2 == g2:
                    tpcount_g2 += 1;
                    allscores_g2.append({"score":ptt[-1],"label":1});
                    continue;
                else:
                    fpcount_g2 += 1;
                    allscores_g2.append({"score":ptt[-1],"label":0});
                    if fpcount_g2 == 1:
                        g2_fp1 = tpcount_g2;
                    if fpcount_g2 == 10:
                        g2_fp10 = tpcount_g2;
        
                if bg3 == g3:
                    tpcount_g3 += 1;
                    allscores_g3.append({"score":ptt[-1],"label":1});
                else:
                    fpcount_g3 += 1;
                    allscores_g3.append({"score":ptt[-1],"label":0});
                    if fpcount_g3 == 1:
                        g3_fp1 = tpcount_g3;
                    if fpcount_g3 == 10:
                        g3_fp10 = tpcount_g3;


            # 全 TP, FP を数えたので ROC-AUC を計算する
            roc_g1 = calc_roc_auc(allscores_g1,reverse=higher_is_better);
            roc_g2 = calc_roc_auc(allscores_g2,reverse=higher_is_better);
            roc_g3 = calc_roc_auc(allscores_g3,reverse=higher_is_better);

            # ave_prec を計算する
            ave_prec_g1 = calc_average_precision_score(allscores_g1,reverse=higher_is_better);
            ave_prec_g2 = calc_average_precision_score(allscores_g2,reverse=higher_is_better);
            ave_prec_g3 = calc_average_precision_score(allscores_g3,reverse=higher_is_better);

            if tpcount_g1 > 0:
                g1_count[g1] += 1;
                g1_score_fp1[g1] += g1_fp1/float(tpcount_g1);
                g1_score_fp10[g1] += g1_fp10/float(tpcount_g1);
                
                assert not math.isnan(roc_g1);
                assert not math.isnan(ave_prec_g1);

                g1_score_roc[g1] += roc_g1;
                g1_score_ave_prec[g1] += ave_prec_g1;

                if g1_fp1 > 0:
                    g1_score_hit1[g1] += 1.0;
                if g1_fp10 > 0:
                    g1_score_hit10[g1] += 1.0;
                    
            if tpcount_g2 > 0:
                g2_count[g2] += 1;
                g2_score_fp1[g2] += g2_fp1/float(tpcount_g2);
                g2_score_fp10[g2] += g2_fp10/float(tpcount_g2);

                assert not math.isnan(roc_g2);
                assert not math.isnan(ave_prec_g2);
                
                g2_score_roc[g2] += roc_g2;
                g2_score_ave_prec[g2] += ave_prec_g2;
                
                if g2_fp1 > 0:
                    g2_score_hit1[g2] += 1.0;
                if g2_fp10 > 0:
                    g2_score_hit10[g2] += 1.0;
                    

            if tpcount_g3 > 0:
                g3_count[g3] += 1;
                g3_score_fp1[g3] += g3_fp1/float(tpcount_g3);
                g3_score_fp10[g3] += g3_fp10/float(tpcount_g3);
                
                assert not math.isnan(roc_g3);
                assert not math.isnan(ave_prec_g3);

                g3_score_roc[g3] += roc_g3;
                g3_score_ave_prec[g3] += ave_prec_g3;
                
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
    score_sum_ave_prec = {}; 

    for (tag,count,score_fp1,score_fp10,score_hit1,score_hit10,score_roc,score_ave_prec) in [
        ("g1",g1_count,g1_score_fp1,g1_score_fp10,g1_score_hit1,g1_score_hit10,g1_score_roc,g1_score_ave_prec),
        ("g2",g2_count,g2_score_fp1,g2_score_fp10,g2_score_hit1,g2_score_hit10,g2_score_roc,g2_score_ave_prec),
        ("g3",g3_count,g3_score_fp1,g3_score_fp10,g3_score_hit1,g3_score_hit10,g3_score_roc,g3_score_ave_prec),
    ]:
        counted_class = 0;
        num_queries[tag]  = 0;
        
        score_sum_fp1[tag] = 0;
        score_sum_fp10[tag] = 0;
        score_sum_hit1[tag] = 0;
        score_sum_hit10[tag] = 0;
        score_sum_roc[tag] = 0;
        score_sum_ave_prec[tag] = 0;

        for gg in list(count.keys()):
            if count[gg] > 0:# 同一クラスに存在するメンバ数で割って平均化する
                counted_class += 1; 
                num_queries[tag] += count[gg];
                score_fp1[gg] /= count[gg];
                score_fp10[gg] /= count[gg];
                score_hit1[gg] /= count[gg];
                score_hit10[gg] /= count[gg];
                score_roc[gg] /= count[gg];
                score_ave_prec[gg] /= count[gg];

                score_sum_fp1[tag] += score_fp1[gg];
                score_sum_fp10[tag] += score_fp10[gg];
                score_sum_hit1[tag] += score_hit1[gg];
                score_sum_hit10[tag] += score_hit10[gg];
                score_sum_roc[tag] += score_roc[gg];
                score_sum_ave_prec[tag] += score_ave_prec[gg];
                
        # 全クラス数で割って更に平均化する
        score_sum_fp1[tag] /= counted_class;
        score_sum_fp10[tag] /= counted_class;
        score_sum_hit1[tag] /= counted_class;
        score_sum_hit10[tag] /= counted_class;
        score_sum_roc[tag] /= counted_class;
        score_sum_ave_prec[tag] /= counted_class;

    def strline(*argg):
        return "\t".join([str(x) for x in argg]);


    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g1","untilfp1","count:",num_queries["g1"],"score:",score_sum_fp1["g1"]));
    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g2","untilfp1","count:",num_queries["g2"],"score:",score_sum_fp1["g2"]));
    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g3","untilfp1","count:",num_queries["g3"],"score:",score_sum_fp1["g3"]));
    
    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g1","hit_at_1","count:",num_queries["g1"],"score:",score_sum_hit1["g1"]));
    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g1","hit_at_10","count:",num_queries["g1"],"score:",score_sum_hit10["g1"]));
    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g2","hit_at_1","count:",num_queries["g2"],"score:",score_sum_hit1["g2"]));
    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g2","hit_at_10","count:",num_queries["g2"],"score:",score_sum_hit10["g2"]));
    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g3","hit_at_1","count:",num_queries["g3"],"score:",score_sum_hit1["g3"]));
    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g3","hit_at_10","count:",num_queries["g3"],"score:",score_sum_hit10["g3"]));

    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g1","roc","count:",num_queries["g1"],"score:",score_sum_roc["g1"]));
    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g2","roc","count:",num_queries["g2"],"score:",score_sum_roc["g2"]));
    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g3","roc","count:",num_queries["g3"],"score:",score_sum_roc["g3"]));
    
    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g1","ave_prec","count:",num_queries["g1"],"score:",score_sum_ave_prec["g1"]));
    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g2","ave_prec","count:",num_queries["g2"],"score:",score_sum_ave_prec["g2"]));
    print(strline(targetdir,score_type,"higher_is_better:"+str(higher_is_better),"g3","ave_prec","count:",num_queries["g3"],"score:",score_sum_ave_prec["g3"]));
    
    
