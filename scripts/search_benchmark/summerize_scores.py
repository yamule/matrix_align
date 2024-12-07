import re,os,sys;
import math;
import copy;

# calc_scores.py の出力をパースしてどのファクターが影響しているか調べる
# 出力の形式が変更されるとこちらも変更する必要がある
# 次何か変更することがあったら汎用的にする

allfiles = re.split(r",",sys.argv[1]);

out_prefix = sys.argv[2];

# 後で使うために書いただけで現在のところあまり有用でない
def line_to_hash(l):
    l = re.sub(r"[\r\n]","",l);
    l = re.sub(r"[\s]*:[\s]*",":",l);
    ptt = re.split(r"[\s]+",l);
    ret = {};
    ret["notag"] = [];
    for pp in ptt:
        mat = re.search(r"^([^:]+):([^:]+)",pp);
        if mat:
            k = mat.group(1);
            v = mat.group(2);
            if k in ret:
                raise Exception("Duplicated tag:"+k);
            ret[k] = v;
        else:
            ret["notag"].append(pp);
    return ret;

col_path = 0;
col_func = 1;
col_group = 3;
col_metrics = 4;
col_score = -1;


category_to_lines = {}; # 一旦カテゴリ (グループ+metrics) だけでまとめる

for ff in list(allfiles):
    with open(ff) as fin:
        for ll in fin:
            if "score" not in ll:
                continue;
            ll = re.sub(r"[\r\n]","",ll);
            ptt = re.split(r"[\s]+",ll);
            catt = ptt[col_group]+"\t"+ptt[col_metrics];
            if catt not in category_to_lines:
                category_to_lines[catt] = [];
            category_to_lines[catt].append(ll);


ranksum = {"all":0};
# sum((maxscore-target_score)/maxscore)
deltaratiosum = {"all":0.0};

positive_zscoresum = {"all":0.0};

plms = set();
poolings = set();
ffuncs = set();
gnorms = set();

sep_tags = [
"plm","pooling","ffunc","gnorm","all"
];
top5s = [];
allscores = [];

with open(out_prefix+"all_scores_tab.dat","wt") as fout:
    fout.write("category\tmetrics\tplm\tper_channel_normalization\tpooling_function\tscoring_function\tscore"+"\n")
    for catt in list(sorted(category_to_lines.keys())):
        catscores = [];
        ptt = re.split(r"[\s]+",catt)
        category = ptt[0];
        metrics = ptt[1];
        for ll in list(category_to_lines[catt]):
            hs = line_to_hash(ll);
            ptt = re.split(r"[\s]+",ll);

            dpath = ptt[col_path];
            dpath = re.sub(r".*/","",dpath);

            xptt = re.split(r"_",dpath);
            plm = xptt[1];
            gnorm = xptt[2];
            
            if gnorm == "split100":
                gnorm = xptt[3];

            xptt = re.split(r"\.",ptt[col_func]);
            pooling = xptt[0];
            ffunc = xptt[1];

            plms.add(plm);
            gnorms.add(gnorm);
            poolings.add(pooling);
            ffuncs.add(ffunc);

            tagmerged =  "all#"+plm+"#"+gnorm+"#"+pooling+"#"+ffunc;
            fout.write("\t".join(
                [category,metrics,plm,gnorm,pooling,ffunc,str(hs["score"])]
            )+"\n");
            for pp in [plm,gnorm,pooling,ffunc,tagmerged]:
                if pp not in ranksum:
                    ranksum[pp] = 0;
                    deltaratiosum[pp] = 0;
                    positive_zscoresum[pp] = 0;

            catscores.append(
                {
                    "category":category,
                    "metrics":metrics,
                    "plm":plm,
                    "gnorm":gnorm,
                    "pooling":pooling,
                    "ffunc":ffunc,
                    "all":tagmerged,
                    "score":float(hs["score"])
                }
            );

        catscores = list(sorted(catscores,key=lambda x:x["score"],reverse=True));
        allscores.extend(catscores);
        ssum = 0;
        for ii in range(len(catscores)):
            ssum += catscores[ii]["score"];
        aave = ssum/float(len(catscores));

        vsum = 0;
        for ii in range(len(catscores)):
            vsum += (catscores[ii]["score"]-aave)*(catscores[ii]["score"]-aave);
        vvar = vsum/float(len(catscores));
        sstd = math.sqrt(vvar);

        maxscore = catscores[0]["score"];
        
        assert sstd != 0.0;

        for ii in range(len(catscores)):
            t = catscores[ii];
            deltaratio = (maxscore-t["score"])/maxscore;
            zsc = (t["score"]-aave)/sstd;

            for ttag in list(sep_tags):
                ranksum[t[ttag]] += ii+1;
                deltaratiosum[t[ttag]] += deltaratio; 
                if zsc > 0.0:
                    positive_zscoresum[t[ttag]] += zsc;



allscores =  list(sorted(allscores,key=lambda x:x["score"],reverse=True));
top10s = {};
for aa in allscores:
    tag = aa["category"]+"#"+aa["metrics"];
    if tag not in top10s:
        top10s[tag] = [];
    if len(top10s[tag]) >= 10:
        continue;
    top10s[tag].append(aa);

top10tags = [
"g1#roc","g2#roc","g3#roc"
,"g1#hit_at_1","g2#hit_at_1","g3#hit_at_1"
,"g1#hit_at_10","g2#hit_at_10","g3#hit_at_10"
,"g1#untilfp1","g2#untilfp1","g3#untilfp1"
,"g1#ave_prec","g2#ave_prec","g3#ave_prec"
];
with open(out_prefix+"top10.dat","wt")as fout:
    fout.write("\t".join(top10tags)+"\n")
    for ii in range(10):
        for ss in list(top10tags):
            atag = [];
            for tt in ["category","metrics","plm","gnorm","pooling","ffunc"]:
                atag.append(top10s[ss][ii][tt]);
            fout.write("#".join(atag)+"({:.4f})".format(top10s[ss][ii]["score"])+"\t")
        fout.write(""+"\n")


basekeys = ["category","plm","gnorm","pooling","ffunc","metrics"];
for basekey in list(basekeys):
    scores_met = {};
    allmet_ = set();
    tag_score = {};

    sepkeys = copy.deepcopy(basekeys);
    sepkeys.remove(basekey);
    for aa in list(allscores):
        allmet_.add(aa[basekey]);
        ak = "#".join([aa[ss] for ss in list(sepkeys)]);
        if ak not in tag_score:
            tag_score[ak] = {};
        assert aa[basekey] not in tag_score[ak],basekey+"\t"+ak;
        tag_score[ak][aa[basekey]] = aa["score"];
        
    allmet = list(sorted(allmet_));
    # 全部表示
    with open(out_prefix+basekey+".table.dat","wt") as fout:
        fout.write(basekey+"\t"+"\t".join(allmet)+"\n");
        for tt in list(sorted(tag_score.keys())):
            for kk in list(allmet):
                if kk not in tag_score[tt]:
                    sys.stderr.write(kk+" was not found in "+tt+".\n"+str(tag_score[tt].keys())+"\n")
                    raise Exception();
            fout.write(tt+"\t"+"\t".join([str(tag_score[tt][kk]) for kk in list(allmet)])+"\n")

assert basekey == "metrics"; # 後の計算で使う
# gnorm あり無し ======================================
g1roc = 0;
g1rocv = 0;
g2roc = 0;
g2rocv = 0;
g3roc = 0;
g3rocv = 0;

g1ap = 0;
g1apv = 0;
g2ap = 0;
g2apv = 0;
g3ap = 0;
g3apv = 0;
with open(out_prefix+"global_per_channel_table.dat","wt") as fout:
    fout.write("\t".join(["tag","nognorm_roc","wgnorm_roc","nognorm_ap","wgnorm_ap"])+"\n")
    for tt in list(sorted(tag_score.keys())):
        if "wgnorm" in tt:
            continue;
        wgtt = re.sub(r"nognorm","wgnorm",tt);
        print(tt+"\t"+"\t".join([
            str(tag_score[tt]["roc"]),str(tag_score[wgtt]["roc"]),str(tag_score[tt]["ave_prec"]),str(tag_score[wgtt]["ave_prec"])
        ]))

        if "g1#" in tt:
            if tag_score[tt]["roc"] < tag_score[wgtt]["roc"]:
                g1roc += 1;
            elif tag_score[tt]["roc"] > tag_score[wgtt]["roc"]:
                g1rocv += 1;

        if "g2#" in tt:
            if tag_score[tt]["roc"] < tag_score[wgtt]["roc"]:
                g2roc += 1;
            elif tag_score[tt]["roc"] > tag_score[wgtt]["roc"]:
                g2rocv += 1;
        
        if "g3#" in tt:
            if tag_score[tt]["roc"] < tag_score[wgtt]["roc"]:
                g3roc += 1;
            elif tag_score[tt]["roc"] > tag_score[wgtt]["roc"]:
                g3rocv += 1;
        
        if "g1#" in tt:
            if tag_score[tt]["ave_prec"] < tag_score[wgtt]["ave_prec"]:
                g1ap += 1;
            elif tag_score[tt]["ave_prec"] > tag_score[wgtt]["ave_prec"]:
                g1apv += 1;
                
        if "g2#" in tt:
            if tag_score[tt]["ave_prec"] < tag_score[wgtt]["ave_prec"]:
                g2ap += 1;
            elif tag_score[tt]["ave_prec"] > tag_score[wgtt]["ave_prec"]:
                g2apv += 1;
                
        if "g3#" in tt:
            if tag_score[tt]["ave_prec"] < tag_score[wgtt]["ave_prec"]:
                g3ap += 1;
            elif tag_score[tt]["ave_prec"] > tag_score[wgtt]["ave_prec"]:
                g3apv += 1;
print("# gnorm is effective =====")
print("roc\t"+"\t".join(str(r) for r in[
    g1roc,g2roc,g3roc
]));
print("rocv\t"+"\t".join(str(r) for r in[
    g1rocv,g2rocv,g3rocv
]));
print("ap\t"+"\t".join(str(r) for r in[
    g1ap,g2ap,g3ap
]));
print("apv\t"+"\t".join(str(r) for r in[
    g1apv,g2apv,g3apv
]));

# ===================



with open(out_prefix+"ranksum.dat","wt") as fout:
    for cat,labelset in [
        ("plm",plms),
        ("pooling",poolings),
        ("global_normalization",gnorms),
        ("scoring_function",ffuncs),
        ("all",None),
    ]:
        r = [];
        d = [];
        z = [];
        for rr in list(ranksum.keys()):      
            if labelset is not None:
                if rr in labelset:
                    r.append({"label":rr,"score":ranksum[rr]});
                    d.append({"label":rr,"score":deltaratiosum[rr]});
                    z.append({"label":rr,"score":positive_zscoresum[rr]});
            else:
                if "all#" in rr:
                    r.append({"label":rr,"score":ranksum[rr]});
                    d.append({"label":rr,"score":deltaratiosum[rr]});
                    z.append({"label":rr,"score":positive_zscoresum[rr]});
                
            r = list(sorted(r,key=lambda x:x["score"],reverse=False));
            d = list(sorted(d,key=lambda x:x["score"],reverse=False));
            z = list(sorted(z,key=lambda x:x["score"],reverse=True));

        fout.write("======="+"\n");
        for rr in r:
            fout.write(cat+"\t"+rr["label"]+"\t"+str(rr["score"])+"\n");
        fout.write("======="+"\n");
        for rr in d:
            fout.write(cat+"\t"+rr["label"]+"\t"+str(rr["score"])+"\n");
        fout.write("======="+"\n");
        for rr in z:
            fout.write(cat+"\t"+rr["label"]+"\t"+str(rr["score"])+"\n");
        fout.write("+++++++++++++"+"\n");
        fout.write("+++++++++++++"+"\n");

