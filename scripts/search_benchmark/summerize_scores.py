import re,os,sys;
import math;

# calc_scores.py の出力をパースしてどのファクターが影響しているか調べる
# 出力の形式が変更されるとこちらも変更する必要がある
# 次何か変更することがあったら汎用的にする


allfiles = re.split(r",",sys.argv[1]);

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

for catt in list(sorted(category_to_lines.keys())):
    allscores = [];
    for ll in list(category_to_lines[catt]):
        hs = line_to_hash(ll);
        ptt = re.split(r"[\s]+",ll);

        dpath = ptt[col_path];
        dpath = re.sub(r".*/","",dpath);

        xptt = re.split(r"_",dpath);
        plm = xptt[1];
        gnorm = xptt[2];
        
        xptt = re.split(r"\.",ptt[col_func]);
        pooling = xptt[0];
        ffunc = xptt[1];

        plms.add(plm);
        gnorms.add(gnorm);
        poolings.add(pooling);
        ffuncs.add(ffunc);

        tagmerged =  "all#"+plm+"#"+gnorm+"#"+pooling+"#"+ffunc;

        for pp in [plm,gnorm,pooling,ffunc,tagmerged]:
            if pp not in ranksum:
                ranksum[pp] = 0;
                deltaratiosum[pp] = 0;
                positive_zscoresum[pp] = 0;

        allscores.append(
            {
                "plm":plm,
                "gnorm":gnorm,
                "pooling":pooling,
                "ffunc":ffunc,
                "all":tagmerged,
                "score":float(hs["score"])
            }
        );

    allscores = list(sorted(allscores,key=lambda x:x["score"],reverse=True));
    
    ssum = 0;
    for ii in range(len(allscores)):
        ssum += allscores[ii]["score"];
    aave = ssum/float(len(allscores));

    vsum = 0;
    for ii in range(len(allscores)):
        vsum += (allscores[ii]["score"]-aave)*(allscores[ii]["score"]-aave);
    vvar = vsum/float(len(allscores));
    sstd = math.sqrt(vvar);

    maxscore = allscores[0]["score"];
    
    assert sstd != 0.0;

    for ii in range(len(allscores)):
        t = allscores[ii];
        deltaratio = (maxscore-t["score"])/maxscore;
        zsc = (t["score"]-aave)/sstd;

        for ttag in list(sep_tags):
            ranksum[t[ttag]] += ii+1;
            deltaratiosum[t[ttag]] += deltaratio; 
            if zsc > 0.0:
                positive_zscoresum[t[ttag]] += zsc;
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

    print("=======");
    for rr in r:
        print(cat,rr["label"],rr["score"]);
    print("=======");
    for rr in d:
        print(cat,rr["label"],rr["score"]);
    print("=======");
    for rr in z:
        print(cat,rr["label"],rr["score"]);
    print("+++++++++++++");
    print("+++++++++++++");

