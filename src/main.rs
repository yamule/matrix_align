use std::collections::*;
use matrix_align::{matrix_process, simple_argparse};
use rayon;
use rayon::prelude::*;
#[allow(unused_imports)]
use matrix_align::gmat::{self, calc_vec_stats, calc_vec_stats_legacy, GMatStatistics};
use matrix_align::aligner::{AlignmentType, GapPenaltyAutoAdjustParam, ProfileAligner, ScoreType, SequenceProfile};
use matrix_align::ioutil::{self, load_multi_gmat, save_lines};
use matrix_align::guide_tree_based_alignment::{self, DistanceBase};
use matrix_align::a3m_pairwise_alignment;
use rand::SeedableRng;
use rand::rngs::StdRng;
use regex::Regex;
use matrix_align::misc::*;

fn main(){
    main_(std::env::args().collect::<Vec<String>>());
}

fn main_(mut args:Vec<String>){

    assert!(args.len() > 0);
    if args[0].contains("matrix_align"){
        let _fst = args.remove(0);
    }
    
    let allowed_arg_:Vec<(&str,Option<&str>,&str,Option<&str>,Vec<&str>,bool)> = vec![
        ("--in",None
        ,"<input file path> : Text based file which contains general profile matrices for multiple sequences. Something \
        like as follows. \n>seq1\n[amino acid letter at position 1]\t[value1]\t[value2]\t[value3]...\n[amino acid letter \
        at position 2]\t[value1]\t[value2]\t[value3]...\n[amino acid letter at position 3]\t[value1]\t[value2]\t[value3]\
        ...\n...\n>seq2\n[amino acid letter at position 1]\t[value1]\t[value2]\t[value3]...\n[amino acid letter at posit\
        ion 2]\t[value1]\t[value2]\t[value3]...\n[amino acid letter at position 3]\t[value1]\t[value2]\t[value3]...\n...\n\
        ...\n"
        ,None,vec![], false),
        
        ("--in_list",None,
        "<list file path> : Ascii text file which contains multiple input files's path. If --in was specified, --in files will be loaded \
        at first and files in --in_list will be loaded the next."
        ,None,vec![],false
        ),

        ("--in_stats",None,
        "<stats file path> : Import statistics from this file instead of calculating them on the fly."
        ,None,vec![],false
        ),

        ("--out",None
        ,"<output file path> : Output file for alignemt result in multi-fasta format."
        ,None,vec![], false),

        ("--out_stats",None
        ,"<output file path> : Output file for statistics which will be used at normalization process. Program stops after calculating statistical values."
        ,None,vec![], false),
        
        ("--num_iter",None
        ,"<int> : Number of alignment iterations. Construct global profile with -1 with this number. Can not be used with tree_guided alignment."
        ,Some("1"),vec![],false),
        
        ("--gap_open_penalty",None
        ,"<float> : Gap open penalty for DP. Must be negative."
        ,Some("-10.0"),vec![],false),
        
        ("--gap_penalty_auto_adjust",None
        ,"<bool or novalue=true> : Adjust gap penalty automaticall. <gap open penalty> = <maximum matching score>*a1*-1.0 + \
        <minimum matching score>*a2; <gap extension penalty> = <gap open penalty>*0.05;"
        ,Some("true"),vec![],false),
        
        ("--gap_penalty_a1_a2",None
        ,"<float>,<float> : Parameters for gap penalty auto adjust."
        ,Some("0.5,0.5"),vec![],false),

        ("--normalize",None
        ,"<bool or novalue=true> : Normalize profile values before alignment."
        ,Some("true"),vec![],false),
        
        ("--alignment_type",None
        ,"<string> : Alignment type."
        ,Some("global"),vec!["global","local"],false),
        
        ("--num_threads",None
        ,"<int> : Number of threads."
        ,Some("4"),vec![],false),
        
        ("--tree_guided",None
        ,"<bool or novalue=true> : Use guide-tree based alignment."
        ,Some("true"),vec![],false),

        ("--a3m_pairwise",None
        ,"<bool or novalue=true> : Perform pairwise alignment with top sequense and other sequences and output alignments as A3M format."
        ,Some("false"),vec![],false),
        
        ("--distance_base",None
        ,"<string> : Estimation method for distance between samples."
        ,Some("averaged_value"),vec!["averaged_value","score_with_seed"],false),

        ("--max_cluster_size",None
        ,"<int> : Limit all-vs-all comparison & tree building with this many number of profiles and align hierarchically. Must be > 10 or -1 (no limit)."
        ,Some("-1"),vec![],false),
        
        ("--random_seed",None
        ,"<int> : Seed for random number generator. Must be positive."
        ,None,vec![],false),
        
        ("--tree_type",None
        ,"<string> : Type of guide tree."
        ,Some("NJ"),vec!["NJ","UPGMA"],false),
        
        ("--score_type",None
        ,"<string> : Type of scoring function. \"distance_zscore\" or \"dot_product\"."
        ,Some("dot_product"),vec![],false),
        
        ("--out_matrix",None
        ,"<string> : Save matrix file of the resulting alignment."
        ,None,vec![],false),

        ("--parallel_load",None
        ,"<bool or novalue=true> : Load multiple files parallely. Set False if each file is too large."
        ,Some("true"),vec![],false),
    ];

    for aa in args.iter(){
        if aa == "--help" || aa == "-h"{
            matrix_process::check_simd();
        }
    }

    let mut argparser = simple_argparse::SimpleArgParse::new(allowed_arg_);
    argparser.parse(args);

    let tree_type = if argparser.get_string("--tree_type").unwrap() == "NJ"{
        guide_tree_based_alignment::TreeType::TreeNj
    }else{
        assert!(argparser.get_string("--tree_type").unwrap() == "UPGMA");
        guide_tree_based_alignment::TreeType::TreeUPGMA
    };

    let mut rngg:StdRng = if let Some(x) = argparser.get_int("--random_seed"){
        assert!(x >  0,"--random_seed must be a positive integer. {}",x);
        StdRng::seed_from_u64(x as u64)
    }else{
        StdRng::from_entropy()
    };

    let max_cluster_size = argparser.get_int("--max_cluster_size").unwrap() as i64;
    if max_cluster_size > -1 {
        assert!(max_cluster_size >  10,"--max_cluster_size must be > 10 or -1. {}",max_cluster_size);
    }// -1 の場合 hierarcical align を行わない

    
    let voidpair:Vec<(&str,&str)> = vec![
        ("--a3m_pairwise","--tree_guided"),
        ("--in_stats","--out_stats"),
        ("--a3m_pairwise","--distance_base"),
        ("--num_iter","--tree_guided"),
        ("--gap_penalty_auto_adjust","--gap_open_penalty"),
        ("--gap_penalty_a1_a2","--gap_open_penalty"),
    ];
    for (a,b) in voidpair{
        if argparser.is_generous_false(a){
            continue;
        }
        if argparser.is_generous_false(b){
            continue;
        }
        if argparser.user_defined(a) && argparser.user_defined(b){
            panic!("{} and {} can not be used at the same time.",a,b);
        }    
    }

    let reqpair:Vec<(&str,&str)> = vec![
        ("--in","--in_list"),
        ("--out","--out_stats"),
    ];
    for (a,b) in reqpair{
        if argparser.is_generous_false(a) && argparser.is_generous_false(b) {
            panic!("One of {} or {} should be set.",a,b);
        }
    }

    argparser.print_items();
    
    let alignment_type;
    let typ = argparser.get_string("--alignment_type").unwrap().to_lowercase();
    if typ.as_str() == "global"{
        alignment_type = AlignmentType::Global;
    }else if typ.as_str() == "local"{
        alignment_type = AlignmentType::Local;
    }else{
        panic!("???");
    }
    
    let distance_base;
    let typ = argparser.get_string("--distance_base").unwrap().to_lowercase();
    if typ.as_str() == "score_with_seed"{
        distance_base = DistanceBase::ScoreWithSeed;
    }else if typ.as_str() == "averaged_value"{
        distance_base = DistanceBase::AveragedValue;
    }else{
        panic!("???");
    }
    
    let score_type;
    let typ = argparser.get_string("--score_type").unwrap().to_lowercase();
    if typ.as_str() == "distance_zscore"{
        score_type = ScoreType::DistanceZscore;
    }else if typ.as_str() == "dot_product"{
        score_type = ScoreType::DotProduct;
    }else{
        panic!("???");
    }

    let mut name_to_res:HashMap<String,String> = HashMap::new();
    
    let mut infiles:Vec<String> =  vec![];
    if let Some(x) = argparser.get_string("--in"){
        for xx in x.split(",").into_iter().map(|m|m.to_string()).into_iter(){
            infiles.push(xx);
        }
    }

    if let Some(x) = argparser.get_string("--in_list"){
        let re = Regex::new(r"[\s]+$").unwrap();
        for xx in x.split(",").into_iter().map(|m|m.to_string()).into_iter(){
            let q:Vec<String> = ioutil::load_lines(&xx,xx.ends_with(".gz"));
            for qq_ in q.into_iter(){
                let qq = re.replace(&qq_,"").to_string();
                if qq.len() == 0{
                    continue;
                }
                if qq.starts_with("#"){
                    continue;
                }
                infiles.push(qq);
            }
        }
    }


    let mut gmatstats:Vec<GMatStatistics>;
    let mut profile_seq:Option<SequenceProfile> = None;
    if let Some(x) = argparser.get_string("--in_stats"){
        let lines = ioutil::load_lines(&x,x.ends_with(".gz"));
        gmatstats = vec![];
        for ll in lines.iter(){
            if ll.len() < 2{
                continue;
            }
            if ll.starts_with("#"){
                continue;
            }
            gmatstats.push(
                GMatStatistics::load_string(ll)
            );
        }
    }else{
        unsafe{
            gmatstats = calc_vec_stats(&infiles);
        }
        if let Some(x) = argparser.get_string("--out_stats"){
            let mut lines:Vec<String> = vec![];
            for gg in gmatstats.iter(){
                lines.push(gg.get_string());
            }
            ioutil::save_lines(&x, lines, x.ends_with(".gz"));
            if !argparser.user_defined("--out"){
                std::process::exit(0);
            }
        }
    }

    let mut name_ordered:Vec<String> = vec![];
    let mut gmat1_:Vec<(String,Vec<char>,Vec<Vec<f32>>,Option<Vec<(f32,f32,f32,f32)>>)> = vec![];
    if argparser.is_generous_false("--a3m_pairwise"){
        for ii in infiles.iter(){
            let mut z = load_multi_gmat(ii,ii.ends_with(".gz"));
            gmat1_.append(&mut z);
        }
    }else{
        // a3m pairwise の場合はまずファイル一つだけのロードでよい
        let fname = infiles.remove(0);
        let mut z = load_multi_gmat(&fname,fname.ends_with(".gz"));
        gmat1_.append(&mut z);
    }

    let veclen = gmat1_[0].2[0].len();
    let mut saligner:ProfileAligner = if argparser.get_bool("--gap_penalty_auto_adjust").unwrap(){
        let gap_penalty_a1_a2:Vec<String> = argparser.get_string("--gap_penalty_a1_a2").unwrap().to_string().split(',').map(|m|m.to_owned()).collect();
        ProfileAligner::new(veclen,300,None
        ,alignment_type,score_type,Some(GapPenaltyAutoAdjustParam{
            a1:gap_penalty_a1_a2[0].parse::<f32>().unwrap_or_else(|e| panic!("{:?} {:?}",gap_penalty_a1_a2,e)),
            a2:gap_penalty_a1_a2[1].parse::<f32>().unwrap_or_else(|e| panic!("{:?} {:?}",gap_penalty_a1_a2,e))
        }))
    }else{
        ProfileAligner::new(veclen,300, Some(argparser.get_float("--gap_open_penalty").unwrap() as f32)
        ,alignment_type,score_type,None)
    };


    let num_threads:usize = argparser.get_int("--num_threads").unwrap() as usize;
    assert!(num_threads > 0);
    match rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global(){
        Ok(_)=>{

        },
        Err(e)=>{
            eprintln!("=====If you are in test process, this can be ignored.=====");
            eprintln!("{:?}",e);
            eprintln!("==========================================================");
        }
    }


    let seqprepare = |mut gg:(String,Vec<char>,Vec<Vec<f32>>,Option<Vec<(f32,f32,f32,f32)>>)
        ,name_to_res:&mut HashMap<String,String>,name_ordered: &mut Vec<String>|
         -> SequenceProfile{
        let n = gg.0.clone();
        if name_to_res.contains_key(&n){
            panic!("Duplicated name {}.",n);
        }
        name_to_res.insert(n.clone(),"".to_owned());
        name_ordered.push(n);
        if argparser.get_bool("--normalize").unwrap(){
            gmat::normalize_seqmatrix(&mut (gg.2), &gmatstats);
        }
        let alen = gg.1.len();
        let seq = SequenceProfile::new(
            vec![(gg.0,gg.1)],alen,gg.2[0].len(),None,Some(gg.2),gg.3
        );
        return seq;
    };

    if argparser.get_bool("--a3m_pairwise").unwrap(){
        
        let mut firstseq_ = gmat1_.remove(0);

        let mut aligners:Vec<ProfileAligner> = vec![];
        for _ in 0..num_threads{
            aligners.push(saligner.clone());
        }

        if argparser.get_bool("--normalize").unwrap(){
            gmat::normalize_seqmatrix(&mut (firstseq_.2), &gmatstats);
        }
        let alen = firstseq_.1.len();
        let firstseq = SequenceProfile::new(
            vec![(firstseq_.0,firstseq_.1)],alen,firstseq_.2[0].len(),None,Some(firstseq_.2),firstseq_.3
        );

        let mut lines:Vec<String> = vec![];
        lines.push(
            ">".to_owned()+firstseq.headers[0].as_str()
        );
        lines.push(
            firstseq.member_sequences[0].iter().filter(|c| **c != '-').map(|m| m.to_string()).collect::<Vec<String>>().join("")
        );

        let mut name_length_order:Vec<(String,usize)> = vec![];
        
        let mut allseqs_:VecDeque<SequenceProfile> = VecDeque::new();
        for gg in gmat1_.into_iter(){
            allseqs_.push_back(seqprepare(gg,&mut name_to_res,&mut name_ordered));
        }

        infiles.reverse();//pop するので逆順にする
        while infiles.len() > 0 || allseqs_.len() > 0{
            while infiles.len() > 0{
                let mut z = vec![];
                if argparser.get_bool("--parallel_load").unwrap(){
                    let mut loader_minibatch:Vec<(usize,String)> = vec![];
                    while infiles.len() >0{ 
                        let ii = infiles.pop().unwrap();
                        loader_minibatch.push((loader_minibatch.len(),ii));
                        if loader_minibatch.len() >= num_threads{
                            break;
                        }
                    }
                    let mut res:Vec<(usize,Vec<(String,Vec<char>,Vec<Vec<f32>>,Option<Vec<(f32,f32,f32,f32)>>)>)> = loader_minibatch.into_par_iter().map(|v|{
                        let fileindex = v.0;
                        let filename = v.1;
                        let pres = load_multi_gmat(&filename,filename.ends_with(".gz"));
                        return (fileindex,pres);
                    }).collect();
                    res.sort_by(|a,b| (a.0.cmp(&b.0)));
                    for mut rr in res.into_iter(){
                        z.append(&mut rr.1);
                    }
                }else{
                    let ii = infiles.pop().unwrap();
                    z = load_multi_gmat(&ii,ii.ends_with(".gz"));
                }
                for seq in z.into_iter(){
                    allseqs_.push_back(
                        seqprepare(seq,&mut name_to_res,&mut name_ordered)
                    );
                }
                if allseqs_.len() >= num_threads{ //用意された Threads 分はメモリに乗る想定
                    break;
                }
            }
            
            let mut allseqs:Vec<SequenceProfile> = vec![];
            
            while allseqs_.len() > 0{
                let ss = allseqs_.pop_front().unwrap();

                assert!(ss.headers.len() == 1);
                assert!(ss.member_sequences.len() == 1);// 初期値ギャップも match_scores として計算される想定
                name_length_order.push(
                    (ss.headers[0].clone(),ss.member_sequences[0].len())
                );
                allseqs.push(ss);
                if allseqs.len() >= num_threads{
                    break;
                }
            }

            let mut res = a3m_pairwise_alignment::create_a3m_pairwise_alignment(&mut aligners,firstseq.clone(),allseqs);

            for nn in name_length_order.iter(){
                if let Some(p) = res.remove(&nn.0){
                    lines.push(">".to_owned()+&nn.0);
                    lines.push(p.0);
                    println!(">{}",nn.0);
                    println!("score:{}",p.1.score);
                    let mut posicount = 0 as usize;
                    for dd in p.1.match_scores.iter(){
                        if *dd > 0.0{
                            posicount += 1;
                        }
                    }
                    println!("positive_count:{}",posicount);
                    println!("profile_length:{}",nn.1);
                }
            }
            assert!(res.len() == 0);
        }
        //結果は全部メモリに乗る想定
        let outfile = argparser.get_string("--out").unwrap();
        save_lines(&outfile, lines,outfile.ends_with(".gz"));
        std::process::exit(0);   
    }
    
    let mut allseqs_:Vec<SequenceProfile> = vec![];
    for mut gg in gmat1_.into_iter(){
        let n = gg.0.clone();
        if name_to_res.contains_key(&n){
            panic!("Duplicated name {}.",n);
        }
        name_to_res.insert(n.clone(),"".to_owned());
        name_ordered.push(n);
        if argparser.get_bool("--normalize").unwrap(){
            gmat::normalize_seqmatrix(&mut (gg.2), &gmatstats);
        }
        let alen = gg.1.len();
        let seq = SequenceProfile::new(
            vec![(gg.0,gg.1)],alen,gg.2[0].len(),None,Some(gg.2),gg.3
        );
        allseqs_.push(seq);
    }

    

    let num_iter:usize = argparser.get_int("--num_iter").unwrap() as usize;
    for ii in 0..num_iter{
        eprintln!("iter: {}",ii);

        let mut seqvec:Vec<SequenceProfile> = vec![];
        
        if let Some(p) = profile_seq{
            seqvec.push(p.create_merged_dummy());
        }

        for ss in allseqs_.iter(){
            seqvec.push(ss.clone());
        }
        
        let mut ans = if argparser.get_bool("--tree_guided").unwrap(){
            if max_cluster_size == -1{
                guide_tree_based_alignment::tree_guided_alignment(seqvec, &distance_base,&mut saligner,false,num_threads,tree_type,&mut rngg)
            }else{
                guide_tree_based_alignment::hierarchical_alignment(seqvec, &distance_base,&mut saligner
                    ,max_cluster_size, &mut rngg,num_threads,tree_type)
            }
        }else{
            saligner.make_msa(seqvec,false)
        };

        assert!(ans.len() == 1);
        println!("score:{}",ans[0].1);
        
        let (alires,_alisc) = ans.pop().unwrap();
        
        //let pstr = alires.gmat_str();
        //for pp in pstr{ //profile 表示
        //    println!("{}",pp);
        //}

        if num_iter == ii+1{
            insert_alinged_string(&alires, &mut name_to_res,true);
            if let Some(x) = argparser.get_string("--out_matrix"){
                ioutil::save_gmat(&x,&vec![(0,&alires)], x.ends_with(".gz"));
            }
            break;
        }
        profile_seq = Some(alires);
    }
    let mut results:Vec<String> = vec![];
    for ii in name_ordered.iter(){
        //println!(">{}",primary_id_to_name.get(ii).unwrap());
        //println!("{}",primary_id_to_res.get(ii).unwrap());
        results.push(
            format!(">{}\n{}",ii,name_to_res.get(ii).unwrap())
        );
    }
    let outfile = argparser.get_string("--out").unwrap();
    save_lines(&outfile, results,outfile.ends_with(".gz"));
    //stattest();
}


#[test]
fn maintest(){
    main_(
        (vec![
            "--in","example_files/test1.gmat,example_files/test2.gmat,example_files/test3.gmat,example_files/test4.gmat","--out","nogit/normalex.dat"
        ]).into_iter().map(|m|m.to_owned()).collect()
    );

    main_(
        (vec![
            "--in","example_files/test1.gmat","--out","nogit/list234ex.dat.gz"
            ,"--in_list","example_files/list234.dat"
        ]).into_iter().map(|m|m.to_owned()).collect()
    );

    main_(
        (vec![
            "--out","nogit/list1234ex.dat.gz"
            ,"--in_list","example_files/list1234.dat"
        ]).into_iter().map(|m|m.to_owned()).collect()
    );
    let v1 = ioutil::load_lines("nogit/normalex.dat",false);
    let v2 = ioutil::load_lines("nogit/list234ex.dat.gz",true);
    let v3 = ioutil::load_lines("nogit/list1234ex.dat.gz",true);
    assert_eq!(v1,v2);
    assert_eq!(v1,v3);
}