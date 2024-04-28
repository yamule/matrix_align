use std::collections::*;
use matrix_align::{matrix_process, simple_argparse};
use rayon;
use matrix_align::gmat::{self, calc_vec_stats, GMatStatistics};
use matrix_align::aligner::{AlignmentType, GapPenaltyAutoAdjustParam, ProfileAligner, ScoreType, SequenceProfile};
use matrix_align::ioutil::{self, load_multi_gmat, save_lines};
use matrix_align::guide_tree_based_alignment::{self, DistanceBase};
use rand::SeedableRng;
use rand::rngs::StdRng;


fn insert_alinged_string(alires:&SequenceProfile,name_to_res:&mut HashMap<String,String>){
    let mut maxpos:usize = 0;
    for mm in alires.alignment_mapping.iter(){
        for mmm in mm.iter(){
            if mmm.1 > -1{
                maxpos = maxpos.max(mmm.1 as usize);
            }
        }
    }
    maxpos += 1;

    for seqidx in 0..alires.member_sequences.len(){
        let mut aseq = alires.get_aligned_seq(seqidx);
        assert!(aseq.len() <= maxpos,"{} {} \n{}",aseq.len(),maxpos,aseq.iter().map(|m| m.to_string()).collect::<Vec<String>>().join(""));
        while aseq.len() < maxpos{
            aseq.push('-');
        }
        let hh = &alires.headers[seqidx];
        if name_to_res.contains_key(hh){
            assert!(name_to_res.get(hh).unwrap().len() == 0,"{}",name_to_res.get(hh).unwrap());
            name_to_res.insert(
                hh.clone(),aseq.into_iter().map(|m|m.to_string()).collect::<Vec<String>>().concat()
            );
        }
    }
    
}


fn main(){
    main_(std::env::args().collect::<Vec<String>>());
}
fn main_(args:Vec<String>){
    let allowed_arg_:Vec<(&str,Option<&str>,&str,Option<&str>,Vec<&str>,bool)> = vec![
        ("--in",None
        ,"<input file path> : Text based file which contains general profile matrices for multiple sequences. Something \
        like as follows. \n>seq1\n[amino acid letter at position 1]\t[value1]\t[value2]\t[value3]...\n[amino acid letter \
        at position 2]\t[value1]\t[value2]\t[value3]...\n[amino acid letter at position 3]\t[value1]\t[value2]\t[value3]\
        ...\n...\n>seq2\n[amino acid letter at position 1]\t[value1]\t[value2]\t[value3]...\n[amino acid letter at posit\
        ion 2]\t[value1]\t[value2]\t[value3]...\n[amino acid letter at position 3]\t[value1]\t[value2]\t[value3]...\n...\n\
        ...\n"
        ,None,vec![], true),
        
        ("--out",None
        ,"<output file path> : Output file for alignemt result in multi-fasta format."
        ,None,vec![], true),
        
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
        ,"<int> : Limit all-vs-all comparison with this many number of profiles and align hierarchically. Must be > 10."
        ,Some("100"),vec![],false),
        
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
        ,None,vec![],false)
    ];

    for aa in args.iter(){
        if aa == "--help" || aa == "-h"{
            unsafe{matrix_process::check_simd();}
        }
    }

    let mut argparser = simple_argparse::SimpleArgParse::new(allowed_arg_);
    argparser.parse(args);

    let tree_type = argparser.get_string("--tree_type").unwrap();
    
    let infiles_ = argparser.get_string("--in").unwrap().clone();
    let outfile = argparser.get_string("--out").unwrap().clone();
    let num_iter:usize = argparser.get_int("--num_iter").unwrap() as usize;
    let gap_open_penalty:f32 = argparser.get_float("--gap_open_penalty").unwrap() as f32;
    let num_threads:usize = argparser.get_int("--num_threads").unwrap() as usize;
    let normalize:bool = argparser.get_bool("--normalize").unwrap();

    let a3m_pairwise:bool = argparser.get_bool("--a3m_pairwise").unwrap();

    let tree_guided:bool = argparser.get_bool("--tree_guided").unwrap();
    let max_cluster_size:i64 = argparser.get_int("--max_cluster_size").unwrap() as i64;
    let out_matrix_path:&Option<String> = &argparser.get_string("--out_matrix");
    
    let gap_penalty_auto_adjust:bool = argparser.get_bool("--gap_penalty_auto_adjust").unwrap();
    let gap_penalty_a1_a2:Vec<String> = argparser.get_string("--gap_penalty_a1_a2").unwrap().to_string().split(',').map(|m|m.to_owned()).collect();
    
    let tree_type = if &tree_type == "NJ"{
        guide_tree_based_alignment::TreeType::TreeNj
    }else{
        assert!(&tree_type == "UPGMA");
        guide_tree_based_alignment::TreeType::TreeUPGMA
    };

    let mut rngg:StdRng = if let Some(x) = argparser.get_int("--random_seed"){
        assert!(x >  0,"--random_seed must be a positive integer. {}",x);
        StdRng::seed_from_u64(x as u64)
    }else{
        StdRng::from_entropy()
    };

    if max_cluster_size > -1 {
        assert!(max_cluster_size >  10,"--max_cluster_size must be > 10 or -1. {}",max_cluster_size);
    }// -1 の場合 hierarcical align を行わない

    
    let voidpair:Vec<(&str,&str)> = vec![
        ("--a3m_pairwise","--tree_guided"),
        ("--a3m_pairwise","--distance_base"),
        ("--num_iter","--tree_guided"),
        ("--gap_penalty_auto_adjust","--gap_open_penalty"),
        ("--gap_penalty_a1_a2","--gap_open_penalty"),
    ];
    for (a,b) in voidpair{
        if argparser.user_defined(a) && argparser.user_defined(b){
            panic!("{} and {} can not be used at the same time.",a,b);
        }    
    }

    argparser.print_items();
    
    let mut alignment_type = AlignmentType::Global;
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
    
    let infiles:Vec<String> = infiles_.split(",").into_iter().map(|m|m.to_string()).collect();
    let gmatstats:Vec<GMatStatistics>;
    let mut profile_seq:Option<SequenceProfile> = None;
    unsafe{
        gmatstats = calc_vec_stats(&infiles );// 統計値のために一回ファイルを読んでいるが後で変更する
    }

    let mut name_ordered:Vec<String> = vec![];
    let mut gmat1_:Vec<(String,Vec<char>,Vec<Vec<f32>>,Option<Vec<(f32,f32,f32,f32)>>)> = vec![];
    for ii in infiles.iter(){
        let mut z = load_multi_gmat(ii,ii.ends_with(".gz"));
        gmat1_.append(&mut z);
    }

    let veclen = gmat1_[0].2[0].len();
    let mut saligner:ProfileAligner = if gap_penalty_auto_adjust{
        ProfileAligner::new(veclen,300,None
        ,alignment_type,score_type,Some(GapPenaltyAutoAdjustParam{
            a1:gap_penalty_a1_a2[0].parse::<f32>().unwrap_or_else(|e| panic!("{:?} {:?}",gap_penalty_a1_a2,e)),
            a2:gap_penalty_a1_a2[1].parse::<f32>().unwrap_or_else(|e| panic!("{:?} {:?}",gap_penalty_a1_a2,e))
        }))
    }else{
        ProfileAligner::new(veclen,300, Some((gap_open_penalty))
        ,alignment_type,score_type,None)
    };
    
    let mut allseqs_:Vec<SequenceProfile> = vec![];
    for mut gg in gmat1_.into_iter(){
        let n = gg.0.clone();
        if name_to_res.contains_key(&n){
            panic!("Duplicated name {}.",n);
        }
        name_to_res.insert(n.clone(),"".to_owned());
        name_ordered.push(n);
        if normalize{
            gmat::normalize_seqmatrix(&mut (gg.2), &gmatstats);
        }
        let alen = gg.1.len();
        let seq = SequenceProfile::new(
            vec![(gg.0,gg.1)],alen,gg.2[0].len(),None,Some(gg.2),gg.3
        );
        allseqs_.push(seq);
    }

    
    if a3m_pairwise{
        allseqs_.reverse();
        let firstseq = allseqs_.pop().unwrap();
        
        let mut lines:Vec<String> = vec![];
        lines.push(
            ">".to_owned()+firstseq.headers[0].as_str()
        );
        lines.push(
            firstseq.member_sequences[0].iter().filter(|c| **c != '-').map(|m| m.to_string()).collect::<Vec<String>>().join("")
        );

        while allseqs_.len() > 0{
            let fst = firstseq.clone();
            let bseq = allseqs_.pop().unwrap();
            let bname = bseq.headers[0].clone();
            let seqvec = vec![fst,bseq];
            let mut ans = saligner.make_msa(seqvec,false);
            assert!(ans.len() == 1);
            let (alires,score) = ans.remove(0);

            println!(">{}",bname);
            println!("score:{}",score);

            insert_alinged_string(&alires, &mut name_to_res);

            let first_aa:Vec<char> = name_to_res.get(&name_ordered[0]).unwrap().chars().into_iter().collect();

            let pres:Vec<char> = name_to_res.get(&bname).unwrap().chars().into_iter().collect();
            assert_eq!(pres.len(),first_aa.len());
            let mut pstr:Vec<String> = vec![];
            for (aa,bb) in first_aa.iter().zip(pres.iter()){
                if aa == &'-'{
                    pstr.push(
                        bb.to_ascii_lowercase().to_string()
                    );
                }else{
                    pstr.push(
                        bb.to_string()
                    );
                }
            }
            lines.push(">".to_owned()+bname.as_str());
            lines.push(pstr.join(""));

            name_to_res.insert(firstseq.headers[0].clone(),"".to_owned());
        }
        save_lines(&outfile, lines,outfile.ends_with(".gz"));
        std::process::exit(0);   
    }

    assert!(num_threads > 0);
    rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global().unwrap();
    for ii in 0..num_iter{
        eprintln!("iter: {}",ii);

        let mut seqvec:Vec<SequenceProfile> = vec![];
        
        if let Some(p) = profile_seq{
            seqvec.push(p.create_merged_dummy());
        }
        profile_seq = None;
        for ss in allseqs_.iter(){
            seqvec.push(ss.clone());
        }
        
        let mut ans = if tree_guided{
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
            insert_alinged_string(&alires, &mut name_to_res);
            if let Some(x) = out_matrix_path{
                ioutil::save_gmat(x,&vec![(0,&alires)], x.ends_with(".gz"));
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
    save_lines(&outfile, results,outfile.ends_with(".gz"));
    //stattest();
}
