use std::collections::*;
use matrix_align::simple_argparse;
use rayon;
use matrix_align::gmat::{self, calc_vec_stats, GMatStatistics};
use matrix_align::aligner::{AlignmentType, GapPenaltyAutoAdjustParam, ProfileAligner, ScoreType, SequenceProfile};
use matrix_align::ioutil::{self, load_multi_gmat, save_lines};
use matrix_align::guide_tree_based_alignment::{self, DistanceBase};
use rand::SeedableRng;
use rand::rngs::StdRng;










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
        ..."
        ,None,vec![], true),
        
        ("--out",None
        ,"<output file path> : Multi-fasta file. Required."
        ,None,vec![], true),
        
        ("--num_iter",None
        ,"<int> : Number of alignment iterations. Construct global profile with -1 times with this number."
        ,Some("2"),vec![],false),
        
        ("--gap_open_penalty",None
        ,"<float> : Gap open penalty for DP. Must be negative."
        ,Some("-10.0"),vec![],false),
        
        ("--gap_extension_penalty",None
        ,"<float> : Gap extension penalty for DP. Must be negative."
        ,Some("-0.5"),vec![],false),

        ("--gap_penalty_auto_adjust",None
        ,"<bool or novalue=true> : Adjust gap penalty automaticall. <gap open penalty> = <maximum matching score>*a1*-1.0 + \
        <minimum matching score>*a2; <gap extension penalty> = <gap open penalty>*0.05;"
        ,Some("true"),vec![],false),
        
        ("--gap_penalty_a1_a2",None
        ,"<float>,<float> : Parameters for gap penalty auto adjust."
        ,Some("0.5,0.5"),vec![],false),

        ("--normalize",None
        ,"<bool or novalue=true> : Normalize profile values before alignment."
        ,Some("false"),vec![],false),
        
        ("--alignment_type",None
        ,"<string> : Alignment type."
        ,Some("global"),vec!["global","local"],false),
        
        ("--num_threads",None
        ,"<int> : Number of threads."
        ,Some("4"),vec![],false),
        
        ("--tree_guided",None
        ,"<bool or novalue=true> : Use guide-tree based alignment."
        ,Some("true"),vec![],false),
        
        ("--distance_base",None
        ,"<string> : Estimation method for distance between samples."
        ,Some("score_with_seed"),vec!["averaged_value","score_with_seed"],false),

        ("--max_cluster_size",None
        ,"<int> : Limit all-vs-all comparison with this many number of profiles and align hierarchically. Must be > 10."
        ,None,vec![],false),
        
        ("--random_seed",None
        ,"<int> : Seed for random number generator. Must be positive."
        ,None,vec![],false),
        
        ("--tree_type",None
        ,"<string> : Type of guide tree."
        ,Some("NJ"),vec!["NJ","UPGMA"],false),
        
        ("--score_type",None
        ,"<string> : Type of scoring function. \"distance_zscore\" or \"dot_product\"."
        ,Some("distance_zscore"),vec![],false),
        
        ("--out_matrix",None
        ,"<string> : Save matrix file of the resulting alignment."
        ,None,vec![],false)
    ];

    let mut argparser = simple_argparse::SimpleArgParse::new(allowed_arg_);
    argparser.parse(args);

    let tree_type = argparser.get_string("--tree_type").unwrap();
    
    let infile = argparser.get_string("--in").unwrap().clone();
    let outfile = argparser.get_string("--out").unwrap().clone();
    let num_iter:usize = argparser.get_int("--num_iter").unwrap() as usize;
    let gap_open_penalty:f32 = argparser.get_float("--gap_open_penalty").unwrap() as f32;
    let gap_extension_penalty:f32 = argparser.get_float("--gap_extension_penalty").unwrap() as f32;
    let num_threads:usize = argparser.get_int("--num_threads").unwrap() as usize;
    let normalize:bool = argparser.get_bool("--normalize").unwrap();

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
    
    let gmatstats:Vec<GMatStatistics>;
    let mut profile_seq:Option<SequenceProfile> = None;
    unsafe{
        gmatstats = calc_vec_stats(& vec![infile.clone()]);// 統計値のために一回ファイルを読んでいるが後で変更する
    }

    let mut name_order:Vec<String> = vec![];
    let gmat1_ = load_multi_gmat(&infile,infile.ends_with(".gz"));

    let veclen = gmat1_[0].2[0].len();
    let mut saligner:ProfileAligner = 
    
    if gap_penalty_auto_adjust{
        ProfileAligner::new(veclen,300,None
        ,alignment_type,score_type,Some(GapPenaltyAutoAdjustParam{
            a1:gap_penalty_a1_a2[0].parse::<f32>().unwrap_or_else(|e| panic!("{:?} {:?}",gap_penalty_a1_a2,e)),
            a2:gap_penalty_a1_a2[1].parse::<f32>().unwrap_or_else(|e| panic!("{:?} {:?}",gap_penalty_a1_a2,e))
        }))
    }else{
        ProfileAligner::new(veclen,300, Some((gap_open_penalty,gap_extension_penalty))
        ,alignment_type,score_type,None)
    };
    
    let mut allseqs_:Vec<SequenceProfile> = vec![];
    for mut gg in gmat1_.into_iter(){
        let n = gg.0.clone();
        if name_to_res.contains_key(&n){
            panic!("Duplicated name {}.",n);
        }
        name_to_res.insert(n.clone(),"".to_owned());
        name_order.push(n);
        if normalize{
            gmat::normalize_seqmatrix(&mut (gg.2), &gmatstats);
        }
        //駄目っぽい。保留。
        //let ss_biased = gmat::ssbias(&mut tt.2,false);
        //tt.2 = ss_biased;
        
        let alen = gg.1.len();
        // ギャップまだ入って無い
        let seq = SequenceProfile::new(
            vec![(gg.0,gg.1)],alen,gg.2[0].len(),None,Some(gg.2),None
        );
        allseqs_.push(seq);
    }

    let similarity_sort = false;
    
    if similarity_sort{
        //先頭配列に近い順でアラインメントする
        let seq1 = allseqs_.swap_remove(0);
        let mut scoresort:Vec<(f32,SequenceProfile)> = vec![];
        for seq2 in allseqs_.into_iter(){
            let dpres = saligner.perform_dp(&seq1,&seq2);
            //let nscore = dpres.1/(seq1.get_alignment_length().min(seq2.get_alignment_length())) as f32;
            //needleman wunsh の合計だと normalize する必要はない気がするが・・・
            //うーん
            //scoresort.push((nscore,seq2));
            scoresort.push((dpres.score,seq2));
        }
        scoresort.sort_by(|a,b|a.0.partial_cmp(&b.0).unwrap());
        scoresort.reverse();
        allseqs_ = vec![seq1];
        for ss in scoresort.into_iter(){
            eprintln!("{} {:?}",ss.0,ss.1.headers);
            allseqs_.push(ss.1);
        }
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
        eprintln!("score:{}",ans[0].1);
        
        let (alires,_alisc) = ans.pop().unwrap();
        
        //let pstr = alires.gmat_str();
        //for pp in pstr{ //profile 表示
        //    println!("{}",pp);
        //}

        if num_iter == ii+1{
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

            if let Some(x) = out_matrix_path{
                ioutil::save_gmat(x,&vec![(0,&alires)], x.ends_with(".gz"));
            }
            break;
        }
        profile_seq = Some(alires);
    }
    let mut results:Vec<String> = vec![];
    for ii in name_order.iter(){
        //println!(">{}",primary_id_to_name.get(ii).unwrap());
        //println!("{}",primary_id_to_res.get(ii).unwrap());
        results.push(
            format!(">{}\n{}",ii,name_to_res.get(ii).unwrap())
        );
    }
    save_lines(&outfile, results,outfile.ends_with(".gz"));
    //stattest();
}
