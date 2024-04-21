use std::collections::*;
use rayon;
use matrix_align::gmat::{self, calc_vec_stats, GMatStatistics};
use matrix_align::aligner::{AlignmentType, GapPenaltyAutoAdjustParam, ProfileAligner, ScoreType, SequenceProfile};
use matrix_align::ioutil::{self, load_multi_gmat, save_lines};
use matrix_align::guide_tree_based_alignment::{self, DistanceBase};
use rand::SeedableRng;
use rand::rngs::StdRng;

fn argparse(mut args:Vec<String>)->HashMap<String,String>{
    assert!(args.len() > 0);
    let fst = args.remove(0);
    assert!(fst.contains("matrix_align"));
    let mut ret:HashMap<String,String> = HashMap::new();
    let mut non_key_count:usize = 0;
    let mut ii:usize = 0;
    while ii < args.len(){
        if args[ii].starts_with("--"){
            assert!(!ret.contains_key(&args[ii]),"{} has already been assigned.",args[ii]);
            if ii == args.len()-1 || args[ii+1].starts_with("--"){
                ret.insert(args[ii].clone(),"true".to_owned());
            }else{
                ret.insert(args[ii].clone(),args[ii+1].clone());
                ii += 1;
            }
        }else{
            ret.insert(format!("--nokey_{}",non_key_count),args[ii].clone());
            non_key_count += 1;
        }
        ii += 1;
    }
    return ret;
}

fn check_bool(a:&str,argkey:&str)-> bool{
    if a.to_lowercase() == "true" || a == "1" {
        return true;
    }
    if a.to_lowercase() == "false" || a == "0" {
        return false;
    }
    panic!("Boolean value must be true, false, 1, or 0. {}:{}",argkey,a);
}

fn check_acceptable(a:&str,acceptable:Vec<&str>)-> bool{ // 10 もないだろうので
    for vv in acceptable.iter(){
        if *vv == a{
            return true;
        }
    }
    panic!("{} is not an acceptable value. {:?}",a,acceptable);
}

fn main(){
    main_(std::env::args().collect::<Vec<String>>());
}
fn main_(args:Vec<String>){
    let argss: HashMap<String,String> = argparse(args);
    let required_arg:Vec<&str> = vec!["--in","--out"];
    let allowed_arg_:Vec<(&str,&str)> = vec![
        ("--in","<input file path> : Text based file which contains general profile matrices for multiple sequences. Something like as follows. \n>seq1\n[amino acid letter at position 1]\t[value1]\t[value2]\t[value3]...\n[amino acid letter at position 2]\t[value1]\t[value2]\t[value3]...\n[amino acid letter at position 3]\t[value1]\t[value2]\t[value3]...\n...\n>seq2\n[amino acid letter at position 1]\t[value1]\t[value2]\t[value3]...\n[amino acid letter at position 2]\t[value1]\t[value2]\t[value3]...\n[amino acid letter at position 3]\t[value1]\t[value2]\t[value3]...\n...\n..."),
        ("--out","<output file path> : Multi-fasta file. Required."),
        ("--num_iter","<int> : Number of alignment iterations. Construct global profile with -1 times with this number. Default '2'."),
        ("--gap_open_penalty","<float> : Gap open penalty for DP. Must be negative. Default '-10.0'"),
        ("--gap_extension_penalty","<float> : Gap extension penalty for DP. Must be negative. Default '-0.5'."),

        ("--gap_penalty_auto_adjust","<bool or novalue=true> : Adjust gap penalty automaticall. <gap open penalty> = <maximum matching score>*a1*-1.0 + <minimum matching score>*a2; <gap extension penalty> = <gap open penalty>*0.05; Default true. "),
        ("--gap_penalty_a1_a2","<float>,<float> : Parameters for gap penalty auto adjust. Default 0.5,0.5."),

        ("--normalize","<bool or novalue=true> : Normalize profile values before alignment. Default 'false'."),
        ("--alignment_type","<global or local> : Alignment type. Default 'global'"),
        ("--num_threads","<int> : Number of threads. Default 4."),
        ("--tree_guided","<bool or novalue=true> : Use guide-tree based alignment. Default true."),
        ("--distance_base","<string> : Estimation procedure for distance between samples. \"averaged_value\" or \"score_with_seed\". Default score_with_seed."),

        ("--max_cluster_size","<int> : Limit all-vs-all comparison with this many number of profiles and align hierarchically. Must be > 10. Default -1."),
        ("--random_seed","<int> : Seed for random number generator."),
        ("--tree_type","<string> : Type of guide tree. \"NJ\" or \"UPGMA\". Default \"NJ\"."),
        ("--score_type","<string> : Type of scoring function. \"distance_zscore\" or \"dot_product\". Default \"distance_zscore\"."),
        ("--out_matrix","<string> : Save matrix file of the resulting alignment."),
        ("--help","<bool or novalue=true> : Print this message."),
    ];
    let allowed_arg:HashSet<&str> = allowed_arg_.clone().into_iter().map(|m|m.0).collect();

    let mut error_message:Vec<String> = vec![];
    for rr in required_arg.iter(){
        if !argss.contains_key(*rr){
            error_message.push(format!("{} is required.",rr));
        }
    }
    let mut argcheck:Vec<String> = vec![];
    for aa in argss.iter(){
        if !allowed_arg.contains(&(aa.0).as_str()){
            if aa.0.starts_with("--nokey"){
                argcheck.push(format!("{}",aa.1));
            }else{
                argcheck.push(format!("{} {}",aa.0,aa.1));
            }
        }
    }
    if argcheck.len() > 0{
        argcheck.sort();
        error_message.push(format!("Unknown arg: {:?}",argcheck));
    }
    
    if !argss.contains_key("--in") || !argss.contains_key("--out"){
        panic!("--in and --out is required.");
    }

    let tree_type = argss.get("--tree_type").unwrap_or(&"NJ".to_owned()).clone();
    let distance_base = argss.get("--distance_base").unwrap_or(&"score_with_seed".to_owned()).clone();
    
    let infile = argss.get("--in").unwrap().clone();
    let outfile = argss.get("--out").unwrap().clone();
    let num_iter:usize = argss.get("--num_iter").unwrap_or(&"2".to_owned()).parse::<usize>().unwrap_or_else(|e|panic!("in --num_iter {:?}",e));
    let gap_open_penalty:f32 = argss.get("--gap_open_penalty").unwrap_or(&"-10.0".to_owned()).parse::<f32>().unwrap_or_else(|e|panic!("in --gap_open_penalty {:?}",e));
    let gap_extension_penalty:f32 = argss.get("--gap_extension_penalty").unwrap_or(&"-0.5".to_owned()).parse::<f32>().unwrap_or_else(|e|panic!("in --gap_extension_penalty {:?}",e));
    let num_threads:usize = argss.get("--num_threads").unwrap_or(&"4".to_owned()).parse::<usize>().unwrap_or_else(|e|panic!("in --num_threads {:?}",e));
    let normalize:bool = check_bool(argss.get("--normalize").unwrap_or(&"false".to_owned()).as_str(),"--normalize");
    let print_help:bool = check_bool(argss.get("--help").unwrap_or(&"false".to_owned()).as_str(),"--help");
    let tree_guided:bool = check_bool(argss.get("--tree_guided").unwrap_or(&"true".to_owned()).as_str(),"--tree_guided");
    let max_cluster_size:i64 = argss.get("--max_cluster_size").unwrap_or(&"-1".to_owned()).parse::<i64>().unwrap_or_else(|e|panic!("in --maximum_cluster_size {:?}",e));
    let out_matrix_path:&Option<&String> = &argss.get("--out_matrix");
    
    let gap_penalty_auto_adjust:bool = check_bool(argss.get("--gap_penalty_auto_adjust").unwrap_or(&"true".to_owned()).as_str(),"--gap_penalty_auto_adjust");
    let gap_penalty_a1_a2:Vec<String> = argss.get("--gap_penalty_a1_a2").unwrap_or(&"0.5,0.5".to_owned()).to_string().split(',').map(|m|m.to_owned()).collect();
    
    if gap_penalty_auto_adjust{
        assert!(!argss.contains_key("--gap_open_penalty") && !argss.contains_key("--gap_extension_penalty"), "To use gap_open_penalty or gap_extension_penalty, set gap_penalty_auto_adjust \"false\".");
    }

    check_acceptable(&tree_type.as_str(),vec!["NJ","UPGMA"] );
    let tree_type = if &tree_type == "NJ"{
        guide_tree_based_alignment::TreeType::TreeNj
    }else{
        guide_tree_based_alignment::TreeType::TreeUPGMA
    };


    if print_help{
        for aa in allowed_arg_.iter(){
            eprintln!("{} {}",aa.0,aa.1);
        }
        std::process::exit(0);
    }

    let mut rngg:StdRng = if let Some(x) = argss.get("--random_seed"){
        let i = x.parse::<i64>().unwrap_or(-1);
        assert!(i >  0,"--random_seed must be a positive integer. {}",x);
        StdRng::seed_from_u64(i as u64)
    }else{
        StdRng::from_entropy()
    };// -1 の場合シードを使わない

    if max_cluster_size > -1 {
        assert!(max_cluster_size >  10,"--max_cluster_size must be > 10 or -1. {}",max_cluster_size);
    }// -1 の場合 hierarcical align を行わない


    let mut alignment_type = AlignmentType::Global;
    if argss.contains_key("--alignment_type"){
        let typ = argss.get("--alignment_type").unwrap().to_lowercase();
        check_acceptable(typ.as_str(),vec!["global","local"] );
        if typ.as_str() == "global"{
            alignment_type = AlignmentType::Global;
        }else if typ.as_str() == "local"{
            alignment_type = AlignmentType::Local;
        }else{
            panic!("???");
        }
    }
    
    let mut distance_base = DistanceBase::ScoreWithSeed;
    if argss.contains_key("--distance_base"){
        let typ = argss.get("--distance_base").unwrap().to_lowercase();
        check_acceptable(typ.as_str(),vec!["score_with_seed","averaged_value"] );
        if typ.as_str() == "score_with_seed"{
            distance_base = DistanceBase::ScoreWithSeed;
        }else if typ.as_str() == "averaged_value"{
            distance_base = DistanceBase::AveragedValue;
        }else{
            panic!("???");
        }
    }
    
    let mut score_type = ScoreType::DistanceZscore;
    if argss.contains_key("--score_type"){
        let typ = argss.get("--score_type").unwrap().to_lowercase();
        check_acceptable(typ.as_str(),vec!["distance_zscore","dot_product"] );
        if typ.as_str() == "distance_zscore"{
            score_type = ScoreType::DistanceZscore;
        }else if typ.as_str() == "dot_product"{
            score_type = ScoreType::DotProduct;
        }else{
            panic!("???");
        }
    }

    if error_message.len() > 0{
        for aa in allowed_arg_.iter(){
            eprintln!("{} {}",aa.0,aa.1);
        }
        eprintln!("Error:");
        for aa in error_message.iter(){
            eprintln!("{}",aa);
        }
        
        panic!();
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
                guide_tree_based_alignment::hierarchical_alignment(seqvec, &distance_base,&mut saligner,max_cluster_size, &mut rngg,num_threads,tree_type)
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
