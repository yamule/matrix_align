use std::collections::*;
use matrix_align::gmat::{self, calc_vec_stats, GMatStatistics};
use matrix_align::aligner::{ScoredSeqAligner,ScoredSequence};
use matrix_align::ioutil::{load_multi_gmat, save_lines};

fn argparse(args:Vec<String>)->HashMap<String,String>{
    let mut ret:HashMap<String,String> = HashMap::new();
    let mut non_key_count:usize = 0;
    let mut ii:usize = 0;
    while(ii < args.len()){
        if args[ii].starts_with("-"){
            assert!(ii < args.len()-1,"{} does not have a value.",args[ii]);
            assert!(!ret.contains_key(&args[ii]),"{} has already been assigned.",args[ii]);
            ret.insert(args[ii].clone(),args[ii+1].clone());
            ii += 1;
        }else{
            ret.insert(format!("nokey_{}",non_key_count),args[ii].clone());
            non_key_count += 1;
        }
        ii += 1;
    }
    return ret;
}


fn main(){
    let argss: HashMap<String,String> = argparse(std::env::args().collect::<Vec<String>>());
    let required_arg:Vec<&str> = vec!["--in","--out"];
    for rr in required_arg.iter(){
        if !argss.contains_key(*rr){
            panic!("{} is required.",rr);
        }
    }
    let infile = argss.get("--in").unwrap().clone();
    let outfile = argss.get("--out").unwrap().clone();
    let num_iter:usize = argss.get("--num_iter").unwrap_or(&"2".to_owned()).parse::<usize>().unwrap_or_else(|e|panic!("in --num_iter {:?}",e));
    
    let mut name_to_res:HashMap<String,String> = HashMap::new();
    
    let gmatstats:Vec<GMatStatistics>;
    let mut profile_seq:Option<ScoredSequence> = None;
    unsafe{
        gmatstats = calc_vec_stats(& vec![infile.clone()]);// 統計値のために一回ファイルを読んでいるが後で変更する
    }
    let mut name_order:Vec<String> = vec![];
    let gmat1_ = load_multi_gmat(&infile,infile.ends_with(".gz"));
    for ii in 0..num_iter{
        eprintln!("iter: {}",ii);

        let veclen = gmat1_[0].2[0].len();
        let mut saligner:ScoredSeqAligner = ScoredSeqAligner::new(veclen,200,100);
        let mut seqvec:Vec<ScoredSequence> = vec![];
        
        if let Some(p) = profile_seq{
            seqvec.push(p);
        }
        profile_seq = None;

        let gmat1 = gmat1_.clone();
        for mut tt in gmat1.into_iter(){
            if ii == 0{
                let n = tt.0.clone();
                if name_to_res.contains_key(&n){
                    panic!("Duplicated name {}.",n);
                }
                name_to_res.insert(n.clone(),"".to_owned());
                name_order.push(n);
            }
            gmat::normalize_seqmatrix(&mut (tt.2), &gmatstats);
            let seq2 = ScoredSequence::new(
                vec![(tt.0,tt.1)],tt.2[0].len(),None,Some(tt.2),None
            );    
            seqvec.push(seq2);
        }

        let mut ans = saligner.make_msa(seqvec,-10.0,-0.5,false);
        assert!(ans.0.len() == 1);
        let alires = ans.0.pop().unwrap();
        if num_iter == ii+1{
            for (ali,hh) in alires.alignments.into_iter().zip(alires.headers.into_iter()){
                if name_to_res.contains_key(&hh){
                    name_to_res.insert(
                        hh,ali.into_iter().map(|m|m.to_string()).collect::<Vec<String>>().concat()
                    );
                }
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
