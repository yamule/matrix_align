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
    
    let mut primary_id_to_name:HashMap<i32,String> = HashMap::new();
    let mut id_order:Vec<i32> = vec![];
    let mut primary_id_to_res:HashMap<i32,String> = HashMap::new();
    
    let gmatstats:Vec<GMatStatistics>;
    let mut profile_seq:Option<(Vec<char>,Vec<Vec<f32>>,Vec<(f32,f32,f32,f32)>)> = None;
    unsafe{
        gmatstats = calc_vec_stats(& vec![infile.clone()]);
    }
    
    let gmat1_ = load_multi_gmat(&infile,infile.ends_with(".gz"));
    for ii in 0..num_iter{
        eprintln!("{}",ii);

        let veclen = gmat1_[0].2[0].len();
        let mut saligner:ScoredSeqAligner = ScoredSeqAligner::new(veclen,200,100);
        let mut seqvec:Vec<ScoredSequence> = vec![];
        let mut head_id:i32 = -1;
        if let Some(p) = profile_seq{
            let dummy_center = ScoredSequence::new(
                p.0,p.1, Some(p.2),&mut saligner,true
            );
            saligner.weights[dummy_center.primary_ids[0] as usize] = gmat1_.len() as f32;
        }
        profile_seq = None;

        let gmat1 = gmat1_.clone();
        for mut tt in gmat1.into_iter(){
            gmat::normalize_seqmatrix(&mut (tt.2), &gmatstats);
            let seq2 = ScoredSequence::new(
                tt.1,tt.2,None,&mut saligner,true
            );
            if num_iter-1 == ii{
                primary_id_to_name.insert(seq2.primary_ids[0],tt.0);
                id_order.push(seq2.primary_ids[0]);    
            }
            seqvec.push(seq2);
        }
        head_id = seqvec[0].primary_ids[0];

        let ans = saligner.make_msa(seqvec,-10.0,-0.5,false);
        let alires = &ans.0[0];
        if num_iter == ii+1{
            for (eii,aa) in alires.alibuff_idx.iter().enumerate(){
                let mut rres:Vec<char> = vec![];
                for ii in 0..alires.alignment_length{
                    rres.push(saligner.ali_get(*aa,ii))
                }
                primary_id_to_res.insert(
                    alires.primary_ids[eii],rres.into_iter().map(|m|m.to_string()).collect::<Vec<String>>().concat()
                );
            }
            break;
        }
        let mut center_char:Vec<char> = vec![];
        let mut center_gmat:Vec<Vec<f32>> = vec![];
        let mut center_gaps:Vec<(f32,f32,f32,f32)> = vec![];
        for (eii,_aa) in alires.alibuff_idx.iter().enumerate(){
            if alires.primary_ids[eii] == head_id{
                center_char = vec![];
                center_gmat = vec![];
                for _ in 0..alires.alignment_length{
                    //center_char.push(saligner.ali_get(*aa,ii));
                    center_char.push('X');
                }
                for ii in 0..=alires.alignment_length{
                    let ppp = saligner.gmat_colget(alires.gmatbuff_id  as usize,ii).clone();
                    center_gaps.push((ppp.match_ratio,ppp.del_ratio,ppp.connected_ratio,ppp.gapped_ratio));
                    if ii < alires.alignment_length{
                        center_gmat.push(ppp.match_vec);
                    }
                }
            }
        }

        profile_seq = Some((
            center_char,center_gmat,center_gaps
        ));
    }
    let mut results:Vec<String> = vec![];
    for ii in id_order.iter(){
        //println!(">{}",primary_id_to_name.get(ii).unwrap());
        //println!("{}",primary_id_to_res.get(ii).unwrap());
        results.push(
            format!(">{}",primary_id_to_name.get(ii).unwrap())
        );
        results.push(
            format!("{}",primary_id_to_res.get(ii).unwrap())
        );
    }
    save_lines(&outfile, results,outfile.ends_with(".gz"));
    //stattest();
}
