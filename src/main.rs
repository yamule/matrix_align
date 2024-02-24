use std::collections::*;
use matrix_align::gmat::{self, calc_vec_stats, calc_vec_stats_, GMatStatistics};
use matrix_align::aligner::{ScoredSeqAligner,ScoredSequence};
use matrix_align::ioutil::{load_multi_gmat, save_lines};

fn argparse(args:Vec<String>)->HashMap<String,String>{
    let mut ret:HashMap<String,String> = HashMap::new();
    let mut non_key_count:usize = 0;
    let mut ii:usize = 0;
    while(ii < args.len()){
        if args[ii].starts_with("--"){
            assert!(ii < args.len()-1,"{} does not have a value.",args[ii]);
            assert!(!ret.contains_key(&args[ii]),"{} has already been assigned.",args[ii]);
            if ii == args.len()-1 || args[ii+1].starts_with("--"){
                ret.insert(args[ii].clone(),"true".to_owned());
            }else{
                ret.insert(args[ii].clone(),args[ii+1].clone());
                ii += 1;
            }
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
    let gap_open_penalty:f32 = argss.get("--gap_open_penalty").unwrap_or(&"-10.0".to_owned()).parse::<f32>().unwrap_or_else(|e|panic!("in --gap_open_penalty {:?}",e));
    let gap_extension_penalty:f32 = argss.get("--gap_extension_penalty").unwrap_or(&"-0.5".to_owned()).parse::<f32>().unwrap_or_else(|e|panic!("in --gap_extension_penalty {:?}",e));
    let normalize:bool = argss.contains_key("--normalize");
    
    let mut name_to_res:HashMap<String,String> = HashMap::new();
    
    let mut gmatstats:Vec<GMatStatistics>;
    let mut profile_seq:Option<ScoredSequence> = None;
    unsafe{
        gmatstats = calc_vec_stats(& vec![infile.clone()]);// 統計値のために一回ファイルを読んでいるが後で変更する
    }

    let mut name_order:Vec<String> = vec![];
    let mut gmat1_ = load_multi_gmat(&infile,infile.ends_with(".gz"));

    let mut column_selector:Option<Vec<usize>> = None;    
    if true{
        let mut diffsum:Vec<f32> = vec![0.0;gmatstats.len()];
        for gg in gmat1_.iter(){
            let vval:Vec<Vec<f32>> = gg.2.clone();
            unsafe{
                let stats = calc_vec_stats_(vval);
                for ss in stats.iter().enumerate(){
                    diffsum[ss.0] += (gmatstats[ss.0].mean-ss.1.mean).abs();
                }
            }
        }
        let mut sumsort:Vec<(f32,usize)> = diffsum.into_iter().enumerate().map(|m|(m.1,m.0)).collect();
        sumsort.sort_by(|a,b|a.0.partial_cmp(&b.0).unwrap());
        println!("{:?}",sumsort);
        //sumsort.reverse();
        let mut sel:Vec<usize> = vec![];
        for ii in 0..sumsort.len()/10{
            sel.push(sumsort[ii].1);
        }

        column_selector = Some(sel);
    }


    if let Some(x) = column_selector{
        
        let hs:HashSet<usize> = x.iter().map(|x|*x).collect();
        let mut newgmat:Vec<GMatStatistics> = vec![];
        for gg in gmatstats.into_iter().enumerate(){
            if !hs.contains(&gg.0){
                newgmat.push(gg.1);
            }
        }
        gmatstats = newgmat;

        for gg in gmat1_.iter_mut(){
            for ff in gg.2.iter_mut(){
                let lnum = ff.len()-hs.len();
                let mut newg:Vec<f32> = vec![0.0;lnum];
                let mut counter  = 0;
                for vv in ff.iter().enumerate(){
                    if !hs.contains(&vv.0){
                        newg[counter] = *vv.1;
                        counter += 1;
                    }
                }
                *ff = newg;
            }
        }           
    }
    let veclen = gmat1_[0].2[0].len();
    let mut saligner:ScoredSeqAligner = ScoredSeqAligner::new(veclen,300,100);
    for ii in 0..num_iter{
        eprintln!("iter: {}",ii);

        let mut seqvec:Vec<ScoredSequence> = vec![];
        
        if let Some(p) = profile_seq{
            seqvec.push(p.create_merged_msa());
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
                if normalize{
                    gmat::normalize_seqmatrix(&mut (tt.2), &gmatstats);
                }
                //駄目っぽい。保留。
                //let ss_biased = gmat::ssbias(&mut tt.2,false);
                //tt.2 = ss_biased;
            }
            let seq2 = ScoredSequence::new(
                vec![(tt.0,tt.1)],tt.2[0].len(),None,Some(tt.2),None
            );    
            seqvec.push(seq2);
        }

        let mut ans = saligner.make_msa(seqvec,gap_open_penalty, gap_extension_penalty,false);
        assert!(ans.0.len() == 1);
        let alires = ans.0.pop().unwrap();

        //let pstr = alires.gmat_str();
        //for pp in pstr{ //profile 表示
        //    println!("{}",pp);
        //}

        if num_iter == ii+1{
            for (ali,hh) in alires.alignments.into_iter().zip(alires.headers.into_iter()){
                if name_to_res.contains_key(&hh){
                    assert!(name_to_res.get(&hh).unwrap().len() == 0,"{}",name_to_res.get(&hh).unwrap());
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
