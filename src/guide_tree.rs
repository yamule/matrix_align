use crate::aligner::ScoredSequence;
use std::sync::{Mutex,Arc};
use std::thread;
use super::matrix_process::calc_euclid_dist;
use super::neighbor_joining::generate_unrooted_tree;
use super::gmat::calc_mean;
use super::aligner::ScoredSeqAligner;


pub fn create_distence_tree(val:&Vec<&Vec<f32>>)-> Vec<(i64,i64,f32)>{
    let vsiz:usize = val.len();
    let mut dist:Vec<f32> = vec![];
    for ii in 0..vsiz{
        for jj in 0..ii{
            dist.push(calc_euclid_dist(val[ii], val[jj]));
        }
        //neigbor_joining モジュールの仕様上
        dist.push(0.0);
    }
    //println!("{:?}",dist);
    return generate_unrooted_tree(&mut dist);
}

pub fn merge_with_weight(aligner:&mut ScoredSeqAligner,aseq:ScoredSequence,bseq:ScoredSequence,aweight:f32,bweight:f32)->(ScoredSequence,f32){
    let dpres = aligner.perform_dp(&aseq,&bseq,aligner.gap_open_penalty,aligner.gap_extension_penalty);
    let res = ScoredSeqAligner::make_alignment(aligner,aseq,bseq,dpres.0,false,Some((aweight,bweight)));
    return (res,dpres.1);
}

pub fn tree_guided_alignment(sequences:Vec<ScoredSequence>,aligner:&mut ScoredSeqAligner,num_threads:usize)-> (Vec<ScoredSequence>,f32){
    
    assert!(sequences.len() > 1);
    if sequences.len() == 2{
        return aligner.make_msa(sequences,false);
    }
    let mut averaged_value:Vec<Vec<f32>> = vec![];
    for ss in sequences.iter(){
        unsafe{
            let vv:Vec<&Vec<f32>> = ss.gmat.iter().map(|m|&m.match_vec).collect();
            averaged_value.push(calc_mean(&vv));
        }
    }
    
    let treenodes = create_distence_tree(&averaged_value.iter().map(|m|m).collect());
    let mut flagcounter:Vec<i64> = vec![0;treenodes.len()]; //0 は子ノードが全て計算されたもの
    let mut parents:Vec<i64> = vec![-1;treenodes.len()];
    for ii in 0..treenodes.len(){
        if treenodes[ii].0 > -1{
            if treenodes[ii].0 as usize == ii{
                assert!(treenodes[ii].1 == -1);
                continue;
            }
            assert!(parents[treenodes[ii].0 as usize] == -1);
            parents[treenodes[ii].0 as usize] = ii as i64;
            flagcounter[ii] -= 1;
        }
        if treenodes[ii].1 > -1{
            assert!(parents[treenodes[ii].1 as usize] == -1);
            parents[treenodes[ii].1 as usize] = ii as i64;
            flagcounter[ii] -= 1;
        }
    }
    let mut profiles:Vec<Option<ScoredSequence>> = vec![None;treenodes.len()];
    let numseq:usize = sequences.len();

    for ss in sequences.into_iter().enumerate(){
        profiles[ss.0] = Some(ss.1);
    }
    let mut updated:Vec<usize> = vec![];
    let leaves = flagcounter.clone();
    for ff in 0..leaves.len(){
        if leaves[ff] == 0{
            if parents[ff] > -1{
                let p = parents[ff] as usize;
                flagcounter[p] += 1;
                if flagcounter[p] == 0{
                    updated.push(p);
                }
            }
        }
    }

    let mut aligners:Vec<ScoredSeqAligner> = vec![];

    for ii in 0..num_threads{
        aligners.push(aligner.clone());
    }

    //println!("{:?}",treenodes);
    //最終的に 3 つノードが残る
    while updated.len() > 0{
        let mut handles = vec![];
        let results_:Arc<Mutex<Vec<(usize,ScoredSequence,ScoredSeqAligner)>>> = Arc::new(Mutex::new(vec![]));
        while updated.len() > 0{

            let uu = updated.pop().unwrap();
            let aa = treenodes[uu].0 as usize;
            let bb = treenodes[uu].1 as usize;
            println!("{}>>>> {} {}",uu,aa,bb);
            let adist = treenodes[aa].2;
            let bdist = treenodes[bb].2;

            profiles.push(None);
            let ap = profiles.swap_remove(aa).unwrap();
            profiles.push(None);
            let bp = profiles.swap_remove(bb).unwrap();
            let mut ali = aligners.pop().unwrap();
            let results = Arc::clone(&results_);
            let handle = thread::spawn(move || {
                let res = merge_with_weight(&mut ali,ap,bp,bdist/(adist+bdist),adist/(adist+bdist));
                results.lock().unwrap_or_else(|e| panic!("multithreaderr1 {:?}",e)).push((
                    uu,res.0,ali
                ));
            });
            handles.push(handle);
            if aligners.len() == 0{
                break;
            }
        }
        
        for handle in handles {
            handle.join().unwrap_or_else(|e| panic!("multithreadingerror ???{:?}",e));
        }
        let results = Arc::try_unwrap(results_).unwrap_or_else(
                |e| panic!("multithreaderr2")).into_inner().unwrap_or_else(|e| panic!("multithreaderr3 {:?}",e));
        for rr in results.into_iter(){
            profiles.push(Some(rr.1));
            profiles.swap_remove(rr.0);
            if parents[rr.0] > -1{
                let p = parents[rr.0] as usize;
                flagcounter[p] += 1;
                if flagcounter[p] == 0{
                    updated.push(p);
                }
            }
            aligners.push(rr.2);
        }

    }
    let mut remained:Vec<(f32,usize)> = vec![];
    for ii in 0..profiles.len(){
        if let Some(x) = & profiles[ii]{
            remained.push((treenodes[ii].2,ii));
        }
    }
    //近いペアを並べる
    remained.sort_by(|a,b| a.0.partial_cmp(&b.0).unwrap());
    let aa = remained[0].1;
    let bb = remained[1].1;
    let adist = treenodes[aa].2;
    let bdist = treenodes[bb].2;

    profiles.push(None);
    let ap = profiles.swap_remove(aa).unwrap();
    profiles.push(None);
    let bp = profiles.swap_remove(bb).unwrap();
    let res = merge_with_weight(aligner,ap,bp,bdist/(adist+bdist),adist/(adist+bdist));
    if remained.len() == 2{
        return (vec![res.0],res.1);
    }
    assert_eq!(remained.len(),3);
    profiles.push(None);
    let cp = profiles.swap_remove(remained[2].1).unwrap();
    let res = merge_with_weight(aligner,res.0,cp,0.5,0.5);
    return (vec![res.0],res.1);
}