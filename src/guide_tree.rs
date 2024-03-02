use crate::aligner::{self, ScoredSequence};
use super::matrix_process::calc_euclid_dist;
use super::neighbor_joining::generate_unrooted_tree;
use super::gmat::calc_mean;
use super::aligner::ScoredSeqAligner;
pub fn create_distence_tree(val:&Vec<&Vec<f32>>)-> Vec<(i64,i64,f32)>{
    let vsiz:usize = val.len();
    let mut dist:Vec<f32> = vec![];
    for ii in 0..vsiz{
        //neigbor_joining モジュールの仕様上
        dist.push(0.0);
        for jj in (ii+1)..vsiz{
            dist.push(calc_euclid_dist(val[ii], val[jj]));
        }
    }
    return generate_unrooted_tree(&mut dist);
}

pub fn tree_guided_alignment(sequences:Vec<ScoredSequence>,mut aligner:ScoredSeqAligner)-> ScoredSequence{
    let mut averaged_value:Vec<Vec<f32>> = vec![];
    for ss in sequences.iter(){
        unsafe{
            let vv:Vec<&Vec<f32>> = ss.gmat.iter().map(|m|&m.match_vec).collect();
            averaged_value.push(calc_mean(&vv));
        }
    }
    
    let mut treenodes = create_distence_tree(averaged_value);
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

    for ss in sequences.into_iter().enumerate(){
        profiles[ss.0] = Some(ss.1);
    }
    let mut updated:Vec<usize> = vec![];
    for ff in 0..flagcounter.len(){
        if flagcounter[ff] == 0{
            if parents[ff] > -1{
                let p = parents[ff] as usize;
                flagcounter[p] += 1;
                if flagcounter[p] == 0{
                    updated.push(p);
                }
            }
        }
    }

    //最終的に 3 つノードが残る
    while updated.len() > 0{
        // そのうち並列化する
        let uu = updated.pop().unwrap();
        
        let aa = treenodes[uu].0 as usize;
        let bb = treenodes[uu].1 as usize;

        profiles.push(None);
        let ap = profiles.swap_remove(aa).unwrap();
        profiles.push(None);
        let bp = profiles.swap_remove(bb).unwrap();
        let res = aligner.make_msa(vec![ap,bp], -10.0, -0.5, false);
        assert!(res.0.len() == 1);
        profiles.push(res.0.pop());
        profiles.swap_remove(uu);
        if parents[uu] > -1{
            let p = parents[uu] as usize;
            flagcounter[p] += 1;
            if flagcounter[p] == 0{
                updated.push(p);
            }
        }
    }
    let mut remained:Vec<(f32,usize)> = vec![];
    for ii in 0..profiles.len(){
        if let Some(x) = profiles[ii]{
            remained.push((treenodes[ii].2,ii));
        }
    }
    remained.sort_by(|a,b| a.0.partial_cmp(&b.0).unwrap());
    let aa = remained[0].1;
    let bb = remained[1].1;
    profiles.push(None);
    let ap = profiles.swap_remove(aa).unwrap();
    profiles.push(None);
    let bp = profiles.swap_remove(bb).unwrap();
    let res = aligner.make_msa(vec![ap,bp], -10.0, -0.5, false);
    if remained.len() == 2{
        return res.0[0];
    }
    assert_eq!(remained.len(),3);
    profiles.push(None);
    let cp = profiles.swap_remove(remained[2].1).unwrap();
    let res = aligner.make_msa(vec![res.0[0],cp], -10.0, -0.5, false);
    ここから
    ToDo
    ウエイト付きマージ
    penalty を ALIGNER のメンバに
    return res.0[0];
}