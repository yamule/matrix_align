use crate::aligner::ScoredSequence;
use super::matrix_process::calc_euclid_dist;
use super::neighbor_joining::generate_unrooted_tree;
use super::gmat::calc_mean;

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

pub fn tree_guided_alignment(sequences:Vec<ScoredSequence>){
    let mut averaged_value:Vec<Vec<f32>> = vec![];
    for ss in sequences.iter(){
        unsafe{
            let vv:Vec<&Vec<f32>> = ss.gmat.iter().map(|m|&m.match_vec).collect();
            averaged_value.push(calc_mean(&vv));
        }
    }
    let mut treenodes = create_distence_tree(averaged_value);
    let mut parents:Vec<i64> = vec![-1;treenodes.len()];
    for ii in 0..treenodes.len(){
        if treenodes[ii].0 > -1{
            assert!(parents[treenodes[ii].0 as usize] == -1);
            parents[treenodes[ii].0 as usize] = ii as i64;
        }
        if treenodes[ii].1 > -1{
            assert!(parents[treenodes[ii].1 as usize] == -1);
            parents[treenodes[ii].1 as usize] = ii as i64;
        }
    }

    let mut profiles:Vec<Option<ScoredSequence>> = vec![None;treenodes.len()];

    for ss in sequences.into_iter().enumerate(){
        profiles[ss.0] = Some(ss.1);
    }
    let mut updated:Vec<usize> = vec![];
    let mut flagcounter:Vec<i64> = vec![0;treenodes.len()]; //子ノードが全て計算されたもの
    //最終的に 3 つノードが残る
    
    

}