use crate::aligner::ScoredSequence;
use crate::neighbor_joining;
use core::num;
use std::collections::HashMap;
use std::sync::{Mutex,Arc};
use std::thread::{self, panicking};
use super::matrix_process::calc_euclid_dist;
use super::neighbor_joining::generate_unrooted_tree;
use super::gmat::calc_weighted_mean;
use super::aligner::ScoredSeqAligner;
use rand::rngs::StdRng;
use rand::Rng;
use rayon::prelude::*;
use rayon;

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

//与えられたベクトルの集合を 4 つのクラスタに分割する
pub fn soft_cluster(val:&Vec<&Vec<f32>>, rngg:&mut StdRng)->Vec<Vec<usize>>{
    assert!(val.len() > 10);
    let numseq = val.len();
    let base:usize = rngg.gen_range(0..numseq);
    let mut dist:Vec<(usize,f32)> = vec![];
    for ii in 0..numseq{
        if ii == base{
            dist.push((ii,-1.0));
            continue;
        }
        dist.push((ii,calc_euclid_dist(val[base], val[ii])));
    }
    dist.sort_by(|a,b|a.1.partial_cmp(&b.1).unwrap());

    let a1 = dist[0].0;
    let goldpoint = (numseq as f64/2.618) as usize;
    let a2 = dist[goldpoint].0;
    let a3 = dist[numseq-1-goldpoint].0;
    let a4 = dist[numseq-1].0;

    let basepoints:Vec<usize> = vec![a1,a2,a3,a4];
    
    let mut clusterids:Vec<i64> = vec![-1;numseq];
    for ii in 0..basepoints.len(){
        clusterids[basepoints[ii]] = ii as i64;
    }

    for ii in 0..numseq{
        if clusterids[ii] > -1{//basepoint
            continue;
        }
        let mut minid = 0_usize;
        let mut minval = calc_euclid_dist(val[basepoints[0]], val[ii]);//前の結果を利用すれば端折れるが、今後のことを考えて計算し直す
        for jj in 1..basepoints.len(){
            let ddis = calc_euclid_dist(val[basepoints[jj]], val[ii]);
            if ddis < minval{
                minval = ddis;
                minid = jj;
            }
        }
        clusterids[ii] = minid as i64;
    }
    let mut ret:Vec<Vec<usize>> = vec![vec![];basepoints.len()];
    for ii in 0..clusterids.len(){
        ret[clusterids[ii] as usize].push(ii);
    }
    return ret;
}

pub fn merge_with_weight(aligner:&mut ScoredSeqAligner,aseq:ScoredSequence,bseq:ScoredSequence,aweight:f32,bweight:f32)->(ScoredSequence,f32){
    let dpres = aligner.perform_dp(&aseq,&bseq,aligner.gap_open_penalty,aligner.gap_extension_penalty);
    let res = ScoredSeqAligner::make_alignment(aligner,aseq,bseq,dpres.0,false,Some((aweight,bweight)));
    return (res,dpres.1);
}

pub fn tree_guided_alignment(sequences:Vec<ScoredSequence>,aligner:&mut ScoredSeqAligner, max_cluster_size:i64, rngg:&mut StdRng,leave_one:bool)-> Vec<(ScoredSequence,f32)>{
    assert!(sequences.len() > 1);
    if sequences.len() == 2{
        return aligner.make_msa(sequences,false);
    }
    let mut averaged_value:Vec<Vec<f32>> = vec![];
    for ss in sequences.iter(){
        unsafe{
            let mut vv:Vec<&Vec<f32>> = ss.gmat.iter().map(|m|&m.match_vec).collect();
            let mut ww:Vec<f32> = ss.gmat.iter().map(|m|m.match_ratio).collect();
            assert_eq!(vv[vv.len()-1].iter().sum::<f32>(), 0.0,"???{:?}",ss.gmat[ss.gmat.len()-1]);
            vv.pop();
            ww.pop();
            averaged_value.push(calc_weighted_mean(&vv,&ww));
        }
        //ToDo 外に出す

    }
    /*
    if max_cluster_size > 0{
        if sequences.len() as i64 > max_cluster_size{
            let clusterids = soft_cluster(&averaged_value.iter().map(|m|m).collect(), rngg);


        }
        return 
    }*/


    /*
    ToDo: 最も長いエッジをとって、その両端から親子関係を作っていくべき
    そうすると最後に二つノードが残ることになると思う
    ・ルートからの三本のうちの一つである場合なのもしないで良い
    ・リーフノードである場合 Outgroup にする
    ・それ以外の場合親を探して親から逆方向に Leaf へパスを作っていく、自分からも Leaf へパスを作っていく
    */
    let mut treenodes = create_distence_tree(&averaged_value.iter().map(|m|m).collect());
    eprintln!("{:?}",treenodes);
    

    let mut node_to_seq:HashMap<usize,usize> = HashMap::new();
    let mut flagcounter:Vec<i64> = vec![0;treenodes.len()]; //0 は子ノードが全て計算されたもの
    let mut parents:Vec<i64> = vec![-1;treenodes.len()];
    let mut maxnode = 0;// 最も長い枝は最後まで残す
    let mut is_leaf = false;
    for ii in 0..treenodes.len(){
        if treenodes[ii].2 < treenodes[maxnode].2{
            maxnode = ii;
        }
        if treenodes[ii].0 > -1{
            if treenodes[ii].0 as usize == ii{
                node_to_seq.insert(ii,treenodes[ii].0 as usize);//Node (branch)には配列順に入っている
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

    let mut is_leaf = false;
    let mut noparent = false;
    if parents[maxnode] == -1{
        noparent = true;
    }
    if treenodes[maxnode].0 == -1 || treenodes[maxnode].1 == -1{
        is_leaf = true;
    }

    if !noparent {
        //元から最終ノードである場合はなにもしない
        
    }else{
        let (outree,oldmap, _) = if is_leaf{
            neighbor_joining::set_outgroup(maxnode, &treenodes,None)
        }else{
            neighbor_joining::change_center_branch(maxnode, &treenodes,None)
        };

        treenodes = outree;
        let mut newseqmap:HashMap<usize,usize> = HashMap::new();
        for ii in 0..oldmap.len(){
            if node_to_seq.contains_key(&ii){
                newseqmap.insert(oldmap[ii] as usize,*node_to_seq.get(&ii).unwrap_or_else(||panic!("???")));
            }
        }
        assert_eq!(node_to_seq.len(), newseqmap.len());
        node_to_seq = newseqmap;
    }


    
    let mut profiles:Vec<Option<ScoredSequence>> = vec![None;treenodes.len()];
    let _numseq:usize = sequences.len();

    let seq_to_node:HashMap<usize,usize> = node_to_seq.iter().map(|m|(*m.1,*m.0)).collect();
    assert_eq!(seq_to_node.len(),node_to_seq.len(),"Error in code");
    for ss in sequences.into_iter().enumerate(){
        profiles[*seq_to_node.get(&ss.0).unwrap()] = Some(ss.1);
    }
    let mut updated_pool:Vec<usize> = vec![];
    let leaves = flagcounter.clone();
    for ff in 0..leaves.len(){
        if leaves[ff] == 0{
            if parents[ff] > -1{
                let p = parents[ff] as usize;
                flagcounter[p] += 1;
                if flagcounter[p] == 0{
                    updated_pool.push(p);
                }
            }
        }
    }

    let mut aligners:Vec<ScoredSeqAligner> = vec![];
    let num_threads = rayon::current_num_threads().max(1);

    for _ in 0..num_threads{
        aligners.push(aligner.clone());
    }

    //println!("{:?}",treenodes);
    //最終的に 3 つノードが残る
    while updated_pool.len() > 0{

        let mut updated_minibatch:Vec<(usize,ScoredSequence,ScoredSequence,f32,f32,ScoredSeqAligner)> = vec![];
        while updated_pool.len() > 0{
            let idx_target = updated_pool.pop().unwrap();
            let idx_child_a = treenodes[idx_target].0 as usize;
            let idx_child_b = treenodes[idx_target].1 as usize;
            //println!("{}>>>> {} {}",uu,aa,bb);
            let adist = treenodes[idx_child_a].2;
            let bdist = treenodes[idx_child_b].2;

            profiles.push(None);
            let profile_child_a = profiles.swap_remove(idx_child_a).unwrap();
            profiles.push(None);
            let profile_child_b = profiles.swap_remove(idx_child_b).unwrap();
            let ali = aligners.pop().unwrap();
            updated_minibatch.push((idx_target,profile_child_a,profile_child_b,adist,bdist,ali));
            if aligners.len() == 0{
                break;
            }
        }
        let results:Vec<(ScoredSeqAligner,usize,ScoredSequence,f32)> = updated_minibatch.into_par_iter().map(|v|{
            let uu = v.0;
            let ap = v.1;
            let bp = v.2;
            let adist = v.3;
            let bdist = v.4;
            let mut ali = v.5;

            let res = merge_with_weight(&mut ali,ap,bp,bdist/(adist+bdist),adist/(adist+bdist));
            (ali,uu,res.0,res.1)            
        }).collect();
        
        for rr in results.into_iter(){
            profiles.push(Some(rr.2));
            profiles.swap_remove(rr.1);
            if parents[rr.1] > -1{
                let p = parents[rr.1] as usize;
                flagcounter[p] += 1;
                if flagcounter[p] == 0{
                    updated_pool.push(p);
                }
            }
            aligners.push(rr.0);
        }

    }
    let mut remained:Vec<(f32,usize)> = vec![];
    for ii in 0..profiles.len(){
        if let Some(x) = & profiles[ii]{
            remained.push((treenodes[ii].2,ii));
        }
    }

    //残りのうち近いペアを並べる
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
        return vec![(res.0,res.1)];
    }
    assert_eq!(remained.len(),3);
    profiles.push(None);
    let cp = profiles.swap_remove(remained[2].1).unwrap();
    let res = merge_with_weight(aligner,res.0,cp,0.5,0.5);
    return vec![(res.0,res.1)];
}