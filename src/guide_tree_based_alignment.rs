
use super::aligner::SequenceProfile;
use super::faster_neighbor_joining;
use super::upgma;
use std::collections::HashMap;
use super::matrix_process::calc_euclid_dist;
use super::gmat::calc_weighted_mean;
use super::aligner::ProfileAligner;
use rayon::prelude::*;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use super::bisecting_kmeans;


#[derive(Copy,Clone)]
pub enum TreeType{
    TreeNj,TreeUPGMA
}
pub fn create_distence_tree(val:&Vec<&Vec<f32>>,num_threads:usize,typ:TreeType)-> Vec<(i64,i64,f32)>{
    let vsiz:usize = val.len();
    let mut dist:Vec<f32> = vec![];
    println!("calc_distance...");
    for ii in 0..vsiz{
        for jj in 0..ii{
            dist.push(calc_euclid_dist(val[ii], val[jj]));
        }
        //neigbor_joining モジュールの仕様上
        dist.push(0.0);
    }
    match typ{
        TreeType::TreeNj =>{
            println!("calc_nj_tree...");
            //println!("{:?}",dist);
            return faster_neighbor_joining::generate_unrooted_tree(&mut dist,num_threads);
        },
        TreeType::TreeUPGMA =>{
            println!("calc_upgma_tree...");
            //println!("{:?}",dist);
            let ret =  upgma::generate_unrooted_tree(&mut dist);
            return ret;
        },
    }
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

pub fn align_and_merge_with_weight(aligner:&mut ProfileAligner,aseq:SequenceProfile,bseq:SequenceProfile,aweight:f32,bweight:f32)->(SequenceProfile,f32){
    let dpres = aligner.perform_dp(&aseq,&bseq);
    let res = ProfileAligner::make_alignment(aligner,aseq,bseq,dpres.0,false,Some((aweight,bweight)));
    return (res,dpres.1);
}

pub fn calc_distance_from(anchors:&Vec<usize>,val:&Vec<Vec<f32>>) -> Vec<Vec<f32>>{
    let mut ret:Vec<Vec<f32>> = vec![];
    let mut anchorpoint:Vec<&Vec<f32>> = vec![];
    for ii in 0..anchors.len(){
        anchorpoint.push(&val[ii]);
    }
    for vv in 0..val.len(){
        let mut xv:Vec<f32> = vec![];
        for ii in 0..anchorpoint.len(){
            let d = calc_euclid_dist(&val[vv],anchorpoint[ii]);
            xv.push(d);
        }
        ret.push(xv);
    }
    return ret;
}

pub fn hierarchical_alignment(sequences:Vec<SequenceProfile>,aligner:&mut ProfileAligner, max_cluster_member_size:i64, rngg:&mut StdRng,num_threads_:usize,tree_type:TreeType)-> Vec<(SequenceProfile,f32)>{
    assert!(sequences.len() > 1);

    if sequences.len() as i64 <= max_cluster_member_size {
        if sequences.len() == 1{
            //なんか mut にするのが気持ち悪いので
            return sequences.into_iter().map(|m| (m,-1.0)).collect();
        }
        return tree_guided_alignment(sequences, aligner, false,num_threads_,tree_type);
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
    }

    let num_threads = num_threads_/20+1;
    let num_mini_threads = num_threads_/num_threads;

    let num_anchors = 5;
    let num_trial = 10;// 何回 anchor の選び直しをするか。
    let mut trial_buff:Vec<Vec<usize>>= vec![];
    let mut minloss = -1.0;
    let mut mincluster:Vec<Vec<usize>> = vec![];
    for _ in 0..num_trial{
        let mut mbuff:Vec<usize> = vec![];
        for _ in 0..num_anchors{
            mbuff.push(rngg.gen_range(0..averaged_value.len()));
        } 
        trial_buff.push(
            mbuff
        );
    }

    while trial_buff.len() > 0{
        let mut trial_minibatch:Vec<(Vec<Vec<f32>>,u64)> = vec![];
        for _ in 0..num_threads{
            let p = trial_buff.pop().unwrap();
            trial_minibatch.push((
                calc_distance_from(&p,&averaged_value),
                rngg.gen::<u64>())
            );
            if trial_buff.len() == 0{
                break;
            }
        }
        println!("clustering minibatch: {}",trial_minibatch.len());
        let res:Vec<Result<(Vec<Vec<usize>>,f32),()>> = trial_minibatch.into_par_iter().map(|m|
        bisecting_kmeans::bisecting_kmeans(&m.0.iter().collect(),max_cluster_member_size as usize,&mut StdRng::seed_from_u64(m.1),num_mini_threads)
        ).collect();
        for rr in res.into_iter(){
            if let Ok(x) = rr{
                let vres = x.0;
                let lloss = x.1;
                if minloss < 0.0 || minloss > lloss{
                    mincluster = vres;
                    minloss = lloss;
                }
            }
        }
        println!("cluster_loss: {}",minloss);
    }
    let mut sequences_hm:HashMap<usize,SequenceProfile> = sequences.into_iter().enumerate().collect();
    let mut clustered:Vec<Vec<SequenceProfile>> = vec![];
    for cluster in mincluster.into_iter(){
        let mut cc:Vec<SequenceProfile> = vec![];
        for clustermem in cluster.into_iter(){
            if !sequences_hm.contains_key(&clustermem){
                panic!("Error in code ???");
            }
            cc.push(sequences_hm.remove(&clustermem).unwrap());
        }
        println!("clustersize: {}",cc.len());
        clustered.push(cc);
    }
    assert!(sequences_hm.len() == 0);
    let mut subres:Vec<SequenceProfile> = vec![];
    for mut cc in clustered.into_iter(){
        if cc.len() > 3{
            let alires = tree_guided_alignment(cc, aligner, true,num_threads_,tree_type);
            for aa in alires.into_iter(){
                subres.push(aa.0);
            }
        }else{
            subres.append(&mut cc);
        }
    }

    if subres.len() as i64 > max_cluster_member_size{
        return hierarchical_alignment(subres, aligner, max_cluster_member_size, rngg, num_threads_,tree_type);
    }
    return tree_guided_alignment(subres, aligner, false,num_threads_,tree_type);

}

pub fn tree_guided_alignment(sequences:Vec<SequenceProfile>,aligner:&mut ProfileAligner,skip_last_merge:bool,num_threads:usize,tree_type:TreeType)-> Vec<(SequenceProfile,f32)>{
    assert!(sequences.len() > 1);
    if sequences.len() == 2{
        return aligner.make_msa(sequences,false);
    }
    println!("align {} sequences.",sequences.len());
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
    }
    /*
    if max_cluster_size > 0{
        if sequences.len() as i64 > max_cluster_size{
            let clusterids = soft_cluster(&averaged_value.iter().map(|m|m).collect(), rngg);


        }
        return 
    }*/


    /*
    最も長いエッジをとって、その両端から親子関係を作っていく
    ・ルートからの三本のうちの一つである場合なのもしない
    ・リーフノードである場合 Outgroup にする
    ・それ以外の場合親を探して親から逆方向に Leaf へパスを作っていく、自分の子は自分へのパスを削除する（最後まで残る）
    */
    let mut treenodes;
    let use_upgma = match tree_type{TreeType::TreeUPGMA=>{true},TreeType::TreeNj=>{false}};
    if use_upgma{
        treenodes = create_distence_tree(&averaged_value.iter().map(|m|m).collect(),num_threads,TreeType::TreeUPGMA);
    }else{
        
        treenodes = create_distence_tree(&averaged_value.iter().map(|m|m).collect(),num_threads,TreeType::TreeNj);
    }
    //println!("{:?}",treenodes);
    


    let mut node_to_seq:HashMap<usize,usize> = HashMap::new();
    let parents:Vec<i64> = faster_neighbor_joining::get_parent_branch(&treenodes);

    if use_upgma{
        for ii in 0..treenodes.len(){
            if treenodes[ii].0 > -1{
                if treenodes[ii].0 as usize == ii{
                    node_to_seq.insert(ii,treenodes[ii].0 as usize);//Node (branch)には配列順に入っている
                    assert!(treenodes[ii].1 == -1);
                    continue;
                }
            }
        }
    }else{
        let mut maxnode = 0;// 最も長い枝は最後まで残す
        let mut is_leaf = false;

        for ii in 0..treenodes.len(){
            if treenodes[ii].2 > treenodes[maxnode].2{
                maxnode = ii;
            }
            if treenodes[ii].0 > -1{
                if treenodes[ii].0 as usize == ii{
                    node_to_seq.insert(ii,treenodes[ii].0 as usize);//Node (branch)には配列順に入っている
                    is_leaf = true;
                    assert!(treenodes[ii].1 == -1);
                    continue;
                }
            }
        }

        //if sequences.len() <= 40{
            //println!("{}",maxnode);
            //println!("{:?}",treenodes);
            //println!("{:?}",parents);
        //}

        let mut noparent = false;
        if parents[maxnode] == -1{
            noparent = true;
        }
        if treenodes[maxnode].0 == -1 || treenodes[maxnode].1 == -1{
            is_leaf = true;
        }

        if noparent {
            //元から最終ノードである場合はなにもしない      
            //println!("nochange");
        }else{
            let (outree,oldmap, _) = if is_leaf{
                //println!("leaf");
                faster_neighbor_joining::set_outgroup(maxnode, &treenodes,None)
            }else{
                //println!("internal");
                faster_neighbor_joining::change_center_branch(maxnode, &treenodes,None)
            };

            treenodes = outree;
            let mut newseqmap:HashMap<usize,usize> = HashMap::new();
            for ii in 0..oldmap.len(){
                if node_to_seq.contains_key(&ii){
                    newseqmap.insert(oldmap[ii] as usize,*node_to_seq.get(&ii).unwrap_or_else(||panic!("???")));
                }
            }
            assert_eq!(node_to_seq.len(), newseqmap.len());
            //println!("{:?} {:?} {:?}",oldmap,node_to_seq,newseqmap);
            node_to_seq = newseqmap;
        }
    }

    let mut flagcounter:Vec<i64> = vec![0;treenodes.len()]; //0 は子ノードが全て計算されたもの
    let parents:Vec<i64> = faster_neighbor_joining::get_parent_branch(&treenodes);
    for ii in 0..treenodes.len(){
        if treenodes[ii].0 > -1{
            if treenodes[ii].0 == ii as i64{
                assert!(treenodes[ii].1 == -1);
                continue;
            }else{
                flagcounter[ii] -= 1;
            }
        }
        if treenodes[ii].1 > -1{
            flagcounter[ii] -= 1;
        }
    }
    
    let mut profiles:Vec<Option<(SequenceProfile,f32)>> = vec![None;treenodes.len()];
    let _numseq:usize = sequences.len();

    let seq_to_node:HashMap<usize,usize> = node_to_seq.iter().map(|m|(*m.1,*m.0)).collect();
    assert_eq!(seq_to_node.len(),node_to_seq.len(),"Error in code");
    for ss in sequences.into_iter().enumerate(){
        profiles[*seq_to_node.get(&ss.0).unwrap()] = Some((ss.1,-1.0));
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

    let mut aligners:Vec<ProfileAligner> = vec![];

    for _ in 0..num_threads{
        aligners.push(aligner.clone());
    }
    println!("align");
    //println!("{:?}",treenodes);
    //unrooted の場合最終的に 3 つノードが残る
    while updated_pool.len() > 0{

        let mut updated_minibatch:Vec<(usize,SequenceProfile,SequenceProfile,f32,f32,ProfileAligner)> = vec![];
        while updated_pool.len() > 0{
            let idx_target = updated_pool.pop().unwrap();
            if let Some(x) = &profiles[idx_target]{
                assert!(idx_target == 0); //root は merge しない
                continue;
            }
            let idx_child_a = treenodes[idx_target].0 as usize;
            let idx_child_b = treenodes[idx_target].1 as usize;
            //println!("{}<< {} {}",idx_target,idx_child_a,idx_child_b);
            let adist = treenodes[idx_child_a].2;
            let bdist = treenodes[idx_child_b].2;

            profiles.push(None);
            let profile_child_a = profiles.swap_remove(idx_child_a).unwrap();
            profiles.push(None);
            let profile_child_b = profiles.swap_remove(idx_child_b).unwrap();
            let ali = aligners.pop().unwrap();
            updated_minibatch.push((idx_target,profile_child_a.0,profile_child_b.0,adist,bdist,ali));
            if aligners.len() == 0{
                break;
            }
        }
        //println!("multi_align_minibatch:{}",updated_minibatch.len());
        //rayon による並列処理
        let results:Vec<(ProfileAligner,usize,SequenceProfile,f32)> = updated_minibatch.into_par_iter().map(|v|{
            let uu = v.0;
            let ap = v.1;
            let bp = v.2;
            let adist = v.3;
            let bdist = v.4;
            let mut ali = v.5;

            let res = align_and_merge_with_weight(&mut ali,ap,bp,bdist/(adist+bdist),adist/(adist+bdist));
            (ali,uu,res.0,res.1)            
        }).collect();
        
        for rr in results.into_iter(){
            profiles.push(Some((rr.2,rr.3)));
            let pst = profiles.swap_remove(rr.1);
            if let Some(x) = pst{
                panic!("This is not expected ????");// None であるはず
            }
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
    if remained.len() == 1{//rooted tree
        let xx = profiles.swap_remove(remained[0].1);
        if let Some(x) =  xx{
            return vec![x];
        }
        panic!("???");
    }
    //残りのうち近いペアを並べる
    remained.sort_by(|a,b| a.0.partial_cmp(&b.0).unwrap());
    let aa = remained[0].1;
    let bb = remained[1].1;
    let adist = treenodes[aa].2;
    let bdist = treenodes[bb].2;

    profiles.push(None);
    let ap = profiles.swap_remove(aa).unwrap().0;
    profiles.push(None);
    let bp = profiles.swap_remove(bb).unwrap().0;
    let res = align_and_merge_with_weight(aligner,ap,bp,bdist/(adist+bdist),adist/(adist+bdist));
    if remained.len() == 2{
        return vec![(res.0,res.1)];
    }
    
    assert_eq!(remained.len(),3);
    profiles.push(None);
    let cp = profiles.swap_remove(remained[2].1).unwrap();
    
    if skip_last_merge{
        return vec![(res.0,res.1),cp];
    }

    let res = align_and_merge_with_weight(aligner,res.0,cp.0,0.5,0.5);
    return vec![res];
}


