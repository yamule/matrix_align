
use super::aligner::SequenceProfile;
use std::collections::HashMap;
use super::aligner::ProfileAligner;
use rayon::prelude::*;
use rand::rngs::StdRng;
use rand::Rng;
use std::collections::HashSet;

pub fn calc_alignment_based_distance_from_seed(sequences:&Vec<SequenceProfile>,aligner:&mut ProfileAligner,
    num_seed_seq_:usize,num_threads:usize,rngg:&mut StdRng
)-> Vec<Vec<f32>>{
    let num_seed_seq = num_seed_seq_.min(sequences.len());
    let num_seqs = sequences.len();

//println!("{:?}",treenodes);
    
    let mut aligners:Vec<ProfileAligner> = vec![];

    for _ in 0..num_threads{
        aligners.push(aligner.clone());
    }

    let mut seedseqs:HashSet<usize> = HashSet::new();
    loop{
        seedseqs.insert(rngg.gen_range(0..num_seqs));
         if seedseqs.len() >= num_seed_seq{
            break;
         }
    }
    let seedseqs:Vec<usize> = seedseqs.into_iter().collect();
    
    let mut seed_to_ret_index:HashMap<usize,usize> = HashMap::new();
    for ii in 0..seedseqs.len(){
        seed_to_ret_index.insert(seedseqs[ii],ii);
    }


    let mut pair_pool:Vec<(usize,usize)> = vec![];

    let mut ret:Vec<Vec<f32>> = vec![vec![0.0;num_seed_seq];num_seqs];
    for ii in seedseqs.iter(){
        for jj in 0..num_seqs{
            pair_pool.push((*ii,jj));
        }
    }
    
    while pair_pool.len() > 0{
        let mut pair_minibatch:Vec<(usize,usize,ProfileAligner)> = vec![];
        while pair_pool.len() > 0{
            let idx_target = pair_pool.pop().unwrap();
            let ali = aligners.pop().unwrap();
            pair_minibatch.push((idx_target.0,idx_target.1,ali));
            if aligners.len() == 0{
                break;
            }
        }
        //println!("multi_align_minibatch:{}",updated_minibatch.len());
        //rayon による並列処理
        let results:Vec<(usize,usize,ProfileAligner,f32)> = pair_minibatch.into_par_iter().map(|v|{
            let s1 = v.0;
            let s2 = v.1;
            let mut ali = v.2;
            let shorter_length = sequences[s1].get_alignment_length().min(sequences[s2].get_alignment_length()) as f32;
            let _longer_length = sequences[s1].get_alignment_length().max(sequences[s2].get_alignment_length()) as f32;
            let res = ali.perform_dp(&sequences[s1],&sequences[s2]);
            let mut positivecount = 0_f32;
            for ss in res.match_scores.iter(){
                if *ss > 0.0{
                    positivecount += *ss;
                }
            }
            (s1,s2,ali,positivecount/(shorter_length as f32))
        }).collect();
        
        for rr in results.into_iter(){
            let seedidx = seed_to_ret_index.get(&rr.0).unwrap_or_else(|| panic!("error in code."));
            ret[rr.1][*seedidx] = rr.3;
            aligners.push(rr.2);
        }

    }

    return ret;
}


