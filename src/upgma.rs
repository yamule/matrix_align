
use std::collections::HashMap;
use rayon;
use rayon::prelude::*;
use super::guide_tree::calc_pos;



//
pub fn get_next_pair(dist:&Vec<f32>,is_dead:&Vec<bool>)->((usize,usize),f32){
    let vlen = is_dead.len();
    let mut kmin:f32 = std::f32::MAX;
    let mut pair:(usize,usize) = (0,0);
    for ii in 0..vlen{
        if is_dead[ii]{
            continue;
        }
        for jj in ii+1..vlen{
            if is_dead[jj]{
                continue;
            }
            let ddist = dist[calc_pos(ii,jj)]; 
            if ddist < kmin {
                kmin = ddist;
                pair = (ii,jj);
            }
        }
    }
    assert!(pair.0 != pair.1);
    return (pair,kmin);
}

//calc_pos でペアの距離が得られるような distance matrix を渡し、UPGMA の系統樹を返す。最終ノード が Root Node;
fn get_upgma_tree(mut dist:Vec<f32>)->Vec<(i64,i64,f32)>{
    let leafnum:usize = ((-1.0+(1.0 as f64 +8.0*dist.len() as f64).sqrt()+0.0001) as usize)/2;
    let mut is_dead:Vec<bool> = vec![false;leafnum];
    let mut min_candidate:Vec<(f32,usize)> = vec![(std::f32::MAX,0);leafnum];
    for ii in 0..leafnum{
        for jj in (ii+1)..leafnum{
            let ppos = calc_pos(ii, jj);
            let ddist = dist[ppos];
            if min_candidate[ii].0 > ddist{
                min_candidate[ii].0 = ddist;
                min_candidate[ii].1 = jj;
            }
            if min_candidate[jj].0 > ddist{
                min_candidate[jj].0 = ddist;
                min_candidate[jj].1 = ii;
            }
        }
    }
    let mut ret:Vec<(i64,i64,f32)> = vec![];
    let mut merged_map:Vec<usize> = (0..leafnum).into_iter().collect();//現在 dist 上の index に割り当てられている node が ret 上のどの要素か
    for ii in 0..leafnum{
        ret.push((ii as i64,-1,-1.0));
    }

    loop{
        let mut nextpair = get_next_pair(&dist,&is_dead);
        //
        let mut minval = std::f32::MAX;
        let mut minpair:(usize,usize) = (0,0);
        for jj in 0..leafnum{
            if is_dead[jj]{
                continue;
            }
            if minval > min_candidate[jj].0{
                minval = min_candidate[jj].0;
                minpair = (jj,min_candidate[jj].1);
            }
        }
        if minval == std::f32::MAX{
            break;
        }
        nextpair = (minpair,minval);

        assert!((nextpair.0).0 < (nextpair.0).1);
        let a = (nextpair.0).0;
        let b = (nextpair.0).1;
        is_dead[b] = true;
        min_candidate[b].0 = std::f32::MAX;

        //最終ノードが Root Node であるはず
        ret.push(
            (merged_map[a] as i64,merged_map[b] as i64,nextpair.1)
        );
        merged_map[a] = ret.len() -1;
        for jj in 0..leafnum{
            if is_dead[jj]{
                continue;
            }
            if jj == a || jj == b{
                continue;
            }
            let ppos1 = calc_pos(a, jj);
            let ppos2 = calc_pos(b, jj);
            let ddist = (dist[ppos1]+dist[ppos2])/2.0;
            dist[ppos1] = ddist;
            if min_candidate[jj].1 == a || min_candidate[jj].1 == b{
                let mut minval = std::f32::MAX;
                let mut minidx = 0_usize;
                for kk in 0..leafnum{
                    if jj == kk || is_dead[kk]{
                        continue;
                    }
                    let qpos = calc_pos(kk, jj);
                    if dist[qpos] > minval{
                        minval = dist[qpos];
                        minidx = kk;
                    }
                }
                min_candidate[jj].0 = minval;
                min_candidate[jj].1 = minidx;
            }
        }
        let mut minval = std::f32::MAX;
        let mut minidx = 0_usize;
        for jj in 0..leafnum{
            if is_dead[jj]{
                continue;
            }
            if a == jj{
                continue;
            }
            let qpos = calc_pos(a, jj);
            if dist[qpos] < minval{
                minval = dist[qpos];
                minidx = jj;
            }
        }
        min_candidate[a].0 = minval;
        min_candidate[a].1 = minidx;
    }
    //id、親から与えられた距離
    let mut updater:Vec<(usize,f32)> = vec![];
    let mut currentindex = ret.len()-1;
    let mut currentval = ret[currentindex].2;
    updater.push((ret[currentindex].0 as usize,currentval/2.0));
    updater.push((ret[currentindex].1 as usize,currentval/2.0));
    ret[currentindex].2 = 0.0;
    while updater.len() > 0{
        let p = updater.pop().unwrap();
        currentindex = p.0;
        if ret[currentindex].1 == -1{
            ret[currentindex].2 = currentval;
            continue;
        }
        assert!(ret[currentindex].0 != -1);
        assert!(ret[currentindex].2 != std::f32::MAX);

        currentval = ret[currentindex].2;
        let pval = p.1 - currentval/2.0;
        ret[currentindex].2 = pval;
        updater.push((ret[currentindex].0 as usize,currentval/2.0));
        updater.push((ret[currentindex].1 as usize,currentval/2.0));
    }
    return ret;
}
