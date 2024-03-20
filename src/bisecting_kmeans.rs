use std::collections::HashSet;

use rand::rngs::StdRng;
use rand::Rng;
use super::matrix_process::*;

static  CLUSTER_LEFT:i8 = 0;
static  CLUSTER_RIGHT:i8 = 1;
static  CLUSTER_NONE:i8 = -1;


pub fn bisecting_kmeans(val:&Vec<&Vec<f32>>,max_member_num:usize,rng:&mut StdRng) -> Result<(Vec<Vec<usize>>,f32),()>{
    let mut sample_bag:Vec<(Vec<usize>,f32)> = vec![((0..val.len()).into_iter().collect(),1000.0)];
    let mut cluster_mapping_buff:Vec<i8> = vec![0;val.len()];
    let mut final_clusters:Vec<(Vec<usize>,f32)> = vec![];
    while sample_bag.len() > 0{
        let (idd,currentloss) = sample_bag.pop().unwrap();
        if idd.len() <= max_member_num{
            final_clusters.push((idd,currentloss));
            continue;
        }
        let mut minloss:f32 = -100.0;
        let mut minloss_lr:(f32,f32) =  (1000.0,1000.0);
        let mut minright:HashSet<usize> = HashSet::new();
        for _ii in 0..4{
            let res = split(val,&idd,rng,&mut cluster_mapping_buff);
            if let Ok((lloss,rloss)) = res{
                let loss = lloss+rloss;
                if minloss < -1.0 || minloss > loss{
                    minloss = loss;
                    minloss_lr = (lloss,rloss);
                    minright.clear();
                    for ss in idd.iter(){
                        if cluster_mapping_buff[*ss] == CLUSTER_RIGHT{
                            minright.insert(*ss);
                        }
                    }
                }
            }
        }
        if minloss > -1.0{
            let mut sample_left:Vec<usize> = vec![];
            for ss in idd.iter(){
                if !minright.contains(ss){
                    sample_left.push(*ss);
                }
            }
            sample_bag.push((sample_left,minloss_lr.0));
            sample_bag.push((minright.into_iter().collect(),minloss_lr.1));
        }else{
            return Err(());
        }
    }
    let mut res:Vec<Vec<usize>> = vec![];
    let mut total_loss:f32 = 0.0;
    for ff in final_clusters.into_iter(){
        res.push(ff.0);
        total_loss += ff.1;
    }
    return Ok((res,total_loss));
}

//dist_to_anchor と名前が付いているが、Kalign のアルゴリズム上の話であり、単なる座標と考えて OK。
pub fn split(dist_to_anchor:&Vec<&Vec<f32>>,sample_ids:&Vec<usize>,rng:&mut StdRng,cluster_result:&mut Vec<i8>)-> Result<(f32,f32),()>{
    let _num_seqs:usize = dist_to_anchor.len();
    let num_anchors:usize = dist_to_anchor[0].len();
    let num_samples:usize = sample_ids.len();
    let mut center:Vec<f32> = vec![0.0;num_anchors];

    for ss in sample_ids.iter(){
        for jj in 0..num_anchors{
            center[jj] += dist_to_anchor[*ss][jj];
        }
    }

    element_multiply(&mut center, 1.0/(num_samples as f32));

    let firstsameple = sample_ids[rng.gen_range(0..num_samples)];
    
    let mut leftcenter:Vec<f32> = dist_to_anchor[firstsameple].iter().map(|m| *m).collect();
    let mut rightcenter:Vec<f32> = vec![0.0;num_anchors];
    let mut has_diff:bool = false;
    for jj in 0..num_anchors{//中心を挟んで反対側の点を取る
        if (leftcenter[jj] - center[jj]).abs() > 1.0e-12{
            has_diff = true;
        }
        rightcenter[jj] = center[jj]*2.0 - leftcenter[jj];
    }

    for ss in sample_ids.iter(){//一応初期化するが基本的に Err が返ってきたら使用しない
        cluster_result[*ss] = CLUSTER_NONE;
    }

    if !has_diff{
        return Err(());
    }

    let mut lloss = 1000000000.0_f32;
    let mut rloss = 1000000000.0_f32;
    let mut loopcount = 0_i64;
    loop{
        loopcount += 1;
        if loopcount > 100000{
            for ss in sample_ids.iter(){
                cluster_result[*ss] = CLUSTER_NONE;
            }
            return Err(());
        }
        let mut left_current:Vec<f32> = vec![0.0;num_anchors];
        let mut right_current:Vec<f32> = vec![0.0;num_anchors];
        let mut num_left:usize = 0;
        let mut num_right:usize = 0;
        lloss = 0.0;
        rloss = 0.0;

        for ss in sample_ids.iter(){
            let ldist = calc_euclid_dist(dist_to_anchor[*ss], &leftcenter);
            let rdist = calc_euclid_dist(dist_to_anchor[*ss], &rightcenter);
            if ldist > rdist{
                vector_add(&mut right_current,dist_to_anchor[*ss]);
                rloss += rdist;
                num_right += 1;
                cluster_result[*ss] = CLUSTER_RIGHT;
            }else{
                vector_add(&mut left_current,dist_to_anchor[*ss]);
                lloss += ldist;
                num_left += 1;
                cluster_result[*ss] = CLUSTER_LEFT;
            }
        }
        if num_left == 0 || num_right == 0{
            for ss in sample_ids.iter(){
                cluster_result[*ss] = CLUSTER_NONE;
            }
            return Err(());
        }

        element_multiply(&mut left_current, 1.0/(num_left as f32));
        element_multiply(&mut right_current, 1.0/(num_right as f32));
        if left_current == leftcenter && right_current == rightcenter{
            break;
        }
        leftcenter = left_current;
        rightcenter = right_current;
    }

    return Ok((lloss,rloss));
}
