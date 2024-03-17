use rand::rngs::StdRng;
use rand::Rng;
use super::matrix_process::*;

static  CLUSTER_LEFT:i8 = 0;
static  CLUSTER_RIGHT:i8 = 1;
static  CLUSTER_NONE:i8 = -1;


pub fn bisecting_kmeans(val:&Vec<&Vec<f32>>,num_clusters:usize){

}

pub fn split(dist_to_anchor:&Vec<&Vec<f32>>,sample_ids:&Vec<usize>,rng:&mut StdRng,cluster_result:&mut Vec<i8>)-> Result<f32,()>{
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

    if !has_diff{
        for ss in sample_ids.iter(){//一応初期化するが基本的に Err が返ってきたら使用しない
            cluster_result[*ss] = CLUSTER_NONE;
        }
        return Err(());
    }

    let mut loss = 1000000000.0_f32;
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
        loss = 0.0;
        for ss in sample_ids.iter(){
            let ldist = calc_euclid_dist(dist_to_anchor[*ss], &leftcenter);
            let rdist = calc_euclid_dist(dist_to_anchor[*ss], &rightcenter);
            if ldist > rdist{
                vector_add(&mut right_current,dist_to_anchor[*ss]);
                loss += rdist;
                num_right += 1;
                cluster_result[*ss] = CLUSTER_RIGHT;
            }else{
                vector_add(&mut left_current,dist_to_anchor[*ss]);
                loss += ldist;
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
    return Ok(loss);
}