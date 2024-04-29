use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon;
use rayon::prelude::*;
use super::matrix_process::*;

static  CLUSTER_LEFT:i8 = 0;
static  CLUSTER_RIGHT:i8 = 1;
static  CLUSTER_NONE:i8 = -1;

//座標にあたるベクトル、クラスタ内の最大メンバ数、乱数生成器を与えて、クラスタに含まれる ID の Vec とスコアを返す
pub fn bisecting_kmeans(val:&Vec<&Vec<f32>>,max_member_num:usize,rng:&mut StdRng, num_threads_:usize) -> Result<(Vec<Vec<usize>>,f32),()>{
    let mut sample_bag:Vec<(Vec<usize>,f32)> = vec![((0..val.len()).into_iter().collect(),1000.0)];
    let mut cluster_mapping_buff:Vec<(Vec<i8>,u64)> = vec![];
    let mut cluster_mapping_best:(Vec<i8>,(f32,f32)) = (vec![-1;val.len()],(-1.0,-1.0));
    assert!(num_threads_ > 0);
    
    let mut num_threads = 1;
    for ii in vec![40,20,10,8,4,2,1]{ //40 回クラスタ中心候補を検討するので割り切れる数にする
        if ii <= num_threads_{
            num_threads = ii;
            break;
        } 
    }
    
    for _ in 0..num_threads{
        cluster_mapping_buff.push((vec![0;val.len()],rng.gen::<u64>()));
    }
    let num_trial = 40/num_threads;

    assert!(num_trial >= 1,"Error in code." );
    let mut final_clusters:Vec<(Vec<usize>,f32)> = vec![];
    while sample_bag.len() > 0{
        let (idd,currentloss) = sample_bag.pop().unwrap();
        if idd.len() <= max_member_num{
            final_clusters.push((idd,currentloss));
            continue;
        }
        cluster_mapping_best.1 = (-1.0,-1.0);
        for _ in 0..num_trial{
            //rayon による並列処理
            let results:Vec<Result<(f32,f32,Vec<i8>),Vec<i8>>> = cluster_mapping_buff.into_par_iter().map(|v|{
                split(val,&idd,&mut StdRng::seed_from_u64(v.1),v.0)
            }).collect();
            let mut minloss_lr:(f32,f32) = (-1.0,-1.0);
            let mut minloss = -1.0;
            cluster_mapping_buff = vec![];
            let mut minidx:i32 = -1;
            for (eii,res) in results.into_iter().enumerate(){
                match res{
                    Ok((lloss,rloss,buff)) =>{
                        let loss = lloss+rloss;
                        if minloss < 0.0 || minloss > loss{
                            minloss = loss;
                            minloss_lr = (lloss,rloss);
                            minidx = eii as i32;
                        }
                        cluster_mapping_buff.push((buff,rng.gen::<u64>()));
                    },
                    Err(x) => {
                        cluster_mapping_buff.push((x,rng.gen::<u64>()));
                    }
                }
            }

            let currentbest = (cluster_mapping_best.1).0+(cluster_mapping_best.1).1;
            if minidx > -1 && (currentbest < 0.0 || minloss < currentbest){
                cluster_mapping_buff.push((cluster_mapping_best.0,rng.gen::<u64>()));
                let b = cluster_mapping_buff.swap_remove(minidx as usize);
                cluster_mapping_best = (b.0,minloss_lr);
                
                /* //デバッグ用
                let mut lcount = 0;
                let mut rcount = 0;
                
                for ii in idd.iter(){
                    if cluster_mapping_best.0[*ii] == CLUSTER_LEFT{
                        lcount += 1;
                    }
                    if cluster_mapping_best.0[*ii] == CLUSTER_RIGHT{
                        rcount += 1;
                    }
                }
                //println!("split_tmp:{} {}",lcount,rcount);
                */
            }
        }
        if (cluster_mapping_best.1).0 >= 0.0{
            let mut sample_left:Vec<usize> = vec![];
            let mut sample_right:Vec<usize> = vec![];
            for ss in idd.iter(){
                if cluster_mapping_best.0[*ss] == CLUSTER_LEFT{
                    sample_left.push(*ss);
                }else if  cluster_mapping_best.0[*ss] == CLUSTER_RIGHT{
                    sample_right.push(*ss);
                }else{
                    panic!("Error in code.");
                }
            }
            //println!("split: {} {} {:?}",sample_left.len(),sample_right.len(),cluster_mapping_best.1);
            sample_bag.push((sample_left,(cluster_mapping_best.1).0));
            sample_bag.push((sample_right,(cluster_mapping_best.1).1));
        }else{
            eprintln!("Can't find any split in {} trial.",num_threads*num_trial);
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
pub fn split(dist_to_anchor:&Vec<&Vec<f32>>,sample_ids:&Vec<usize>,rng:&mut StdRng,mut cluster_result:Vec<i8>)-> Result<(f32,f32,Vec<i8>),Vec<i8>>{
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
        return Err(cluster_result);
    }

    let mut lloss;
    let mut rloss;
    let mut loopcount = 0_i64;
    loop{
        loopcount += 1;
        if loopcount > 100000{
            for ss in sample_ids.iter(){
                cluster_result[*ss] = CLUSTER_NONE;
            }
            eprintln!("Split was not converged.");
            return Err(cluster_result);
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
            return Err(cluster_result);
        }
        
        element_multiply(&mut left_current, 1.0/(num_left as f32));
        element_multiply(&mut right_current, 1.0/(num_right as f32));
    
        if left_current == leftcenter && right_current == rightcenter{
            break;
        }
        leftcenter = left_current;
        rightcenter = right_current;
    }

    return Ok((lloss,rloss,cluster_result));
}

#[test]
fn kmeantest(){
    use rand::SeedableRng;
    use std::collections::HashSet;
    let clustercenter:Vec<Vec<f32>> = vec![
        vec![0.0,0.0,0.0],
        vec![100.0,0.0,100.0],
        vec![-100.0,-100.0,0.0],
        vec![50.0,0.0,-100.0],
    ];
    
    let mut rng = StdRng::from_entropy();
    let mut points:Vec<Vec<f32>> = vec![];
    let mut ans:Vec<usize> = vec![];
    for ii in 0..clustercenter.len(){
        for _jj in 0..100{
            let a:f32 = rng.gen_range(0.0..1.0)*20.0-10.0+clustercenter[ii][0];
            let b:f32 = rng.gen_range(0.0..1.0)*20.0-10.0+clustercenter[ii][1];
            let c:f32 = rng.gen_range(0.0..1.0)*20.0-10.0+clustercenter[ii][2];
            points.push(vec![a,b,c]);
            ans.push(ii);
        }
    }
    let val:Vec<&Vec<f32>> = points.iter().collect(); 
    let res_ = bisecting_kmeans(&val,110,&mut rng, 4);
    let mut res = res_.unwrap();
    let mut processed:HashSet<usize> = HashSet::new();
    for rr in res.0.iter_mut(){
        let cid = ans[rr[0]];
        for rrr in rr.iter(){
            assert_eq!(ans[*rrr],cid,"{} in {} vs {} in {}",rr[0],cid,*rrr,ans[*rrr]);
            processed.insert(*rrr);
            //eprintln!("{:?}",val[*rrr]);
        }
        rr.sort();
        eprintln!("{:?}",rr);
    }
    assert_eq!(processed.len(),points.len());
}