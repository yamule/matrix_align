use std::collections::HashMap;

use self::{matrix_process::*, misc::line_to_hash};
use super::*;

#[derive(Debug, PartialEq)]
pub struct GMatStatistics{
    pub max:f32,
    pub min:f32,
    pub sum:f32,
    pub mean:f32,
    pub var:f32,
    pub count:usize
}

impl GMatStatistics{
    pub fn get_string(&self)->String{
        return format!(
            "max:\t{}\tmin:\t{}\tsum:\t{}\tmean:\t{}\tvar:\t{}\tcount:\t{}"
            ,self.max,self.min,self.sum,self.mean,self.var,self.count
        );
    }
    pub fn load_string(strr:&str) -> GMatStatistics{
        let hs:HashMap<String,String> = misc::line_to_hash(strr);
        for kk in vec!["max","min","sum","mean","var","count"]{
            if !hs.contains_key(kk){
                panic!("The line must field {} .",kk);
            }
        }
        let mmax:f32 = hs.get("max").unwrap().parse::<f32>().unwrap();
        let mmin:f32 = hs.get("min").unwrap().parse::<f32>().unwrap();
        let ssum:f32 = hs.get("sum").unwrap().parse::<f32>().unwrap();
        let mmean:f32 = hs.get("mean").unwrap().parse::<f32>().unwrap();
        let vvar:f32 = hs.get("var").unwrap().parse::<f32>().unwrap();
        let ccount:usize = hs.get("count").unwrap().parse::<usize>().unwrap();
        return GMatStatistics{
           max:mmax,
           min:mmin,
           sum:ssum,
           mean:mmean,
           var:vvar,
           count:ccount
        }
    }
}
//calc_vec_stats_ と一部重複しているが、これだけ使うこともあるので
pub unsafe fn calc_mean(allval:&Vec<&Vec<f32>>)->Vec<f32>{
    let mut ssum:Vec<f32> = vec![];
    if ssum.len() == 0{
        ssum = vec![0.0;allval[0].len()];
    }
    for pp in allval.iter(){
        matrix_process::vector_add(&mut ssum, pp);
    }
    matrix_process::element_multiply(&mut ssum,1.0/(allval.len() as f32));
    return ssum;
}

pub unsafe fn calc_weighted_mean(allval:&Vec<&Vec<f32>>,allweight:&Vec<f32>)->Vec<f32>{
    assert_eq!(allval.len(),allweight.len());
    let mut ssum:Vec<f32> = vec![];
    if ssum.len() == 0{
        ssum = vec![0.0;allval[0].len()];
    }
    let mut wsum = 0.0_f32;
    for pp in allval.iter().zip(allweight.iter()){
        let mut vp = (*pp.0).clone();
        wsum += *pp.1;
        matrix_process::element_multiply(&mut vp, *pp.1);
        matrix_process::vector_add(&mut ssum, & vp);
    }
    assert!(wsum > 0.0);
    matrix_process::element_multiply(&mut ssum,1.0/wsum);
    return ssum;
}

pub unsafe fn calc_vec_stats_(allval:Vec<Vec<f32>>)->Vec<GMatStatistics>{
    let mut ssum:Vec<f32> = vec![];
    let mut smax:Vec<f32> = vec![];
    let mut smin:Vec<f32> = vec![];
    let mut counter = 0_usize;
    
    if ssum.len() == 0{
        ssum = vec![0.0;allval[0].len()];
        smax = allval[0].clone();
        smin = allval[0].clone();
    }
    for pp in allval.iter(){
        matrix_process::vector_add(&mut ssum, &pp);
        matrix_process::vector_max(&mut smax, &pp);
        matrix_process::vector_min(&mut smin, &pp);
        counter += 1;
    }

    let mut smean = ssum.clone();
    matrix_process::element_multiply(&mut smean,1.0/(counter as f32));

    let mut svar:Vec<f32> = vec![0.0;ssum.len()];
    let mut smeanneg = smean.clone();
    matrix_process::element_multiply(&mut smeanneg,-1.0);
    for mut pp in allval.into_iter(){
        matrix_process::vector_add(&mut pp,&smeanneg);
        matrix_process::vector_square(&mut pp);
        matrix_process::vector_add(&mut svar,&pp);
    }
    for pp in svar.iter_mut(){
        *pp /= counter as f32;
    }
    let mut ret:Vec<GMatStatistics> = vec![];
    for ii in 0..ssum.len(){
        ret.push(
            GMatStatistics{
                max:smax[ii],min:smin[ii],sum:ssum[ii],mean:smean[ii],var:svar[ii],count:counter
            }
        );
        //println!("mean:\t{}\tmax:\t{}\tmin:\t{}\tvar:{}",smean[ii],smax[ii],smin[ii],svar[ii]/(counter as f32));
    }
    return ret;
}

//各カラムの Max とか Min とか計算して返す
pub unsafe fn calc_vec_stats(filenames:&Vec<String>)->Vec<GMatStatistics>{
    let mut allval:Vec<Vec<f32>> = vec![];//まあ多分メモリ上に乘るだろう。。。
    for fname in filenames.into_iter(){
        let gmat_ = ioutil::load_multi_gmat(fname,fname.ends_with(".gz"));
        for mut gmat1 in gmat_.into_iter(){
            allval.append(&mut gmat1.2);
        }
    }
    return calc_vec_stats_(allval);
}


pub fn normalize(vec:&mut Vec<f32>,gmatstats:&Vec<GMatStatistics>){
    assert_eq!(vec.len(),gmatstats.len());
    for vv in vec.iter_mut().zip(gmatstats.iter()){
        assert!(vv.1.var > 0.0);
        *vv.0 = (*vv.0 -vv.1.mean)/vv.1.var;
    }
}

pub fn normalize_seqmatrix(vec:&mut Vec<Vec<f32>>, gmatstats:&Vec<GMatStatistics>){
    let vlen = vec.len();
    for ii in 0..vlen{
        normalize(&mut vec[ii], gmatstats);
    }
}


//secondary structure を表現する場所があると数残基シフトの Alternative Alignment が出来てしまうので調整する
pub fn ssbias(vec:&mut Vec<Vec<f32>>,ignore_last:bool) -> Vec<Vec<f32>>{
    let mut vlen = vec.len();
    if ignore_last {
        vlen -= 1;
    }
    let vec_size = vec[0].len();
    let mut ret = vec.clone();
    for ii in 0..vlen{
        for kk in 0..vec_size{// 効果があるなら SIMD にする
            let mut ccount = 0;
            let mut ssum = 0_f32;
            
            for jj in -4_i32..=4{
                if jj.abs() < 2{
                    continue;
                }
                let pos = ii as i32+jj;
                if pos < 0{
                    continue;
                }
                if pos as usize > vlen-1{
                    continue;
                }
                ssum += vec[pos as usize][kk];
                ccount += 1;
            }
            ret[ii][kk] = vec[ii][kk]-ssum/(ccount as f32)*0.5;
        }
    }
    return ret;
}


// Match State で加算される、スコアの行列を、(avec.len(),bvec.len()) の長さの行列として返す。
// 要素は一方の行列のある行と他方の行列の全行のユークリッド距離を計算し、統計値を取る。
// その統計値をもとに ZSCORE を計算し
// avec, bvec 両方で同様のことを行って平均をとって要素とする。
//Pantolini, Lorenzo, et al. "Embedding-based alignment: combining protein language models with dynamic programming alignment to detect structural similarities in the twilight-zone." Bioinformatics 40.1 (2024): btad786.
pub fn calc_dist_zscore_matrix(avec:& Vec<&Vec<f32>>,bvec:& Vec<&Vec<f32>>,aweight:Option<&Vec<f32>>,bweight:Option<&Vec<f32>>)-> Vec<Vec<f32>>{
    let alen = avec.len();
    let blen = bvec.len();
    let mut ret:Vec<Vec<f32>> = vec![];
    let vecsiz = avec[0].len() as f32;
    for rr in 0..alen{
        let mut evec:Vec<f32> = vec![];
        for cc in 0..blen{
            let distt = calc_euclid_dist(avec[rr], bvec[cc])/vecsiz;
            evec.push((-1.0*distt).exp());
        }
        let sstat:VectorStats = if let Some(x) = &bweight{
            calc_weighted_stats(&evec,x)
        }else{
            calc_stats(&evec)
        };
        let stdev = sstat.var.sqrt()+0.0000001;
        unsafe{
            element_add(&mut evec,-1.0*sstat.mean);
            element_multiply(&mut evec, 1.0/stdev/2.0);
        }
        ret.push(evec);
    }
    
    for cc in 0..blen{
        let mut evec:Vec<f32> = vec![];    
        for rr in 0..alen{
            let distt = calc_euclid_dist(&avec[rr], &bvec[cc])/vecsiz;//すぐに 0 とか INF になるのでベクトルのサイズで割る
            evec.push((-1.0*distt).exp());
        }
        let sstat:VectorStats = if let Some(x) = &aweight{
            calc_weighted_stats(&evec,x)
        }else{
            calc_stats(&evec)
        };
        let stdev = sstat.var.sqrt()+0.0000001;
        unsafe{
            element_add(&mut evec,-1.0*sstat.mean);
           element_multiply(&mut evec, 1.0/stdev/2.0);
        }
        for rr in 0..alen{
            ret[rr][cc] += evec[rr];
        }
    }

    return ret;
}

 
pub fn calc_dot_product_matrix(avec:&Vec<&Vec<f32>>,bvec:&Vec<&Vec<f32>>)->Vec<Vec<f32>>{
    let alen = avec.len();
    let blen = bvec.len();
    let mut ret:Vec<Vec<f32>> = vec![];
    for rr in 0..alen{
        let mut evec:Vec<f32> = vec![];
        for cc in 0..blen{
            unsafe{
                let score = matrix_process::dot_product(avec[rr], bvec[cc]);
                evec.push(score);
            }
        }
        ret.push(evec);
    }
    return ret;
}
#[cfg(test)]
mod tests{
    use super::gmat::{calc_vec_stats,GMatStatistics};

    #[test]
    fn stattest(){
        let filenames:Vec<String> = vec![
            "./example_files/test1.gmat".to_owned(),
            "./example_files/test2.gmat".to_owned(),
        ];
        unsafe{
            let res = calc_vec_stats(& filenames);
            let mut sstr:Vec<String> = vec![];
            for rr in res.iter(){
                sstr.push(rr.get_string());
            }
            let mut chk:Vec<GMatStatistics> = vec![];
            chk.push( GMatStatistics{  mean: 10.533,  max: 20.000,  var: 104.382,  min: -5.000,  sum: 158.0,  count: 15 });
            chk.push( GMatStatistics{  mean: 8.400,  max: 20.000,  var: 40.640,  min: -6.000,  sum: 126.0,  count: 15 });
            chk.push( GMatStatistics{  mean: 10.867,  max: 30.000,  var: 46.782,  min: 3.000,  sum: 163.0,  count: 15 });
            chk.push( GMatStatistics{  mean: 12.267,  max: 40.000,  var: 78.329,  min: 2.000,  sum: 184.0,  count: 15 });
            for ii in 0..chk.len(){
                assert!(
                    (chk[ii].max - res[ii].max).abs() < 0.003
                );
                assert!(
                    (chk[ii].min - res[ii].min).abs() < 0.003
                );
                assert!(
                    (chk[ii].mean - res[ii].mean).abs() < 0.003
                );
                assert!(
                    (chk[ii].sum - res[ii].sum).abs() < 0.003
                );
                assert!(
                    (chk[ii].var - res[ii].var).abs() < 0.003
                );
                assert!(
                    (chk[ii].count as f32 - res[ii].count as f32).abs() < 0.003
                );
                
            }
            let mut res:Vec<GMatStatistics> = vec![];
            for ii in 0..sstr.len(){
                res.push(GMatStatistics::load_string(&sstr[ii]));
            }
            for ii in 0..chk.len(){
                assert!(
                    (chk[ii].max - res[ii].max).abs() < 0.003
                );
                assert!(
                    (chk[ii].min - res[ii].min).abs() < 0.003
                );
                assert!(
                    (chk[ii].mean - res[ii].mean).abs() < 0.003
                );
                assert!(
                    (chk[ii].sum - res[ii].sum).abs() < 0.003
                );
                assert!(
                    (chk[ii].var - res[ii].var).abs() < 0.003
                );
                assert!(
                    (chk[ii].count as f32 - res[ii].count as f32).abs() < 0.003
                );
                
            }
        }
    }
}