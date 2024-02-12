use self::matrix_process::{calc_euclid_dist, calc_stats, element_add, element_multiply, VectorStats};

use super::*;

#[derive(Debug, PartialEq)]
pub struct PssmStatistics{
    pub max:f32,
    pub min:f32,
    pub sum:f32,
    pub mean:f32,
    pub var:f32,
    pub count:usize
}

pub unsafe fn calc_vec_stats(filenames:&Vec<String>)->Vec<PssmStatistics>{
    let mut ssum:Vec<f32> = vec![];
    let mut smax:Vec<f32> = vec![];
    let mut smin:Vec<f32> = vec![];
    let mut counter = 0_usize;
    
    let mut allval:Vec<Vec<f32>> = vec![];//まあ多分メモリ上に乘るだろう。。。
    for fname in filenames.into_iter(){
        let mut pssm1 = ioutil::load_pssm_matrix(fname,fname.ends_with(".gz"));
        if ssum.len() == 0{
            ssum = vec![0.0;pssm1.1[0].len()];
            smax = pssm1.1[0].clone();
            smin = pssm1.1[0].clone();
        }
        for pp in pssm1.1.iter(){
            matrix_process::vector_add(&mut ssum, &pp);
            matrix_process::vector_max(&mut smax, &pp);
            matrix_process::vector_min(&mut smin, &pp);
            counter += 1;
        }
        allval.append(&mut pssm1.1);
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
    let mut ret:Vec<PssmStatistics> = vec![];
    for ii in 0..ssum.len(){
        ret.push(
            PssmStatistics{
                max:smax[ii],min:smin[ii],sum:ssum[ii],mean:smean[ii],var:svar[ii],count:counter
            }
        );
        //println!("mean:\t{}\tmax:\t{}\tmin:\t{}\tvar:{}",smean[ii],smax[ii],smin[ii],svar[ii]/(counter as f32));
    }
    return ret;
}


pub fn normalize(vec:&mut Vec<f32>,pssmstats:&Vec<PssmStatistics>){
    assert_eq!(vec.len(),pssmstats.len());
    for vv in vec.iter_mut().zip(pssmstats.iter()){
        assert!(vv.1.var > 0.0);
        *vv.0 = (*vv.0 -vv.1.mean)/vv.1.var;
    }
}

pub fn normalize_seqmatrix(vec:&mut Vec<Vec<f32>>, pssmstats:&Vec<PssmStatistics>){
    let vlen = vec.len();
    for ii in 0..vlen{
        normalize(&mut vec[ii], pssmstats);
    }
}

// Match State で加算される、スコアの行列を、(avec.len(),bvec.len()) の長さの行列として返す。
// 要素は一方の行列のある行と他方の行列の全行のユークリッド距離を計算し、統計値を取る。
// その統計値をもとに ZSCORE を計算し
// avec, bvec 両方で同様のことを行って平均をとって要素とする。
//Pantolini, Lorenzo, et al. "Embedding-based alignment: combining protein language models with dynamic programming alignment to detect structural similarities in the twilight-zone." Bioinformatics 40.1 (2024): btad786.
pub fn calc_dist_zscore_matrix(avec:& Vec<&Vec<f32>>,bvec:& Vec<&Vec<f32>>)-> Vec<Vec<f32>>{
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
        let sstat:VectorStats = calc_stats(&evec);
        let stdev = sstat.var.sqrt()+0.0000001;
        element_add(&mut evec,-1.0*sstat.mean);
        element_multiply(&mut evec, 1.0/stdev/2.0);// 逆方向の値との平均を取るので 2.0 で割る
        ret.push(evec);
    }
    
    for cc in 0..blen{
        let mut evec:Vec<f32> = vec![];    
        for rr in 0..alen{
            let distt = calc_euclid_dist(&avec[rr], &bvec[cc])/vecsiz;
            evec.push((-1.0*distt).exp());
        }
        let sstat:VectorStats = calc_stats(&evec);
        let stdev = sstat.var.sqrt()+0.0000001;
        element_add(&mut evec,-1.0*sstat.mean);
        element_multiply(&mut evec, 1.0/stdev);
        for rr in 0..alen{
            ret[rr][cc] += evec[rr]/2.0;
        }
    }
    return ret;
}


#[cfg(test)]
mod tests{
    use super::pssm::{calc_vec_stats,PssmStatistics};

    #[test]
    fn stattest(){
        let filenames:Vec<String> = vec![
            "./example_files/test1.pssm".to_owned(),
            "./example_files/test2.pssm".to_owned(),
        ];
        unsafe{
            let res = calc_vec_stats(& filenames);
            let mut chk:Vec<PssmStatistics> = vec![];
            chk.push( PssmStatistics{  mean: 10.533,  max: 20.000,  var: 104.382,  min: -5.000,  sum: 158.0,  count: 15 });
            chk.push( PssmStatistics{  mean: 8.400,  max: 20.000,  var: 40.640,  min: -6.000,  sum: 126.0,  count: 15 });
            chk.push( PssmStatistics{  mean: 10.867,  max: 30.000,  var: 46.782,  min: 3.000,  sum: 163.0,  count: 15 });
            chk.push( PssmStatistics{  mean: 12.267,  max: 40.000,  var: 78.329,  min: 2.000,  sum: 184.0,  count: 15 });
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