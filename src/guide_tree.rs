


// 対象な距離行列を半分の三角形のみで表現するための関数
// 自分自身もあり
pub fn calc_pos(aa:usize,bb:usize) -> usize{
    let mut v1:f64 = aa as f64;
    let mut v2:f64 = bb as f64;
    if v2 < v1{// v2 の方が常に大きい
        let c = v1;
        v1 = v2;
        v2 = c;
    }
    return (v2*v2/2.0+v2/2.0+v1 +0.000001) as usize;
}
