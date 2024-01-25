


//Henikoff & Henikoff 1994 の全長にわたるウエイトを計算する
//余り頻繁に呼ばれることは想定していない
pub fn calc_henikoff_henikoff_weight(seqs:& Vec<Vec<char>>,ids:& Vec<usize>,alilen:usize)->Vec<f32>{
    let mut ret:Vec<f32> = vec![0.0;alilen];
    for cc in 0..alilen{
        let mut counter:Vec<f32> = vec![0.0;256];
        for rr in ids.iter(){
            let pp = seqs[*rr][cc] as i8;
            assert!( pp > -1,"invalid char found {}",seqs[*rr][cc]);
            counter[pp as usize] += 1.0;
        }
        for rr in ids.iter(){
            assert!(counter[seqs[*rr][cc] as usize] != 0.0);
            ret[*rr] += 1.0/counter[seqs[*rr][cc] as usize];
        }
    }
    let wsum:f32 = ret.iter().sum();
    for rr in 0..seqs.len(){
        ret[rr] /= wsum;
    }
    return ret;
}