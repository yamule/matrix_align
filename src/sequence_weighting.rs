


//Henikoff & Henikoff 1994 の全長にわたるウエイトを計算する
//余り頻繁に呼ばれることは想定していない
pub fn calc_henikoff_henikoff_weight(seqs:& Vec<Vec<u8>>)->Vec<f32>{
    let mut ret:Vec<f32> = vec![0.0;seqs.len()];
    for cc in 0..seqs[0].len(){
        let mut counter:Vec<f32> = vec![0.0;256];
        for rr in 0..seqs.len(){
            counter[seqs[rr][cc] as usize] += 1.0;
        }
        for rr in 0..seqs.len(){
            assert!(counter[seqs[rr][cc] as usize] != 0.0);
            ret[rr] += 1.0/counter[seqs[rr][cc] as usize];
        }
    }
    for rr in 0..seqs.len(){
        ret[rr] /= seqs[0].len() as f32;
    }
    return ret;
}