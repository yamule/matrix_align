


//Henikoff & Henikoff 1994 の全長にわたるウエイトを計算する
//余り頻繁に呼ばれることは想定していない
pub fn calc_henikoff_henikoff_weight(seqs:& Vec<Vec<char>>)->Vec<f32>{
    let alilen = seqs[0].len();
    for aa in seqs.iter(){
        if alilen != aa.len(){
            panic!("Sequences must have been aligned.\n{:?}\n{:?}",seqs[0],aa);
        }
    }
    
    let mut ret:Vec<f32> = vec![0.0;alilen];
    for cc in 0..alilen{
        let mut counter:Vec<f32> = vec![0.0;256];
        for rr in 0..seqs.len(){
            let pp = seqs[rr][cc] as i8;
            assert!( pp > -1,"invalid char found {}",seqs[rr][cc]);
            counter[pp as usize] += 1.0;
        }
        for rr in 0..seqs.len(){
            assert!(counter[seqs[rr][cc] as usize] != 0.0);
            ret[rr] += 1.0/counter[seqs[rr][cc] as usize];
        }
    }
    let wsum:f32 = ret.iter().sum();
    for rr in 0..seqs.len(){
        ret[rr] /= wsum;
    }
    return ret;
}