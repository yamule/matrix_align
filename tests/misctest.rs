
//use matrix_align;


#[cfg(test)]
mod tests{
    //extern crate matrix_align;
    use matrix_align::pssm::{self, calc_vec_stats, PssmStatistics};
    use matrix_align::aligner::{ScoredSeqAligner,ScoredSequence};
    use matrix_align::ioutil::load_pssm_matrix;


    #[test]
    fn aligntest(){
        
        let filenames:Vec<String> = vec![
            "./example_files/esm2_650m_example_output/d6iyia_.res.gz".to_owned(),
            "./example_files/esm2_650m_example_output/d7diha_.res.gz".to_owned()
        ];
        let pssmstats:Vec<PssmStatistics>;
        unsafe{
            pssmstats = calc_vec_stats(& filenames);
        }
        let mut pssm1 = load_pssm_matrix("./example_files/esm2_650m_example_output/d6iyia_.res.gz",true);
        //let pssm2 = ioutil::load_pssm_matrix("./example_files/esm2_650m_example_output/d6iyia_.res.gz",true);
        let mut pssm2 = load_pssm_matrix("./example_files/esm2_650m_example_output/d7diha_.res.gz",true);
        pssm::normalize_seqmatrix(&mut (pssm1.1), &pssmstats);
        pssm::normalize_seqmatrix(&mut (pssm2.1), &pssmstats);

        //let pssm1 = ioutil::load_pssm_matrix("./example_files/test1.pssm",false);
        //let pssm2 = ioutil::load_pssm_matrix("./example_files/test2.pssm",false);
        let mut saligner:ScoredSeqAligner = ScoredSeqAligner::new(pssm1.1[0].len(),200,100);
        let seq1 = ScoredSequence::new(
            pssm1.0,pssm1.1,&mut saligner,true
        );

        let seq2 = ScoredSequence::new(
            pssm2.0,pssm2.1,&mut saligner,true
        );

        let res = saligner.perform_dp(
            &seq1,&seq2,-10.0,-0.5
        );

        /*
        println!("{:?}",res);

        for cc in 0..3{
            for ii in 0..11{
                for jj in 0..6{
                    print!("{:>.1} ",aligner.dp_matrix[ii][jj][cc]);
                }
                println!("");
            }
            println!("");
        }
        */
        let plen = res.0.len();
        let apos = seq1.alibuff_idx[0];
        let bpos = seq2.alibuff_idx[0];
        let mut aseq = "".to_owned();
        let mut bseq = "".to_owned();
        for rr in res.0.iter(){
            if rr.0 > -1{
                aseq += &saligner.alignment_buffer[apos][rr.0 as usize].to_string();
            }else{
                aseq += "-";
            }
            if rr.1 > -1{
                bseq += &saligner.alignment_buffer[bpos][rr.1 as usize].to_string();
            }else{
                bseq += "-";
            }
        }
        println!("{}",aseq);
        println!("{}",bseq);
    }
}