
//use matrix_align;


#[cfg(test)]
mod tests{
    //extern crate matrix_align;
    use matrix_align::pssm::{self, calc_vec_stats, PssmStatistics};
    use matrix_align::aligner::{ScoredSeqAligner,ScoredSequence};
    use matrix_align::ioutil::load_pssm_matrix;
    use matrix_align::misc::*;


    #[test]
    fn aligntest(){
        let filenames:Vec<String> = vec![
            "./example_files/esm2_650m_example_output/d6iyia_.res.gz".to_owned(),
            "./example_files/esm2_650m_example_output/d7diha_.res.gz".to_owned(), //fam
            "./example_files/esm2_650m_example_output/d4f6ia_.res.gz".to_owned(), //supfam
            //"./example_files/esm2_650m_example_output/d6lumb2.res.gz".to_owned() //fold        
        ];
        let pssmstats:Vec<PssmStatistics>;
        unsafe{
            pssmstats = calc_vec_stats(& filenames);
        }
        let pssm1_ = load_pssm_matrix(&filenames[0],true);
        let mut saligner:ScoredSeqAligner = ScoredSeqAligner::new(pssm1_.1[0].len(),200,100);
        let mut seqvec:Vec<ScoredSequence> = vec![];
        for ii in 0..filenames.len(){
            let mut pssm2 = load_pssm_matrix(&filenames[ii],true);
            pssm::normalize_seqmatrix(&mut (pssm2.1), &pssmstats);
            let seq2 = ScoredSequence::new(
                pssm2.0,pssm2.1,None,&mut saligner,true
            );
            seqvec.push(seq2);
        }
        let ans = saligner.make_msa(seqvec,-10.0,-0.5,false);
        let alires = &ans.0[0];
        let mut center_char:Vec<char> = vec![];
        let mut center_pssm:Vec<Vec<f32>> = vec![];
        let mut center_gaps:Vec<(f32,f32,f32,f32)> = vec![];
        for (eii,aa) in alires.alibuff_idx.iter().enumerate(){
            for ii in 0..alires.alignment_length{
                print!("{}",saligner.ali_get(*aa,ii));
            }
            if alires.primary_ids[eii] == 0{
                center_char = vec![];
                center_pssm = vec![];
                for _ in 0..alires.alignment_length{
                    //center_char.push(saligner.ali_get(*aa,ii));
                    center_char.push('X');
                }
                for ii in 0..=alires.alignment_length{
                    let ppp = saligner.pssm_colget(alires.pssmbuff_id  as usize,ii).clone();
                    center_gaps.push((ppp.match_weight,ppp.del_weight,ppp.connected_weight,ppp.gapped_weight));
                    if ii < alires.alignment_length{
                        center_pssm.push(ppp.match_vec);
                    }
                }
            }
            println!("");
        }
        let mut saligner:ScoredSeqAligner = ScoredSeqAligner::new(pssm1_.1[0].len(),200,100);
        let dummy_center = ScoredSequence::new(
            center_char,center_pssm, Some(center_gaps),&mut saligner,true
        );
        saligner.weights[dummy_center.primary_ids[0] as usize] = 1.0;
        let didx = dummy_center.alibuff_idx[0];
        let mut seqvec:Vec<ScoredSequence> = vec![dummy_center];
        for ii in 0..filenames.len(){
            let mut pssm2 = load_pssm_matrix(&filenames[ii],true);
            pssm::normalize_seqmatrix(&mut (pssm2.1), &pssmstats);
            let seq2 = ScoredSequence::new(
                pssm2.0,pssm2.1, None,&mut saligner,true
            );
            seqvec.push(seq2);
        }

        let ans = saligner.make_msa(seqvec,-10.0,-0.5,false);
        let alires = &ans.0[0];
        
        for (eii,aa) in alires.alibuff_idx.iter().enumerate(){
            for ii in 0..alires.alignment_length{
                print!("{}",saligner.ali_get(*aa,ii));
            }
            println!("");
        }
        

        /*
            let res = saligner.perform_dp(
                &seq1,&seq2,-10.0,-0.5
            );


            println!("file:{}\nScore: {}",filenames[ii],res.1);

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
            let apos = seq1.alibuff_idx[0];
            let bpos = seq2.alibuff_idx[0];
            let mut aseq:Vec<char> = vec![];
            let mut bseq:Vec<char> = vec![];
            for rr in res.0.iter(){
                if rr.0 > -1{
                    aseq.push(saligner.alignment_buffer[apos][rr.0 as usize]);
                }else{
                    aseq.push('-');
                }
                if rr.1 > -1{
                    bseq.push(saligner.alignment_buffer[bpos][rr.1 as usize]);
                }else{
                    bseq.push('-');
                }
            }
            
            let alires = saligner.make_alignment(seq1,seq2,res.0 ,false);
            for aa in alires.alibuff_idx.iter(){
                for ii in 0..alires.alignment_length{
                    print!("{}",saligner.ali_get(*aa,ii));
                }
                println!("");
            }
            
            let afas:String =aseq.iter().collect::<String>();
            let araw:String = aseq.iter().filter(|c| **c != '-').collect::<String>();
            println!(">seq1\n{}",afas);
            println!(">seq2\n{}",bseq.iter().collect::<String>());
            
            let ba3m:String = subject_ali_to_a3m(&aseq,&bseq).concat();
            /*
            println!("{}",araw);
            println!("{}",ba3m);
            */
            a3mres.push(ba3m.clone());
            assert_eq!(
            query_a3m_to_ali(&araw.chars().into_iter().collect::<Vec<char>>(),&ba3m.chars().into_iter().collect::<Vec<char>>()).concat()
            ,afas);
        }
        */
        /*
        for ii in 0..filenames.len(){
            println!(">{}",filenames[ii]);
            println!("{}",a3mres[ii]);
        }
        */
    }
}