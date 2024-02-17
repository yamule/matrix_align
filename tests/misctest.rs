
//use matrix_align;


#[cfg(test)]
mod tests{
    use std::collections::HashMap;

    //extern crate matrix_align;
    use matrix_align::gmat::{self, calc_vec_stats, GMatStatistics};
    use matrix_align::aligner::{ScoredSeqAligner,ScoredSequence};
    use matrix_align::ioutil::{load_gmat, load_multi_gmat};
    use matrix_align::misc::*;


    #[test]
    fn aligntest(){
        let filename = "example_files/esm2_650m_example_output/a1_mat.dat".to_owned();
        let gmatstats:Vec<GMatStatistics>;
        unsafe{
            gmatstats = calc_vec_stats(& vec![filename.clone()]);
        }
        let gmat1 = load_multi_gmat(&filename,filename.ends_with(".gz"));
        let veclen = gmat1[0].2[0].len();
        let mut saligner:ScoredSeqAligner = ScoredSeqAligner::new(veclen,200,100);
        let mut seqvec:Vec<ScoredSequence> = vec![];
        for mut tt in gmat1.into_iter(){
            gmat::normalize_seqmatrix(&mut (tt.2), &gmatstats);
            let seq2 = ScoredSequence::new(
                tt.1,tt.2,None,&mut saligner,true
            );
            seqvec.push(seq2);
        }

        let ans = saligner.make_msa(seqvec,-10.0,-0.5,false);
        let alires = &ans.0[0];
        let mut center_char:Vec<char> = vec![];
        let mut center_gmat:Vec<Vec<f32>> = vec![];
        let mut center_gaps:Vec<(f32,f32,f32,f32)> = vec![];
        for (eii,aa) in alires.alibuff_idx.iter().enumerate(){
            for ii in 0..alires.alignment_length{
                print!("{}",saligner.ali_get(*aa,ii));
            }
            if alires.primary_ids[eii] == 0{
                center_char = vec![];
                center_gmat = vec![];
                for _ in 0..alires.alignment_length{
                    //center_char.push(saligner.ali_get(*aa,ii));
                    center_char.push('X');
                }
                for ii in 0..=alires.alignment_length{
                    let ppp = saligner.gmat_colget(alires.gmatbuff_id  as usize,ii).clone();
                    center_gaps.push((ppp.match_weight,ppp.del_weight,ppp.connected_weight,ppp.gapped_weight));
                    if ii < alires.alignment_length{
                        center_gmat.push(ppp.match_vec);
                    }
                }
            }
            println!("");
        }
        let mut saligner:ScoredSeqAligner = ScoredSeqAligner::new(veclen,200,100);
        let dummy_center = ScoredSequence::new(
            center_char,center_gmat, Some(center_gaps),&mut saligner,true
        );
        saligner.weights[dummy_center.primary_ids[0] as usize] = 1.0;
        let mut seqvec:Vec<ScoredSequence> = vec![dummy_center];
        
        let mut primary_id_to_name:HashMap<i32,String> = HashMap::new();
        let mut id_order:Vec<i32> = vec![];
        let gmat1 = load_multi_gmat(&filename,filename.ends_with(".gz"));
        for mut tt in gmat1.into_iter(){
            gmat::normalize_seqmatrix(&mut (tt.2), &gmatstats);
            let seq2 = ScoredSequence::new(
                tt.1,tt.2,None,&mut saligner,true
            );
            assert!(seq2.primary_ids.len() == 1);
            primary_id_to_name.insert(seq2.primary_ids[0],tt.0);
            id_order.push(seq2.primary_ids[0]);
            seqvec.push(seq2);
        }
        
        let ans = saligner.make_msa(seqvec,-10.0,-0.5,false);
        let alires = &ans.0[0];
        let mut primary_id_to_res:HashMap<i32,String> = HashMap::new();
        for (eii,aa) in alires.alibuff_idx.iter().enumerate(){
            let mut rres:Vec<char> = vec![];
            for ii in 0..alires.alignment_length{
                rres.push(saligner.ali_get(*aa,ii))
            }
            primary_id_to_res.insert(
                alires.primary_ids[eii],rres.into_iter().map(|m|m.to_string()).collect::<Vec<String>>().concat()
            );
        }
        for ii in id_order.iter(){
            println!(">{}",primary_id_to_name.get(ii).unwrap());
            println!("{}",primary_id_to_res.get(ii).unwrap());
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