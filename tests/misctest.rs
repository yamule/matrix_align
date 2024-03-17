
//use matrix_align;


#[cfg(test)]
mod tests{
    use std::collections::HashMap;

    //extern crate matrix_align;
    use matrix_align::gmat::{self, calc_vec_stats, GMatStatistics};
    use matrix_align::aligner::{AlignmentType, ScoredSeqAligner, ScoredSequence};
    use matrix_align::ioutil::{load_gmat, load_multi_gmat,save_lines};
    use matrix_align::misc::*;


    #[test]
    fn aligntest(){
        let infile = "example_files/esm2_650m_example_output/a1_mat.dat".to_owned();
        let outfile = "nogit/testout.test.dat".to_owned();
        let num_iter:usize = 3;
        
        let mut name_to_res:HashMap<String,String> = HashMap::new();
        
        let gmatstats:Vec<GMatStatistics>;
        let mut profile_seq:Option<ScoredSequence> = None;
        unsafe{
            gmatstats = calc_vec_stats(& vec![infile.clone()]);
        }
        
        let mut name_order:Vec<String> = vec![];
        let gmat1_ = load_multi_gmat(&infile,infile.ends_with(".gz"));
        for ii in 0..num_iter{
            eprintln!("{}",ii);

            let veclen = gmat1_[0].2[0].len();
            let mut saligner:ScoredSeqAligner = ScoredSeqAligner::new(veclen,200,-10.0,-0.5,AlignmentType::Global);
            let mut seqvec:Vec<ScoredSequence> = vec![];
            
            if let Some(p) = profile_seq{
                seqvec.push(p);
            }
            profile_seq = None;

            let gmat1 = gmat1_.clone();
            for mut tt in gmat1.into_iter(){
                if ii == 0{
                    let n = tt.0.clone();
                    if name_to_res.contains_key(&n){
                        panic!("Duplicated name {}.",n);
                    }
                    name_to_res.insert(n.clone(),"".to_owned());
                    name_order.push(n);
                }
                gmat::normalize_seqmatrix(&mut (tt.2), &gmatstats);
                let tlen = tt.1.len();
                let seq2 = ScoredSequence::new(
                    vec![(tt.0,tt.1)],tlen,tt.2[0].len(),None,Some(tt.2),None
                );    
                seqvec.push(seq2);
            }

            let mut ans = saligner.make_msa(seqvec,false);
            assert!(ans.len() == 1);
            let (alires,_alisc) = ans.pop().unwrap();
            if num_iter == ii+1{
                for (ali,hh) in alires.alignments.into_iter().zip(alires.headers.into_iter()){
                    if name_to_res.contains_key(&hh){
                        name_to_res.insert(
                            hh,ali.into_iter().map(|m|m.to_string()).collect::<Vec<String>>().concat()
                        );
                    }
                }
                break;
            }
            profile_seq = Some(alires);
        }
        let mut results:Vec<String> = vec![];
        for ii in name_order.iter(){
            //println!(">{}",primary_id_to_name.get(ii).unwrap());
            //println!("{}",primary_id_to_res.get(ii).unwrap());
            results.push(
                format!(">{}\n{}",ii,name_to_res.get(ii).unwrap())
            );
            println!(">{}\n{}",ii,name_to_res.get(ii).unwrap())
        }
        save_lines(&outfile, results,outfile.ends_with(".gz"));
    }

}