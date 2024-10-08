#![allow(unused_unsafe)]

#[cfg(test)]
mod tests{
    use std::collections::HashMap;

    //extern crate matrix_align;
    use matrix_align::gmat::{self, calc_vec_stats, GMatStatistics};
    use matrix_align::aligner::{AlignmentType, ProfileAligner, SequenceProfile};
    use matrix_align::ioutil::{load_multi_gmat,save_lines};
    use matrix_align::guide_tree_based_alignment::{self, DistanceBase};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    #[test]
    fn aligntest(){
        let infile = "example_files/esm2_650m_example_output/a1_mat.dat".to_owned();
        let outfile = "nogit/testout.test.dat".to_owned();
        let num_iter:usize = 3;
        
        let mut name_to_res:HashMap<String,String> = HashMap::new();
        
        let gmatstats:Vec<GMatStatistics>;
        let mut profile_seq:Option<SequenceProfile> = None;
        unsafe{
            gmatstats = calc_vec_stats(& vec![infile.clone()],4,None);
        }
        
        let mut name_order:Vec<String> = vec![];
        let gmat1_ = load_multi_gmat(&infile,infile.ends_with(".gz"));
        for ii in 0..num_iter{
            eprintln!("{}",ii);

            let veclen = gmat1_[0].2[0].len();
            let mut saligner:ProfileAligner = ProfileAligner::new(veclen,200,Some(-10.0),AlignmentType::Global,matrix_align::aligner::ScoreType::DistanceZscore,None
            ,true);
            let mut seqvec:Vec<SequenceProfile> = vec![];
            
            if let Some(p) = profile_seq{
                seqvec.push(p);
            }

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
                let seq2 = SequenceProfile::new(
                    vec![(tt.0,tt.1)],tlen,tt.2[0].len(),None,Some(tt.2),None
                );    
                seqvec.push(seq2);
            }

            let mut ans = saligner.make_msa(seqvec,false);
            assert!(ans.len() == 1);
            let (mut alires,_alisc) = ans.pop().unwrap();
            if num_iter == ii+1{
                for aii in 0..alires.member_sequences.len(){
                    let hh = &alires.headers[aii];
                    let res = alires.get_aligned_seq(aii);
                    if name_to_res.contains_key(hh){
                        name_to_res.insert(
                            hh.clone(),res.into_iter().map(|m|m.to_string()).collect::<Vec<String>>().concat()
                        );
                    }
                }
                break;
            }
            println!("{}",alires.member_sequences.len());
            alires.member_sequences.clear();
            alires.headers.clear();
            alires.alignment_mapping.clear();
            alires.alignment_mapping_ids.clear();
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

    #[test]
    fn hierarchical_aligntest(){
        let infile = "example_files/esm2_650m_example_output/a1_mat.dat".to_owned();
        let outfile = "nogit/testout_hierarchical.test.dat".to_owned();
        
        let mut name_to_res:HashMap<String,String> = HashMap::new();
        
        let mut rngg:StdRng = StdRng::seed_from_u64(123);
        let gmatstats:Vec<GMatStatistics>;

        unsafe{
            gmatstats = calc_vec_stats(& vec![infile.clone()],2,None);
        }
        
        let mut name_order:Vec<String> = vec![];
        let gmat1__ = load_multi_gmat(&infile,infile.ends_with(".gz"));

        let mut gmat1_:Vec<(String,Vec<char>,Vec<Vec<f32>>,Option<Vec<(f32,f32,f32,f32)>>)> = vec![];
        for gg in gmat1__.into_iter(){

            for ii in 0..200{//200 で 20 秒くらい
                let mut gp = gg.clone();
                for jj in gp.2.iter_mut(){
                    for kk in jj.iter_mut().enumerate(){
                        let r = rngg.gen_range((gmatstats[kk.0].min/4.0)..(gmatstats[kk.0].max/4.0));
                        *kk.1 += r;
                    }
                }
                gp.0 = gp.0+"_"+ii.to_string().as_str();
                gmat1_.push(gp);
            }
        }

        let veclen = gmat1_[0].2[0].len();
        let mut saligner:ProfileAligner = ProfileAligner::new(veclen,200,Some(-10.0),AlignmentType::Global,matrix_align::aligner::ScoreType::DistanceZscore,None,true);
        let mut seqvec:Vec<SequenceProfile> = vec![];
        
        let gmat1 = gmat1_.clone();
        for mut tt in gmat1.into_iter(){
            let n = tt.0.clone();
            if name_to_res.contains_key(&n){
                panic!("Duplicated name {}.",n);
            }
            name_to_res.insert(n.clone(),"".to_owned());
            name_order.push(n);
        
            gmat::normalize_seqmatrix(&mut (tt.2), &gmatstats);
            let tlen = tt.1.len();
            let seq2 = SequenceProfile::new(
                vec![(tt.0,tt.1)],tlen,tt.2[0].len(),None,Some(tt.2),None
            );    
            seqvec.push(seq2);
        }

        let mut ans = guide_tree_based_alignment::hierarchical_alignment(
            seqvec, &DistanceBase::AveragedValue, &mut saligner,500,&mut rngg,8,guide_tree_based_alignment::TreeType::TreeNj);
        assert!(ans.len() == 1);
        let (alires,_alisc) = ans.pop().unwrap();
        for aii in 0..alires.member_sequences.len(){
            let hh = &alires.headers[aii];
            let res = alires.get_aligned_seq(aii);
            if name_to_res.contains_key(hh){
                name_to_res.insert(
                    hh.clone(),res.into_iter().map(|m|m.to_string()).collect::<Vec<String>>().concat()
                );
            }
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