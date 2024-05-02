use crate::aligner::DPResult;
use super::aligner::SequenceProfile;
use std::collections::HashMap;
use super::aligner::ProfileAligner;
use rayon::prelude::*;
use super::misc::*;

pub fn create_a3m_pairwise_alignment(aligners:&mut Vec<ProfileAligner>,firstseq:SequenceProfile,mut allseqs_:Vec<SequenceProfile>)->HashMap<String,(String,DPResult)>{
    
    let mut ret:HashMap<String,(String,DPResult)> = HashMap::new();

    println!("Create a3m alignment with {} sequences...",allseqs_.len());
    
    while allseqs_.len() > 0{
        let mut minibatch:Vec<
        (SequenceProfile,SequenceProfile,ProfileAligner)
        > = vec![];
        while allseqs_.len() > 0{
            let bseq = allseqs_.pop().unwrap();
            let ali = aligners.pop().unwrap();
            let fst = firstseq.clone();
            minibatch.push((fst,bseq,ali));
            if aligners.len() == 0{
                break;
            }
        }


        let results:Vec<(ProfileAligner,String,String,DPResult)> = minibatch.into_par_iter().map(|v|{
            let fst = v.0;
            let bseq = v.1;

            let fname = fst.headers[0].clone();
            let bname = bseq.headers[0].clone();
            let mut ali = v.2;
            let mut dpres = ali.perform_dp(&fst,&bseq);
            
            let newgroup = ali.make_alignment(fst,bseq,dpres.alignment.clone(),false,None);
            
            //スコア以外は捨てる
            dpres.alignment.clear();
            
            let mut name_to_res:HashMap<String, String> = HashMap::new();
            insert_alinged_string(&newgroup, &mut name_to_res, false);
            let first_aa:Vec<char> = name_to_res.get(&fname).unwrap().chars().into_iter().collect();
            let pres:Vec<char> = name_to_res.get(&bname).unwrap().chars().into_iter().collect();
            assert_eq!(pres.len(),first_aa.len());
            let mut pstr:Vec<String> = vec![];
            for (aa,bb) in first_aa.iter().zip(pres.iter()){
                if aa == &'-'{
                    pstr.push(
                        bb.to_ascii_lowercase().to_string()
                    );
                }else{
                    pstr.push(
                        bb.to_string()
                    );
                }
            }
            
            (ali,bname,pstr.join(""),dpres)
        }).collect();
        for rr in results.into_iter(){
            aligners.push(rr.0);
            ret.insert(
                rr.1,(rr.2,rr.3)
            );
        }
    }
    return ret;
}