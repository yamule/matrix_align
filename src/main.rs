use std::collections::*;
use matrix_align::gmat::{self, calc_vec_stats, calc_vec_stats_, GMatStatistics};
use matrix_align::aligner::{ScoredSeqAligner,ScoredSequence};
use matrix_align::ioutil::{load_multi_gmat, save_lines};
use matrix_align::matrix_process;
use matrix_align::misc::FloatWrap;

fn argparse(args:Vec<String>)->HashMap<String,String>{
    let mut ret:HashMap<String,String> = HashMap::new();
    let mut non_key_count:usize = 0;
    let mut ii:usize = 0;
    while(ii < args.len()){
        if args[ii].starts_with("--"){
            assert!(ii < args.len()-1,"{} does not have a value.",args[ii]);
            assert!(!ret.contains_key(&args[ii]),"{} has already been assigned.",args[ii]);
            if ii == args.len()-1 || args[ii+1].starts_with("--"){
                ret.insert(args[ii].clone(),"true".to_owned());
            }else{
                ret.insert(args[ii].clone(),args[ii+1].clone());
                ii += 1;
            }
        }else{
            ret.insert(format!("nokey_{}",non_key_count),args[ii].clone());
            non_key_count += 1;
        }
        ii += 1;
    }
    return ret;
}

//ToDo argparser 作る

fn main(){
    let argss: HashMap<String,String> = argparse(std::env::args().collect::<Vec<String>>());
    let required_arg:Vec<&str> = vec!["--in","--out"];
    for rr in required_arg.iter(){
        if !argss.contains_key(*rr){
            panic!("{} is required.",rr);
        }
    }
    let infile = argss.get("--in").unwrap().clone();
    let outfile = argss.get("--out").unwrap().clone();
    let num_iter:usize = argss.get("--num_iter").unwrap_or(&"2".to_owned()).parse::<usize>().unwrap_or_else(|e|panic!("in --num_iter {:?}",e));
    let gap_open_penalty:f32 = argss.get("--gap_open_penalty").unwrap_or(&"-10.0".to_owned()).parse::<f32>().unwrap_or_else(|e|panic!("in --gap_open_penalty {:?}",e));
    let gap_extension_penalty:f32 = argss.get("--gap_extension_penalty").unwrap_or(&"-0.5".to_owned()).parse::<f32>().unwrap_or_else(|e|panic!("in --gap_extension_penalty {:?}",e));
    let normalize:bool = argss.contains_key("--normalize");
    
    let mut name_to_res:HashMap<String,String> = HashMap::new();
    
    let gmatstats:Vec<GMatStatistics>;
    let mut profile_seq:Option<ScoredSequence> = None;
    unsafe{
        gmatstats = calc_vec_stats(& vec![infile.clone()]);// 統計値のために一回ファイルを読んでいるが後で変更する
    }

    let mut name_order:Vec<String> = vec![];
    let gmat1_ = load_multi_gmat(&infile,infile.ends_with(".gz"));

    let veclen = gmat1_[0].2[0].len();
    let mut saligner:ScoredSeqAligner = ScoredSeqAligner::new(veclen,300,100);
    
    let mut allseqs_:Vec<ScoredSequence> = vec![];
    for mut gg in gmat1_.into_iter(){
        let n = gg.0.clone();
        if name_to_res.contains_key(&n){
            panic!("Duplicated name {}.",n);
        }
        name_to_res.insert(n.clone(),"".to_owned());
        name_order.push(n);
        if normalize{
            gmat::normalize_seqmatrix(&mut (gg.2), &gmatstats);
        }
        //駄目っぽい。保留。
        //let ss_biased = gmat::ssbias(&mut tt.2,false);
        //tt.2 = ss_biased;
        
        // ギャップまだ入って無い
        let seq = ScoredSequence::new(
            vec![(gg.0,gg.1)],gg.2[0].len(),None,Some(gg.2),None
        );
        allseqs_.push(seq);
    }

    let similarity_sort = false;
    
    if similarity_sort{
        //先頭配列に近い順でアラインメントする
        let seq1 = allseqs_.swap_remove(0);
        let mut scoresort:Vec<(f32,ScoredSequence)> = vec![];
        for seq2 in allseqs_.into_iter(){
            let dpres = saligner.perform_dp(&seq1,&seq2,gap_open_penalty,gap_extension_penalty);
            //let nscore = dpres.1/(seq1.get_alignment_length().min(seq2.get_alignment_length())) as f32;
            //needleman wunsh の合計だと normalize する必要はない気がするが・・・
            //うーん
            //scoresort.push((nscore,seq2));
            scoresort.push((dpres.1,seq2));
        }
        scoresort.sort_by(|a,b|a.0.partial_cmp(&b.0).unwrap());
        scoresort.reverse();
        allseqs_ = vec![seq1];
        for ss in scoresort.into_iter(){
            eprintln!("{} {:?}",ss.0,ss.1.headers);
            allseqs_.push(ss.1);
        }
    }

    for ii in 0..num_iter{
        eprintln!("iter: {}",ii);

        let mut seqvec:Vec<ScoredSequence> = vec![];
        
        if let Some(p) = profile_seq{
            seqvec.push(p.create_merged_msa());
        }
        profile_seq = None;
        for ss in allseqs_.iter(){
            seqvec.push(ss.clone());
        }
        
        let softtree = true;
        
        let mut ans = if softtree{
            let numseq = seqvec.len();
            let mut meanval:Vec<Vec<f32>> = vec![];
            for jj in 0..numseq{
                unsafe{
                    let mut vv:Vec<Vec<f32>> = vec![];
                    for gg in seqvec[jj].gmat.iter(){
                        vv.push(gg.match_vec.clone());
                    }
                    let vstat = calc_vec_stats_(vv);
                    let mut mm:Vec<f32> = vec![];
                    for ss in vstat.into_iter(){
                        mm.push(ss.mean);
                    }
                    meanval.push(mm);
                }
            }
            
            let mut ssorter:Vec<(f32,usize,usize)> = vec![];
            for rr in 0..numseq{
                for cc in (rr+1)..numseq{
                    let distt = matrix_process::calc_euclid_dist(&meanval[rr],& meanval[cc]);
                    ssorter.push((distt,rr,cc));
                }
            }
            ssorter.sort_by(|a,b| a.partial_cmp(&b).unwrap_or_else(||panic!("???")));
            let mut edges:Vec<(usize,usize)> = vec![];
            for ii in 0..ssorter.len()/2{
                edges.push((ssorter[ii].1,ssorter[ii].2));
            }
            println!("{:?}",ssorter);
            println!("{:?}",edges);
            let mut ress = saligner.make_msa_with_edge(seqvec,gap_open_penalty, gap_extension_penalty,false,edges,true);
            assert!(ress.len() == 1);
            ress.pop().unwrap()
        }else{
            saligner.make_msa(seqvec,gap_open_penalty, gap_extension_penalty,false)
        };

        assert!(ans.0.len() == 1);
        eprintln!("score:{}",ans.1);
        let alires = ans.0.pop().unwrap();

        //let pstr = alires.gmat_str();
        //for pp in pstr{ //profile 表示
        //    println!("{}",pp);
        //}

        if num_iter == ii+1{
            for (ali,hh) in alires.alignments.into_iter().zip(alires.headers.into_iter()){
                if name_to_res.contains_key(&hh){
                    assert!(name_to_res.get(&hh).unwrap().len() == 0,"{}",name_to_res.get(&hh).unwrap());
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
    }
    save_lines(&outfile, results,outfile.ends_with(".gz"));
    //stattest();
}
