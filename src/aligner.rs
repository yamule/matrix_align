use std::collections::HashMap;
use std::collections::HashSet;
use std::vec;

use rand::distributions::weighted::alias_method::Weight;

use self::matrix_process::element_multiply;
use self::matrix_process::vector_add;

use super::*;

#[allow(unused_imports)]
use std::time::Instant;


const NUM_CHARTYPE:usize = 28;
const ACCEPTABLE_CHARS:&str = "ACDEFGHIKLMNPQRSTVWXY-";
const GAP_CHAR:char = '-';
const DIREC_UPLEFT:u8 = 0;
const DIREC_UP:u8 = 1;
const DIREC_LEFT:u8 = 2;

#[derive(Debug,Clone)]
pub struct GMatColumn{
    pub match_vec:Vec<f32>,//match 時に使用されるベクトル
    pub match_ratio:f32,
    pub del_ratio:f32,
    pub connected_ratio:f32, // 前の残基と連続している重み合計
    pub gapped_ratio:f32,// 前の残基と連続していない重み合計
}
impl GMatColumn{
    pub fn new(vecsize:usize,v:Option<Vec<f32>>,exx:Option<(f32,f32,f32,f32)>)->GMatColumn{
        if let Some(x) = exx{
            return GMatColumn{
                match_vec:v.unwrap_or_else(||vec![0.0;vecsize]),
                match_ratio:x.0,
                del_ratio:x.1,
                connected_ratio:x.2,
                gapped_ratio:x.3
            }
        }else{
            return GMatColumn{
                match_vec:v.unwrap_or_else(||vec![0.0;vecsize]),
                match_ratio:1.0,
                del_ratio:0.0,
                connected_ratio:1.0,
                gapped_ratio:0.0
            }
        }
    }
    pub fn set(&mut self,vvec:&Vec<f32>,match_ratio:f32,del_ratio:f32,connected_weight:f32,gapped_weight:f32){
        self.match_vec = vvec.clone();
        self.match_ratio = match_ratio;
        self.del_ratio = del_ratio;
        self.connected_ratio = connected_weight;
        self.gapped_ratio = gapped_weight;
    }
}

pub struct ScoredSeqAligner{
    pub dp_matrix:Vec<Vec<Vec<f32>>>,
    pub path_matrix:Vec<Vec<Vec<u8>>>,
    pub charmap:Vec<usize>,
    pub vec_size:usize,
    pub alen:usize,
    pub blen:usize,
    pub penalty_warning:bool
}
impl ScoredSeqAligner {
    pub fn new(vec_size:usize,buff_len:usize,buff_seqnum:usize)->ScoredSeqAligner{
        let dp_matrix:Vec<Vec<Vec<f32>>> = vec![vec![vec![];1];1];
        let path_matrix:Vec<Vec<Vec<u8>>> = vec![vec![vec![];1];1];
        let mut charmap:Vec<usize> = vec![NUM_CHARTYPE;256];
        let weights:Vec<f32> = vec![1.0;buff_seqnum];
        let cc:Vec<char> = ACCEPTABLE_CHARS.chars().into_iter().collect();
        for ee in cc.into_iter().enumerate(){
            charmap[ee.1 as usize] = ee.0;
        }

        let mut ret = ScoredSeqAligner{
            dp_matrix:dp_matrix
            ,path_matrix:path_matrix
            ,charmap
            ,vec_size:vec_size
            ,alen:0
            ,blen:0
            ,penalty_warning:false
        };
        ret.reconstruct_matrix(buff_len, buff_len);
        return ret;
    }
    
    pub fn reconstruct_matrix(&mut self,amax:usize,bmax:usize){
        self.dp_matrix = vec![vec![vec![0.0;3];bmax+1];amax+1];
        self.path_matrix  = vec![vec![vec![0;3];bmax+1];amax+1];
    }

    pub fn perform_dp(&mut self,a:&ScoredSequence,b:&ScoredSequence,gap_open_penalty:f32,gap_extension_penalty:f32)->(Vec<(i32,i32)>,f32) {
        assert!(gap_extension_penalty <= 0.0);
        assert!(gap_open_penalty <= 0.0);
        if gap_extension_penalty < gap_open_penalty && !self.penalty_warning{
            eprintln!("Gap open penalty is larger than gap extension penalty. open :{}, extension: {}",gap_open_penalty,gap_extension_penalty);
            self.penalty_warning = true;
        }
        let aalen = a.get_alignment_length();
        let bblen = b.get_alignment_length();
        let recflag = self.dp_matrix.len() <= aalen || self.dp_matrix[0].len() <= bblen;
        if recflag{
            self.reconstruct_matrix(aalen+25, bblen+25);
        }
        
        let mut currentpenal:f32;
        
        // B 側 N 末にギャップを入れる
        for ii in 0..=aalen{
            if ii != 0{
                if ii == 1{
                    currentpenal = a.gmat[ii-1].connected_ratio*gap_open_penalty+a.gmat[ii-1].gapped_ratio*gap_extension_penalty;
                }else{
                    currentpenal = a.gmat[ii-1].connected_ratio*gap_extension_penalty+a.gmat[ii-1].gapped_ratio*gap_extension_penalty;
                }
                self.dp_matrix[ii][0][DIREC_LEFT as usize] = self.dp_matrix[ii-1][0][DIREC_LEFT as usize] + currentpenal;
                self.dp_matrix[ii][0][DIREC_UPLEFT as usize] = self.dp_matrix[ii][0][DIREC_LEFT as usize]-1000.0;
                self.dp_matrix[ii][0][DIREC_UP as usize] = self.dp_matrix[ii][0][DIREC_LEFT as usize]-1000.0;
            }
            self.path_matrix[ii][0][0] = DIREC_LEFT;
            self.path_matrix[ii][0][1] = DIREC_LEFT;
            self.path_matrix[ii][0][2] = DIREC_LEFT;
        }

        // A 側 N 末にギャップを入れる
        let mut currentpenal;
        for ii in 0..=bblen{
            if ii != 0{
                if ii == 1{
                    currentpenal = b.gmat[ii-1].connected_ratio*gap_open_penalty+b.gmat[ii-1].gapped_ratio*gap_extension_penalty;
                }else{
                    currentpenal = b.gmat[ii-1].connected_ratio*gap_extension_penalty+b.gmat[ii-1].gapped_ratio*gap_extension_penalty;
                }

                self.dp_matrix[0][ii][DIREC_UP as usize] = self.dp_matrix[0][ii-1][DIREC_UP as usize]+currentpenal;
                self.dp_matrix[0][ii][DIREC_UPLEFT as usize] = self.dp_matrix[0][ii][DIREC_UP as usize]-1000.0;
                self.dp_matrix[0][ii][DIREC_LEFT as usize] = self.dp_matrix[0][ii][DIREC_UP as usize]-1000.0;
            }
            self.path_matrix[0][ii][0] = DIREC_UP;
            self.path_matrix[0][ii][1] = DIREC_UP;
            self.path_matrix[0][ii][2] = DIREC_UP;
        }
        self.dp_matrix[0][0][DIREC_UP as usize] = 0.0;
        self.dp_matrix[0][0][DIREC_UPLEFT as usize] = 0.0;
        self.dp_matrix[0][0][DIREC_LEFT as usize] = 0.0;
        self.path_matrix[0][0][0] = 0;
        self.path_matrix[0][0][1] = 0;
        self.path_matrix[0][0][2] = 0;

        let mut aavec:Vec<&Vec<f32>> = vec![];
        let mut aweight:Vec<f32> = vec![];
        for ii in 0..aalen{
            aavec.push(&a.gmat[ii].match_vec);
            aweight.push(a.gmat[ii].match_ratio);
        }
        
        let mut bbvec:Vec<&Vec<f32>> = vec![];
        let mut bweight:Vec<f32> = vec![];
        for ii in 0..bblen{
            bbvec.push(&b.gmat[ii].match_vec);
            bweight.push(b.gmat[ii].match_ratio);
        }
        //バッファに入れようかと思ったが、結局新しく領域を確保していたのでやめた
        let match_score:Vec<Vec<f32>> = gmat::calc_dist_zscore_matrix(&aavec, &bbvec,None,None);
        //let match_score:Vec<Vec<f32>> = gmat::calc_dot_product_matrix(&aavec, &bbvec,Some(&aweight),Some(&bweight));
        
        /*
        println!("#{:?} vs {:?}",a.headers,b.headers);
        for ii in 0..aalen{
            for jj in 0..bblen{
                print!("{}\t",match_score[ii][jj]);
            }
            println!("");
        }
        println!("===");
        */
        
        for ii in 1..=aalen{
            for jj in 1..=bblen{
                let acol = &a.gmat[ii-1];
                let bcol = &b.gmat[jj-1];
                //let sc:f32 = ScoredSeqAligner::calc_match_score(&acol.0,&bcol.0);
                let sc:f32 = match_score[ii-1][jj-1];
                

                let diag_m:f32 = self.dp_matrix[ii-1][jj-1][DIREC_UPLEFT as usize] + sc;
                let diag_l:f32 = self.dp_matrix[ii-1][jj-1][DIREC_LEFT as usize] + sc;
                let diag_u:f32 = self.dp_matrix[ii-1][jj-1][DIREC_UP as usize] + sc;


                let lef_m:f32 = self.dp_matrix[ii-1][jj][DIREC_UPLEFT as usize] + acol.connected_ratio*gap_open_penalty + acol.gapped_ratio*gap_extension_penalty;
                let lef_l:f32 = self.dp_matrix[ii-1][jj][DIREC_LEFT as usize] + acol.connected_ratio*gap_extension_penalty + acol.gapped_ratio*gap_extension_penalty;
                let lef_u:f32 = self.dp_matrix[ii-1][jj][DIREC_UP as usize] + acol.connected_ratio*gap_open_penalty + acol.gapped_ratio*gap_extension_penalty;


                let up_m:f32 = self.dp_matrix[ii][jj-1][DIREC_UPLEFT as usize] + bcol.connected_ratio*gap_open_penalty + bcol.gapped_ratio*gap_extension_penalty;
                let up_l:f32 = self.dp_matrix[ii][jj-1][DIREC_LEFT as usize] + bcol.connected_ratio*gap_open_penalty + bcol.gapped_ratio*gap_extension_penalty;
                let up_u:f32 = self.dp_matrix[ii][jj-1][DIREC_UP as usize] + bcol.connected_ratio*gap_extension_penalty + bcol.gapped_ratio*gap_extension_penalty;
                

                let px = vec![
                    (DIREC_UPLEFT,(diag_m,diag_l,diag_u)),
                    (DIREC_LEFT,(lef_m,lef_l,lef_u)),
                    (DIREC_UP,(up_m,up_l,up_u)),
                ];
                //println!("{} {} {}",diag,leff,upp);
                for pp in px.iter(){
                    let poss = pp.0;
                    let (_m,_l,_u) = pp.1;
                    if _m >= _l && _m >= _u{
                        self.dp_matrix[ii][jj][poss as usize] = _m;
                        self.path_matrix[ii][jj][poss as usize] = DIREC_UPLEFT;
                    }else{
                        if _l >= _u{
                            self.dp_matrix[ii][jj][poss as usize] = _l;
                            self.path_matrix[ii][jj][poss as usize] = DIREC_LEFT;
                        }else{
                            self.dp_matrix[ii][jj][poss as usize] = _u;
                            self.path_matrix[ii][jj][poss as usize] = DIREC_UP;
                        }
                    }
                }
            }
        }
        //panic!("{:?}",self.dp_matrix);
        let mut currentx = aalen;
        let mut currenty = bblen;
        let mut currentpos = DIREC_UPLEFT;
        let mut maxscore = self.dp_matrix[currentx][currenty][currentpos as usize];
        for ii in 0..3{
            if maxscore < self.dp_matrix[currentx][currenty][ii]{
                currentpos = ii as u8;
                maxscore = self.dp_matrix[currentx][currenty][ii];
            }
        }
        let mut nexpos = self.path_matrix[currentx][currenty][currentpos as usize];
        let mut aligned_tuple:Vec<(i32,i32)> = vec![];
        while currentx > 0 || currenty > 0{
            if currentpos == DIREC_UPLEFT{
                currentx -= 1;
                currenty -= 1;
                aligned_tuple.push((currentx as i32,currenty as i32));
            }else if currentpos == DIREC_UP{
                currenty -= 1;
                aligned_tuple.push((-1,currenty as i32));
            }else if currentpos == DIREC_LEFT{
                currentx -= 1;
                aligned_tuple.push((currentx as i32,-1));
            }else{
                panic!("???");
            }
            currentpos = nexpos;
            nexpos = self.path_matrix[currentx][currenty][currentpos as usize];
            
        }
        //println!("{:?}",self.dp_matrix);
        aligned_tuple.reverse();
        return (aligned_tuple,maxscore);
    }

    pub fn make_alignment(&mut self
        ,mut a:ScoredSequence
        ,mut b:ScoredSequence
        ,alignment:Vec<(i32,i32)>
        ,profile_only:bool // true にすると alignment の文字は a に関してのみ保持する
    )->ScoredSequence{
        assert!(!profile_only);//後で追加する
        let anumaliseq = a.get_num_seq();
        let bnumaliseq = b.get_num_seq();
        let numallseq = anumaliseq+bnumaliseq;
        let vec_size = a.get_vec_size();
        let alignment_length = alignment.len();

        let mut new_alignments:Vec<Vec<char>> = vec![vec![];numallseq];
        let mut gapper:Vec<Vec<char>> = vec![vec![],vec![]];//全体ギャップ計算
        for aa in 0..alignment.len(){
            let ppos = alignment[aa];
            if ppos.0 > -1{
                let poss = ppos.0 as usize;
                for ap in 0..anumaliseq{
                    new_alignments[ap].push(a.alignments[ap][poss as usize]);
                }
                gapper[0].push('X');
            }else{
                for ap in 0..anumaliseq{
                    new_alignments[ap].push('-');
                }
                gapper[0].push('-');
            }
            if ppos.1 > -1{
                let poss = ppos.1 as usize;
                for ap in 0..bnumaliseq{
                    new_alignments[ap+anumaliseq].push(b.alignments[ap][poss as usize]);
                }
                gapper[1].push('X');
            }else{
                for ap in 0..bnumaliseq{
                    new_alignments[ap+anumaliseq].push('-');
                }
                gapper[1].push('-');
            }
        }
        let mut headers:Vec<String> = vec![];
        headers.append(&mut a.headers);
        headers.append(&mut b.headers);

        let mut ret:ScoredSequence = ScoredSequence::new(headers.into_iter().zip(new_alignments.into_iter()).collect()
        , vec_size,None,None,None);

        let aweight = a.get_weight_sum();
        let bweight = b.get_weight_sum();
        let mut ex_weights:Vec<(f32,f32,f32,f32)> = vec![(0.0,0.0,0.0,0.0);alignment_length+1];

        for (wei,alichar,sprof) in vec![(aweight,&gapper[0],&a),(bweight,&gapper[1],&b)]{
            let mut poscount = 0_usize;
            for alipos in 0..alignment_length{
                let mut sum_weight = 0.0;
                let mut sum_weight_del = 0.0;
                let mut ungapratio = 0.0;
                let mut gapratio = 0.0;
                // 前後で繋がっていて GAPOPEN が必要なもののウェイトを取る
                if alipos == 0{
                    if alichar[alipos] != GAP_CHAR{
                        ungapratio += wei*sprof.gmat[poscount].connected_ratio;
                        gapratio += wei*sprof.gmat[poscount].gapped_ratio;
                    }else{
                        gapratio += wei*sprof.gmat[poscount].gapped_ratio+ wei*sprof.gmat[poscount].connected_ratio;
                    }
                }else{
                    if alichar[alipos-1] != GAP_CHAR && alichar[alipos] != GAP_CHAR {
                        ungapratio += wei*sprof.gmat[poscount].connected_ratio;
                        gapratio += wei*sprof.gmat[poscount].gapped_ratio;
                    }else{
                        gapratio += wei*sprof.gmat[poscount].gapped_ratio+ wei*sprof.gmat[poscount].connected_ratio;
                    }
                }

                if alichar[alipos] != GAP_CHAR{
                    for vv in 0..vec_size{
                        ret.gmat[alipos].match_vec[vv] += sprof.gmat[poscount].match_vec[vv]
                        *wei*sprof.gmat[poscount].match_ratio;
                    }
                    sum_weight += wei*sprof.gmat[poscount].match_ratio;
                    sum_weight_del += wei*sprof.gmat[poscount].del_ratio;
                    poscount += 1;
                }else{
                    sum_weight_del += wei;
                }
                ex_weights[alipos].0 += sum_weight;
                ex_weights[alipos].1 += sum_weight_del;
                ex_weights[alipos].2 += ungapratio;
                ex_weights[alipos].3 += gapratio;
            }
        }
        for alipos in 0..alignment_length{
            
            let sum_weight = ex_weights[alipos].0;
            let sum_weight_del = ex_weights[alipos].1;
            let ungapratio = ex_weights[alipos].2;
            let gapratio = ex_weights[alipos].3;

            assert!((ungapratio+gapratio) > 0.0);

            ret.gmat[alipos].connected_ratio = ungapratio/(ungapratio+gapratio);
            ret.gmat[alipos].gapped_ratio = gapratio/(ungapratio+gapratio);

            if sum_weight > 0.0{
                unsafe{
                    element_multiply(&mut ret.gmat[alipos].match_vec,1.0/sum_weight);
                }
            }
            
            assert!(sum_weight+sum_weight_del > 0.0);

            ret.gmat[alipos].match_ratio = sum_weight/(sum_weight+sum_weight_del);
            ret.gmat[alipos].del_ratio = sum_weight_del/(sum_weight+sum_weight_del);
        }

        // 最後の残基以降のギャップ
        let mut ungapratio = 0.0;
        let mut gapratio = 0.0;
        for (wei,alichar,sprof) in vec![(aweight,&gapper[0],&a),(bweight,&gapper[1],&b)]{
            if alichar[alignment_length-1] != GAP_CHAR{
                ungapratio += wei*sprof.gmat[sprof.gmat.len()-1].connected_ratio;
                gapratio += wei*sprof.gmat[sprof.gmat.len()-1].gapped_ratio;
            }else{
                gapratio += wei*sprof.gmat[sprof.gmat.len()-1].gapped_ratio+ wei*sprof.gmat[sprof.gmat.len()-1].connected_ratio;
            }
        }
        
        assert!(ungapratio+gapratio > 0.0);
        ret.gmat[alignment_length].connected_ratio = ungapratio/(ungapratio+gapratio);
        ret.gmat[alignment_length].gapped_ratio = gapratio/(ungapratio+gapratio);
        ret.seq_weights.clear();
        ret.seq_weights.append(&mut a.seq_weights);
        ret.seq_weights.append(&mut b.seq_weights);

        return ret;
    }
    pub fn make_msa(&mut self,mut sequences: Vec<ScoredSequence>,gap_open_penalty:f32,gap_extension_penalty:f32,profile_only:bool)
    -> (Vec<ScoredSequence>,f32){
        let mut final_score:f32 = 0.0;
        sequences.reverse();
        let mut center_seq = sequences.pop().unwrap();
        let mut firstrun = true;
        while sequences.len() > 0{
            let bseq = sequences.pop().unwrap();
            let dpres;
            let mut newgroup;
            if firstrun{
                dpres = self.perform_dp(&center_seq,&bseq,gap_open_penalty,gap_extension_penalty);
                newgroup = ScoredSeqAligner::make_alignment(self,center_seq,bseq,dpres.0,profile_only);
                firstrun = false;
            }else{
                dpres = self.perform_dp(&bseq,&center_seq,gap_open_penalty,gap_extension_penalty);
                newgroup = ScoredSeqAligner::make_alignment(self,bseq,center_seq,dpres.0,profile_only);
            };
            final_score = dpres.1;
            
            center_seq = newgroup;
        }
        return (vec![center_seq],final_score);
    }

    //配列のウェイトを計算する
    //全配列の MSA が作られている前提
    pub fn calc_weights(alignments:&Vec<Vec<char>>)->Vec<f32>{
        return sequence_weighting::calc_henikoff_henikoff_weight(&alignments);
    }

}

pub struct ScoredSequence{
    pub gmat:Vec<GMatColumn>,
    pub alignments:Vec<Vec<char>>,
    pub headers:Vec<String>,//Fasta Header 行
    pub seq_weights:Vec<f32>
}

impl ScoredSequence{
    pub fn new(alignments_:Vec<(String,Vec<char>)>,vec_size:usize,seq_weights_:Option<Vec<f32>>
        ,gmat_:Option<Vec<Vec<f32>>>,gap_state_:Option<Vec<(f32,f32,f32,f32)>>)-> ScoredSequence{
        let seqnum = alignments_.len();
        let mut gmat:Vec<GMatColumn> = vec![];
        if let Some(mut x) = gmat_{
            if let Some(y) = gap_state_{
                assert!(x.len()+1 == y.len());
                x.push(vec![0.0;vec_size]);// ギャップ情報だけあるカラムがあるので
                for xx in x.into_iter().zip(y.into_iter()){
                    gmat.push(GMatColumn::new(vec_size,Some(xx.0),Some(xx.1)));
                }
            }else{
                for xx in x.into_iter(){
                    gmat.push(GMatColumn::new(vec_size,Some(xx),None));
                }
                gmat.push(GMatColumn::new(vec_size,None,None));
            }
        }else{
            for xx in 0..alignments_[0].1.len(){
                gmat.push(GMatColumn::new(vec_size,None,None));
            }
            gmat.push(GMatColumn::new(vec_size,None,None));// ギャップ情報だけあるカラム
        }

        let mut alignments:Vec<Vec<char>> = vec![];
        let mut headers:Vec<String> = vec![];
        for aa in alignments_.into_iter(){
            alignments.push(aa.1);
            headers.push(aa.0);
        }
        return ScoredSequence{
            alignments:alignments,
            headers:headers,
            gmat:gmat,
            seq_weights:seq_weights_.unwrap_or_else(||vec![1.0;seqnum])
        }
    }

    pub fn gmat_str(&self)->Vec<String>{
        let mut ret:Vec<String> = vec![];
        for ii in 0..self.gmat.len(){
            let gg = &self.gmat[ii];
            ret.push(format!("@\t{}\t{}\t{}\t{}",gg.match_ratio,gg.del_ratio,gg.connected_ratio,gg.gapped_ratio));
            let mut pret:Vec<String> = vec![];
            if ii < self.gmat.len()-1{
                pret.push(format!("{}",self.alignments[0][ii]));
                for ff in self.gmat[ii].match_vec.iter(){
                    pret.push(format!("{}",ff));
                }
                ret.push(pret.join("\t"));
            }
        }
        return ret;
    }

    pub fn create_merged_msa(&self)->ScoredSequence{
        let alignments:Vec<Vec<char>> = vec![vec!['X';self.alignments[0].len()]];
        let headers:Vec<String> = vec!["dummy".to_owned()];
        let gmat = self.gmat.clone();
        let seq_weights = vec![self.get_weight_sum()];
        return ScoredSequence{
            alignments:alignments,
            headers:headers,
            gmat:gmat,
            seq_weights:seq_weights
        };
    }

    pub fn get_alignment_length(&self)->usize{// 使い回す方針を取るならこの辺変える
        return self.alignments[0].len();
    }

    pub fn get_num_seq(&self)-> usize{
        return self.alignments.len();
    }

    pub fn get_vec_size(&self)-> usize{
        return self.gmat[0].match_vec.len();
    }

    pub fn get_weight_sum(&self)->f32{
        return self.seq_weights.iter().sum();
    }
    /*
    pub fn split_one_with_id(&mut self,id:i32,aligner:&ScoredSeqAligner)->ScoredSequence{
        let mut targetidx:usize = self.get_index_of(id);
        return self.split_one(targetidx, aligner);
    }
    */
    
}
