use std::collections::HashMap;
use std::collections::HashSet;

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
    pub match_weight:f32,
    pub del_weight:f32,
    pub connected_weight:f32, // 前の残基と連続している重み合計
    pub gapped_weight:f32,// 前の残基と連続していない重み合計
}
impl GMatColumn{
    pub fn new(vecsize:usize)->GMatColumn{
        return GMatColumn{
            match_vec:vec![0.0;vecsize],
            match_weight:1.0,
            del_weight:0.0,
            connected_weight:1.0,
            gapped_weight:0.0
        }
    }
    pub fn set(&mut self,vvec:&Vec<f32>,match_weight:f32,del_weight:f32,connected_weight:f32,gapped_weight:f32){
        self.match_vec = vvec.clone();
        self.match_weight = match_weight;
        self.del_weight = del_weight;
        self.connected_weight = connected_weight;
        self.gapped_weight = gapped_weight;
    }
}

pub struct ScoredSeqAligner{
    pub dp_matrix:Vec<Vec<Vec<f32>>>,
    pub path_matrix:Vec<Vec<Vec<u8>>>,
    pub charmap:Vec<usize>,
    pub vec_size:usize,
    pub alen:usize,
    pub blen:usize,
    pub weights:Vec<f32>,
    pub alignment_buffer:Vec<Vec<char>>,//全配列がアラインメントされた状態で保持されている
    pub gmat_buffer:Vec<Vec<GMatColumn>>,
    pub alibuff_used:Vec<bool>,
    pub gmat_used:Vec<bool>,
    pub alibuff_next:usize,
    pub gmat_next:usize,
    pub penalty_warning:bool
}
impl ScoredSeqAligner {
    pub fn new(vec_size:usize,buff_len:usize,buff_seqnum:usize)->ScoredSeqAligner{
        let dp_matrix:Vec<Vec<Vec<f32>>> = vec![vec![vec![];1];1];
        let path_matrix:Vec<Vec<Vec<u8>>> = vec![vec![vec![];1];1];
        let mut charmap:Vec<usize> = vec![NUM_CHARTYPE;256];
        let alignment_buffer = vec![vec!['-';buff_len];buff_seqnum];
        let gmat_buffer:Vec<Vec<GMatColumn>> = vec![vec![GMatColumn::new(vec_size);buff_len+1];buff_seqnum];
        let buff_used:Vec<bool> = vec![false;buff_seqnum];
        let gmat_used:Vec<bool> = vec![false;buff_seqnum];
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
            ,weights:weights
            ,alignment_buffer:alignment_buffer
            ,alibuff_used:buff_used
            ,gmat_buffer
            ,gmat_used
            ,alibuff_next:0
            ,gmat_next:0
            ,penalty_warning:false
        };
        ret.reconstruct_matrix(buff_len, buff_len);
        return ret;
    }
    
    
    pub fn get_unused_alibuffid(&mut self)->usize{
        for ii in self.alibuff_next..self.alibuff_used.len(){
            if !self.alibuff_used[ii]{
                self.alibuff_next = ii+1;
                return ii;
            }
        }
        for ii in 0..self.alibuff_next{
            if !self.alibuff_used[ii]{
                self.alibuff_next = ii+1;
                return ii;
            }
        }
        let newid = self.alibuff_used.len();
        self.alibuff_used.push(false);
        self.alignment_buffer.push(vec!['-';self.alignment_buffer[0].len()]);
        self.alibuff_next = self.alibuff_used.len();
        return newid;
    }
    
    pub fn register_alibuff(&mut self,id:usize,len:usize){
        assert!(!self.alibuff_used[id]);
        self.check_alibufflen(id, len);
        self.alibuff_used[id] = true;
    }

    pub fn release_alibuff(&mut self,id:usize){
        assert!(self.alibuff_used[id]);
        self.alibuff_used[id] = false;
    }

    pub fn check_alibufflen(&mut self,id:usize,len:usize){
        if self.alignment_buffer[id].len() < len{
            self.alignment_buffer[id] = vec!['-';len+5];
        }
    }

    pub fn get_unused_gmatbuffid(&mut self)->usize{
        for ii in self.gmat_next..self.gmat_used.len(){
            if !self.gmat_used[ii]{
                self.gmat_next = ii+1;
                return ii;
            }
        }
        for ii in 0..self.gmat_next{
            if !self.gmat_used[ii]{
                self.gmat_next = ii+1;
                return ii;
            }
        }
        let newid = self.gmat_used.len();
        self.gmat_used.push(false);
        self.gmat_buffer.push(vec![GMatColumn::new(self.vec_size);self.gmat_buffer[0].len()]);
        return newid;
    }
    
    pub fn register_gmatbuff(&mut self,id:usize,len:usize){
        assert!(!self.gmat_used[id]);
        self.format_gmatbuff(id, len);
        self.check_gmatbuff_length(id, len);
        self.gmat_used[id] = true;
    }

    pub fn copy_gmat_a_to_b(&mut self,a:usize,b:usize,len:usize){
        assert!(self.gmat_buffer[a].len() >= len);
        self.check_gmatbuff_length(b,len);
        self.gmat_buffer[b] = self.gmat_buffer[a].clone();
    }
    
    pub fn copy_ali_a_to_b(&mut self,a:usize,b:usize,len:usize){
        assert!(self.alignment_buffer[a].len() >= len);
        self.check_alibufflen(b,len);

        self.alignment_buffer[b] = self.alignment_buffer[a].clone();

    }

    pub fn release_gmatbuff(&mut self,id:usize){
        assert!(self.gmat_used[id]);
        self.gmat_used[id] = false;
    }

    //文字列を保持しておくバッファのサイズを確認し、小さい場合は冗長性を少し加えて大きくする
    pub fn check_gmatbuff_length(&mut self,id:usize,len:usize){
        if self.gmat_buffer[id].len() < len+1{
            self.gmat_buffer[id] = vec![GMatColumn::new(self.vec_size);len+6];
        }
    }
    
    pub fn format_gmatbuff(&mut self,id:usize,len:usize){
        self.check_gmatbuff_length(id, len);
        for vv in 0..=len{
            self.gmat_buffer[id][vv].match_vec.fill(0.0);
            self.gmat_buffer[id][vv].del_weight = 0.0;
            self.gmat_buffer[id][vv].del_weight = 0.0;
            self.gmat_buffer[id][vv].del_weight = 0.0;
        }
    }

    //呼び出し元から直接 dot_product に飛ばしてもいいかも
    pub fn calc_match_score(a:&Vec<f32>,b:&Vec<f32>)->f32{
        unsafe{
            return matrix_process::dot_product(a,b);
        }
    }
    
    pub fn gmat_colget(&self,lid:usize,pos:usize)->&GMatColumn{
        return &self.gmat_buffer[lid][pos];
    }

    //別関数になるかも
    pub fn gmat_set(&mut self,lid:usize,pos:usize,vvec:&Vec<f32>,match_weight:f32,del_weight:f32,connected_weight:f32,gapped_weight:f32){
        self.gmat_buffer[lid][pos].set(vvec,match_weight,del_weight,connected_weight,gapped_weight);
    }


    //この辺違う関数になるはず
    /*
    pub fn gmat_add(&mut self,lid:usize,pos:usize,c:char,num:i32){
        self.gmat_buffer[lid][pos][self.charmap[c as usize]] += num;
        assert!(self.gmat_buffer[lid][pos][self.charmap[c as usize]] >= 0);
    }

    pub fn gmat_add_i(&mut self,lid:usize,pos:usize,c:usize,num:i32){
        self.gmat_buffer[lid][pos][c] += num;
        assert!(self.gmat_buffer[lid][pos][c] >= 0);
    }
    */

    pub fn ali_set(&mut self,aid:usize,pos:usize,c:char){
        self.alignment_buffer[aid][pos] = c;
    }

    pub fn ali_get(&self,aid:usize,pos:usize)->char{
        return self.alignment_buffer[aid][pos];
    }
    pub fn ali_get_i(&self,aid:usize,pos:usize)->usize{
        return self.charmap[self.alignment_buffer[aid][pos] as usize];
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
        let aalen = a.alignment_length;
        let bblen = b.alignment_length;
        let recflag = self.dp_matrix.len() <= aalen || self.dp_matrix[0].len() <= bblen;
        if recflag{
            self.reconstruct_matrix(aalen+25, bblen+25);
        }
        
        let alid = a.gmatbuff_id as usize;
        let blid = b.gmatbuff_id as usize;
        
        let mut currentpenal:f32;
        
        // B 側 N 末にギャップを入れる
        for ii in 0..=aalen{
            if ii != 0{
                if ii == 1{
                    currentpenal = self.gmat_colget(alid,ii-1).connected_weight*gap_open_penalty+self.gmat_colget(alid,ii-1).gapped_weight*gap_extension_penalty;
                }else{
                    currentpenal = self.gmat_colget(alid,ii-1).connected_weight*gap_extension_penalty+self.gmat_colget(alid,ii-1).gapped_weight*gap_extension_penalty;
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
                    currentpenal = self.gmat_colget(blid,ii-1).connected_weight*gap_open_penalty+self.gmat_colget(blid,ii-1).gapped_weight*gap_extension_penalty;
                }else{
                    currentpenal = self.gmat_colget(blid,ii-1).connected_weight*gap_extension_penalty+self.gmat_colget(blid,ii-1).gapped_weight*gap_extension_penalty;
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
            aavec.push(&self.gmat_buffer[alid][ii].match_vec);
            aweight.push(self.gmat_buffer[alid][ii].match_weight);
        }
        
        let mut bbvec:Vec<&Vec<f32>> = vec![];
        let mut bweight:Vec<f32> = vec![];
        for ii in 0..bblen{
            bbvec.push(&self.gmat_buffer[blid][ii].match_vec);
            bweight.push(self.gmat_buffer[blid][ii].match_weight);
        }
        //バッファに入れようかと思ったが、結局新しく領域を確保していたのでやめた
        let match_score:Vec<Vec<f32>> = gmat::calc_dist_zscore_matrix(&aavec, &bbvec,Some(&aweight),Some(&bweight));
        //let match_score:Vec<Vec<f32>> = gmat::calc_dist_zscore_matrix(&aavec, &bbvec,None,None);

        for ii in 1..=aalen{
            for jj in 1..=bblen{
                let acol = self.gmat_colget(alid,ii-1);
                let bcol = self.gmat_colget(blid,jj-1);
                //let sc:f32 = ScoredSeqAligner::calc_match_score(&acol.0,&bcol.0);
                let sc:f32 = match_score[ii-1][jj-1];
                

                let diag_m:f32 = self.dp_matrix[ii-1][jj-1][DIREC_UPLEFT as usize] + sc;
                let diag_l:f32 = self.dp_matrix[ii-1][jj-1][DIREC_LEFT as usize] + sc;
                let diag_u:f32 = self.dp_matrix[ii-1][jj-1][DIREC_UP as usize] + sc;


                let lef_m:f32 = self.dp_matrix[ii-1][jj][DIREC_UPLEFT as usize] + acol.connected_weight*gap_open_penalty + acol.gapped_weight*gap_extension_penalty;
                let lef_l:f32 = self.dp_matrix[ii-1][jj][DIREC_LEFT as usize] + acol.connected_weight*gap_extension_penalty + acol.gapped_weight*gap_extension_penalty;
                let lef_u:f32 = self.dp_matrix[ii-1][jj][DIREC_UP as usize] + acol.connected_weight*gap_open_penalty + acol.gapped_weight*gap_extension_penalty;


                let up_m:f32 = self.dp_matrix[ii][jj-1][DIREC_UPLEFT as usize] + bcol.connected_weight*gap_open_penalty + bcol.gapped_weight*gap_extension_penalty;
                let up_l:f32 = self.dp_matrix[ii][jj-1][DIREC_LEFT as usize] + bcol.connected_weight*gap_open_penalty + bcol.gapped_weight*gap_extension_penalty;
                let up_u:f32 = self.dp_matrix[ii][jj-1][DIREC_UP as usize] + bcol.connected_weight*gap_extension_penalty + bcol.gapped_weight*gap_extension_penalty;
                

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
        aligned_tuple.reverse();
        return (aligned_tuple,maxscore);
    }

    pub fn make_alignment(&mut self
        ,mut a:ScoredSequence
        ,mut b:ScoredSequence
        ,alignment:Vec<(i32,i32)>
        ,profile_only:bool // true にすると alignment の文字は a に関してのみ保持する
    )->ScoredSequence{

        assert!(!profile_only);//全体の Profile を iterative に作成して、その間はアラインメントの更新は行わず、最終的な MSA とする際だけ更新する予定だが、今は全体 Alignment を作らないと計算できない
        // 後で変える
        if a.id >= b.id{//なぜこうしているのか忘れてしまった・・・
            panic!("a.id must be smaller than b.id! {} {} ",a.id,b.id);
        }

        let anumaliseq = a.primary_ids.len();//alignment にだけ使う
        let bnumaliseq = b.primary_ids.len();
        let numallseq = anumaliseq+bnumaliseq;
        let alignment_length = alignment.len();

        let mut new_aids = vec![];

        if !profile_only{
            for _ii in 0..numallseq{
                let aid = self.get_unused_alibuffid();
                new_aids.push(aid);
                self.register_alibuff(aid,alignment_length);
            }
        }else{
            for _ii in 0..anumaliseq{
                let aid = self.get_unused_alibuffid();
                new_aids.push(aid);
                self.register_alibuff(aid,alignment_length);
            }
        }


        for aa in 0..alignment.len(){
            let ppos = alignment[aa];
            if ppos.0 > -1{
                let poss = ppos.0 as usize;
                //所属する配列について、アラインメントされた状態の配列を新しいバッファ上に構築する
                for ap in 0..anumaliseq{
                    self.ali_set(new_aids[ap],aa, self.ali_get(a.alibuff_idx[ap], poss as usize));
                }
            }else{
                //所属する配列について、アラインメントされた状態の配列を新しいバッファ上に構築する
                for ap in 0..anumaliseq{
                    self.ali_set(new_aids[ap],aa,'-');
                }
            }
            if ppos.1 > -1{
                let poss = ppos.1 as usize;
                if !profile_only{
                    //所属する配列について、アラインメントされた状態の配列を新しいバッファ上に構築する
                    for ap in 0..bnumaliseq{
                        self.ali_set(new_aids[ap+anumaliseq],aa,
                            self.ali_get(b.alibuff_idx[ap], poss as usize));
                    }
                }
            }else{
                if !profile_only{
                    //所属する配列について、アラインメントされた状態の配列を新しいバッファ上に構築する
                    for ap in 0..bnumaliseq{
                        self.ali_set(new_aids[ap+anumaliseq],aa,'-');
                    }
                }
            }
        }
        let mut primary_ids:Vec<i32> = vec![];
        
        for ap in 0..anumaliseq{
            primary_ids.push(a.primary_ids[ap]);
        }
        if !profile_only{
            for ap in 0..bnumaliseq{
                primary_ids.push(b.primary_ids[ap]);
            }
        }
        let sumweight = a.weight_sum+b.weight_sum;

        

        a.release_buffs(self);
        b.release_buffs(self);

        let lid = self.get_unused_gmatbuffid();
        self.register_gmatbuff(lid,alignment_length);
        self.format_gmatbuff(lid,alignment_length);
        let mut ret =  ScoredSequence{
            id:-1,
            num_sequences:a.num_sequences+b.num_sequences,
            primary_ids:primary_ids,
            gmatbuff_id:lid as i32,
            alibuff_idx:new_aids,
            alignment_length:alignment_length,
            weight_sum:sumweight
        };
        ret.calc_gmat( self);
        return ret;
    }
    pub fn make_msa(&mut self,mut sequences: Vec<ScoredSequence>,gap_open_penalty:f32,gap_extension_penalty:f32,profile_only:bool)
    -> (Vec<ScoredSequence>,f32){
        let mut next_id = sequences.iter().fold(0,|s,a| s.max(a.id)) as i32 + 10000;
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
            newgroup.set_id(next_id);
            next_id += 1;
            center_seq = newgroup;
        }
        return (vec![center_seq],final_score);
    }

    

    //primid_map は、ScoredSequence 内の primary_ids → distmat の index へのマップ
    //複数 Seq は 全 vs 全 の平均
    pub fn make_msa_with_distmat(&mut self, mut sequences:Vec<ScoredSequence>,distmat:Vec<Vec<f32>>,primid_map:HashMap<i32,usize>,gap_open_penalty:f32,gap_extension_penalty:f32)
    ->(Vec<ScoredSequence>,f32){
        let mut next_id = sequences.iter().fold(0,|s,a| s.max(a.id)) as i32 + 10000;

        while sequences.len() > 1{
            let mut nearest_pair_index = (-1,-1);
            let numseq = sequences.len();
            let mut mindist:f32 = 0.0;
            for ii in 0..numseq{
                for jj in (ii+1)..numseq{
                    let mut sumdist = 0.0;
                    let mut sumcount = 0.0;
                    for pp in sequences[ii].primary_ids.iter(){
                        for qq in sequences[jj].primary_ids.iter(){
                            sumdist += distmat[*primid_map.get(pp).unwrap()][*primid_map.get(qq).unwrap()];
                            sumcount += 1.0;
                        }
                    }
                    if sumcount != 0.0{
                        sumdist /= sumcount;
                    }
                    if nearest_pair_index.0 == -1 || sumdist < mindist{
                        nearest_pair_index = (ii as i32,jj as i32);
                        mindist = sumdist;
                    }
                }
            }
            assert!(nearest_pair_index.0 > -1);
            let mut nearest_pair:(i32,i32) = (sequences[nearest_pair_index.0 as usize].id,sequences[nearest_pair_index.1 as usize].id);
            if nearest_pair.0 > nearest_pair.1{
                let pe = nearest_pair_index.0;
                nearest_pair_index.0 = nearest_pair_index.1;
                nearest_pair_index.1 = pe;
                
                let pe = nearest_pair.0;
                nearest_pair.0 = nearest_pair.1;
                nearest_pair.1 =pe;
            }
            let dpres =  self.perform_dp(&sequences[nearest_pair_index.0 as usize],&sequences[nearest_pair_index.1 as usize]
                ,gap_open_penalty,gap_extension_penalty);

            let mut aseq:Option<ScoredSequence> = None;
            let mut bseq:Option<ScoredSequence> = None;
            for ii in 0..numseq{
                if sequences[ii].id == nearest_pair.0{
                    aseq = Some(sequences.remove(ii));
                    if let Some(_) = bseq{
                        break;
                    }
                    if sequences[ii].id == nearest_pair.1{
                        bseq = Some(sequences.remove(ii));
                        break;
                    }
                }
                if sequences[ii].id == nearest_pair.1{
                    bseq = Some(sequences.remove(ii));
                    if let Some(_) = aseq{
                        break;
                    }
                    if sequences[ii].id == nearest_pair.0{
                        aseq = Some(sequences.remove(ii));
                        break;
                    }
                }
            }

            let mut newgroup = ScoredSeqAligner::make_alignment(self,aseq.unwrap(),bseq.unwrap(),dpres.0,false);
            newgroup.set_id(next_id);
            next_id += 1;
            sequences.insert(0,newgroup);
        }
        return (sequences,-1.0);
    }


    //配列のウェイトを計算する
    //全配列の MSA が作られている前提
    pub fn calc_weights(&mut self,alibuff_ids:&Vec<usize>,alilen:usize){
        for aa in alibuff_ids.iter(){
            if alilen < self.alignment_buffer.len()+1{
                assert!(self.alignment_buffer[*aa][alilen+1] == '-',"All sequences must have been aligned. Found {} at sequence {}.",self.alignment_buffer[*aa][alilen+1],*aa);
            }
        }
        let weights = sequence_weighting::calc_henikoff_henikoff_weight(&self.alignment_buffer,alibuff_ids,alilen);
        for aa in alibuff_ids.iter().enumerate(){
            self.weights[*aa.1] = weights[aa.0];
        }
    }

    pub fn calc_alignment_dist(&self,seq:&ScoredSequence)->(Vec<Vec<f32>>,HashMap<i32,usize>){
        let mut pids = seq.primary_ids.clone();
        pids.sort();
        let primaryid_map:HashMap<i32,usize> = pids.iter().enumerate().map(|m|(*m.1,m.0)).collect();
        assert_eq!(primaryid_map.len(),pids.len());
        
        let mut dist:Vec<Vec<f32>> = vec![vec![0.0;pids.len()];pids.len()];
        for ii in 0..seq.alibuff_idx.len(){
            for jj in (ii+1)..seq.alibuff_idx.len(){
                let mut diff:i32 = 0;
                for kk in 0..seq.alibuff_idx.len(){
                    if self.ali_get(seq.alibuff_idx[ii],kk)
                        != self.ali_get(seq.alibuff_idx[jj],kk){
                        diff += 1;
                    }
                }
                dist[ii][jj] = diff as f32;
            }
        }

        return (dist,primaryid_map);
    }
    
}

pub struct ScoredSequence{
    pub id:i32,
    //pub letters:Vec<Vec<usize>>,//カラムごとの文字の出現数
    //pub alignments:Vec<Vec<char>>,
    pub num_sequences:usize,
    pub primary_ids:Vec<i32>, //weights のインデクスにも対応している
    pub gmatbuff_id:i32,
    pub alibuff_idx:Vec<usize>,
    pub alignment_length:usize,
    pub weight_sum:f32
}

impl ScoredSequence{
    pub fn new(a:Vec<char>,a_gmat:Vec<Vec<f32>>,gap_state:Option<Vec<(f32,f32,f32,f32)>>,aligner:&mut ScoredSeqAligner,register_id:bool)-> ScoredSequence{
        let lid = aligner.get_unused_gmatbuffid();
        let aid = aligner.get_unused_alibuffid();
        aligner.register_gmatbuff(lid,a.len());
        aligner.register_alibuff(aid,a.len());
        let alen = a.len();
        
        for aa in a.into_iter().zip(a_gmat.into_iter()).enumerate(){
            if aligner.charmap[(aa.1).0 as usize] != NUM_CHARTYPE{
                aligner.gmat_set(lid,aa.0,&(aa.1).1,1.0,0.0
                ,1.0,0.0); 
                aligner.ali_set(aid,aa.0,(aa.1).0);
            }else{
                panic!("can not process {}!",(aa.1).0);
            }
        }
        if let Some(x) = gap_state{
            assert_eq!(alen+1,x.len());
            for ii in 0..alen{
                aligner.gmat_buffer[lid][ii].match_weight = x[ii].0;
                aligner.gmat_buffer[lid][ii].del_weight = x[ii].1;
                aligner.gmat_buffer[lid][ii].connected_weight = x[ii].2;
                aligner.gmat_buffer[lid][ii].gapped_weight = x[ii].3;
            }
            aligner.gmat_buffer[lid][alen].connected_weight = x[alen].2;
            aligner.gmat_buffer[lid][alen].gapped_weight = x[alen].3;
        }
        
        aligner.gmat_set(lid,alen,&(vec![0.0;aligner.vec_size]),1.0,0.0
        ,1.0,0.0);
        
        let mut idd = -1_i32;
        let mut pidd:Vec<i32> = vec![];
        if register_id{
            idd = aid as i32;
            pidd.push(aid as i32);
        }
        return ScoredSequence{
            id:idd,
            gmatbuff_id:lid as i32,
            alibuff_idx:vec![aid],
            num_sequences:1,
            primary_ids:pidd,
            alignment_length:alen,
            weight_sum:1.0
        }
    }

    pub fn get_alignment_length(&self)->usize{
        return self.alignment_length;
    }

    pub fn check_buff(&mut self,aligner:&mut ScoredSeqAligner){
        self.set_alignment_length(self.alignment_length, aligner);
    }

    pub fn set_alignment_length(&mut self,len:usize,aligner:&mut ScoredSeqAligner){
        self.alignment_length = len;
        if self.gmatbuff_id < 0{
            let lid = aligner.get_unused_gmatbuffid();
            aligner.register_gmatbuff(lid,self.alignment_length);
            self.gmatbuff_id = lid as i32;
        }
        while self.alibuff_idx.len() < self.primary_ids.len(){
            let aid = aligner.get_unused_alibuffid();
            aligner.register_alibuff(aid,self.alignment_length); 
            self.alibuff_idx.push(aid);
        }
    }

    pub fn get_num_seq(&self)-> usize{
        return self.num_sequences;
    }

    pub fn release_buffs(&mut self,aligner:&mut ScoredSeqAligner){
        aligner.release_gmatbuff(self.gmatbuff_id as usize);
        self.gmatbuff_id = -1;
        for pp in self.alibuff_idx.iter(){
            aligner.release_alibuff(*pp);
        }
        self.alibuff_idx.clear();
    }

    pub fn set_id(&mut self,idd:i32){
        if self.primary_ids.len() == 0{
            self.id = idd;
            self.primary_ids.push(idd); 
        }else{
            self.id = idd;
        }
    }

    pub fn get_index_of(& self,id:i32)->Option<usize>{
        let mut targetidx:Option<usize> = None;
        for (ee,ii) in self.primary_ids.iter().enumerate(){
            if *ii == id{
                targetidx = Some(ee);
                break;
            }
        }
        return targetidx;
    }


    /*
    pub fn split_one_with_id(&mut self,id:i32,aligner:&ScoredSeqAligner)->ScoredSequence{
        let mut targetidx:usize = self.get_index_of(id);
        return self.split_one(targetidx, aligner);
    }
    */
    
    pub fn calc_gmat(&mut self,aligner:&mut ScoredSeqAligner){
        
        let lid:usize = self.gmatbuff_id as usize;
        let alilen = self.alignment_length;
        let veclen = aligner.vec_size;
        let mut all_weights:f32 = 0.0;
        for ppid in self.primary_ids.iter(){
            all_weights += aligner.weights[*ppid as usize];
        }
        //let weights = sequence_weighting::calc_henikoff_henikoff_weight(&aligner.alignment_buffer,&self.alibuff_ids,alilen);
        let mut aacount:Vec<usize> = vec![0;alilen]; // 配列毎の次に参照するべき位置
        for alipos in 0..alilen{
            (aligner.gmat_buffer[lid][alipos].match_vec).fill(0.0);
            let mut sum_weight = 0.0;
            let mut sum_weight_del = 0.0;
            let mut ungapratio = 0.0;
            let mut gapratio = 0.0;
            for (eii,(aidx,ppid)) in self.alibuff_idx.iter().zip(self.primary_ids.iter()).enumerate(){
                if *ppid < 0{
                    panic!("???");
                }
                // 前後で繋がっていて GAPOPEN が必要なもののウェイトを取る
                if alipos == 0{
                    if aligner.ali_get(*aidx,alipos) != GAP_CHAR{
                        ungapratio += aligner.weights[*ppid as usize]*aligner.gmat_buffer[lid][aacount[eii]].connected_weight;
                        gapratio +=  aligner.weights[*ppid as usize]*aligner.gmat_buffer[lid][aacount[eii]].gapped_weight;
                    }else{
                        gapratio +=  aligner.weights[*ppid as usize]*(aligner.gmat_buffer[lid][aacount[eii]].gapped_weight
                            +aligner.gmat_buffer[lid][aacount[eii]].connected_weight);
                    }
                }else{
                    if aligner.alignment_buffer[*aidx][alipos-1] != GAP_CHAR && aligner.alignment_buffer[*aidx][alipos] != GAP_CHAR{
                        ungapratio += aligner.weights[*ppid as usize]*aligner.gmat_buffer[lid][aacount[eii]].connected_weight;
                        gapratio +=  aligner.weights[*ppid as usize]*aligner.gmat_buffer[lid][aacount[eii]].gapped_weight;
                    }else{
                        gapratio +=  aligner.weights[*ppid as usize]*(aligner.gmat_buffer[lid][aacount[eii]].gapped_weight
                            +aligner.gmat_buffer[lid][aacount[eii]].connected_weight);
                    }
                }

                if aligner.ali_get(*aidx,alipos) != GAP_CHAR{
                    for vv in 0..veclen{
                        aligner.gmat_buffer[lid][alipos].match_vec[vv] += aligner.gmat_buffer[*ppid as usize][aacount[eii]].match_vec[vv]
                        *aligner.weights[*ppid as usize]*aligner.gmat_buffer[*ppid as usize][aacount[eii]].match_weight;
                    }
                    aacount[eii] += 1;
                    sum_weight += aligner.weights[*ppid as usize]*aligner.gmat_buffer[*ppid as usize][aacount[eii]].match_weight;
                    sum_weight_del += aligner.weights[*ppid as usize]*aligner.gmat_buffer[*ppid as usize][aacount[eii]].del_weight;
                }else{
                    sum_weight_del += aligner.weights[*ppid as usize];
                }
            }
            aligner.gmat_buffer[lid][alipos].connected_weight = ungapratio/all_weights;
            aligner.gmat_buffer[lid][alipos].gapped_weight = gapratio/all_weights;

            if sum_weight > 0.0{
                element_multiply(&mut aligner.gmat_buffer[lid][alipos].match_vec,1.0/sum_weight);
            }
            let pallweight = sum_weight+sum_weight_del;
            aligner.gmat_buffer[lid][alipos].match_weight = sum_weight/pallweight;
            aligner.gmat_buffer[lid][alipos].del_weight = sum_weight_del/pallweight;
        }

        // 最後の残基以降のギャップ
        let mut ungapratio = 0.0;
        let mut gapratio = 0.0;
        for (eii,(aidx,ppid)) in self.alibuff_idx.iter().zip(self.primary_ids.iter()).enumerate(){
            if aligner.alignment_buffer[*aidx][alilen-1] != GAP_CHAR{
                ungapratio += aligner.weights[*ppid as usize]*aligner.gmat_buffer[lid][aacount[eii]].connected_weight;
                gapratio +=  aligner.weights[*ppid as usize]*aligner.gmat_buffer[lid][aacount[eii]].gapped_weight;
            }else{
                gapratio +=  aligner.weights[*ppid as usize]*(aligner.gmat_buffer[lid][aacount[eii]].gapped_weight
                    +aligner.gmat_buffer[lid][aacount[eii]].connected_weight);
            }
        }
        aligner.gmat_buffer[lid][alilen].connected_weight = ungapratio/all_weights;
        aligner.gmat_buffer[lid][alilen].gapped_weight = gapratio/all_weights;
        self.num_sequences = self.primary_ids.len();
    }
    
}
