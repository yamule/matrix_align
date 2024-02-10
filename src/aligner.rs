use std::collections::HashMap;
use std::collections::HashSet;

use super::*;

#[allow(unused_imports)]
use std::time::Instant;


const NUM_CHARTYPE:usize = 28;
const ACCEPTABLE_CHARS:&str = "ACDEFGHIKLMNPQRSTVWXY-";
const GAP_CHAR:char = '-';
const DIREC_UPLEFT:u8 = 0;
const DIREC_UP:u8 = 1;
const DIREC_LEFT:u8 = 2;



pub struct ScoredSeqAligner{
    pub dp_matrix:Vec<Vec<Vec<f32>>>,
    pub path_matrix:Vec<Vec<Vec<u8>>>,
    pub charmap:Vec<usize>,
    pub vec_size:usize,
    pub alen:usize,
    pub blen:usize,
    pub weights:Vec<f32>,
    pub alignment_buffer:Vec<Vec<char>>,//全配列がアラインメントされた状態で保持されている
    pub pssm_buffer:Vec<Vec<(Vec<f32>,f32,f32)>>,//その内積がスコア計算に使われるベクトル。第二要素は前の要素と結合している配列のウェイト、最後の残基についてはベクトルの要素は使わない
    pub alibuff_used:Vec<bool>,
    pub pssmbuff_used:Vec<bool>,
    pub alibuff_next:usize,
    pub pssmbuff_next:usize,
    pub penalty_warning:bool
}
impl ScoredSeqAligner {
    pub fn new(vec_size:usize,buff_len:usize,buff_seqnum:usize)->ScoredSeqAligner{
        let dp_matrix:Vec<Vec<Vec<f32>>> = vec![vec![vec![];1];1];
        let path_matrix:Vec<Vec<Vec<u8>>> = vec![vec![vec![];1];1];
        let mut charmap:Vec<usize> = vec![NUM_CHARTYPE;256];
        let alignment_buffer = vec![vec!['-';buff_len];buff_seqnum];
        let pssm_buffer:Vec<Vec<(Vec<f32>,f32,f32)>> = vec![vec![(vec![0.0;vec_size],0.0,0.0);buff_len+1];buff_seqnum];
        let buff_used:Vec<bool> = vec![false;buff_seqnum];
        let pssmbuff_used:Vec<bool> = vec![false;buff_seqnum];
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
            ,pssm_buffer:pssm_buffer
            ,pssmbuff_used:pssmbuff_used
            ,alibuff_next:0
            ,pssmbuff_next:0
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

    pub fn get_unused_pssmbuffid(&mut self)->usize{
        for ii in self.pssmbuff_next..self.pssmbuff_used.len(){
            if !self.pssmbuff_used[ii]{
                self.pssmbuff_next = ii+1;
                return ii;
            }
        }
        for ii in 0..self.pssmbuff_next{
            if !self.pssmbuff_used[ii]{
                self.pssmbuff_next = ii+1;
                return ii;
            }
        }
        let newid = self.pssmbuff_used.len();
        self.pssmbuff_used.push(false);
        self.pssm_buffer.push(vec![(vec![0.0;NUM_CHARTYPE],0.0,0.0);self.pssm_buffer[0].len()]);
        return newid;
    }
    
    pub fn register_pssmbuff(&mut self,id:usize,len:usize){
        assert!(!self.pssmbuff_used[id]);
        self.format_pssmbuff(id, len);
        self.check_pssmbuff_length(id, len);
        self.pssmbuff_used[id] = true;
    }

    pub fn copy_pssm_a_to_b(&mut self,a:usize,b:usize,len:usize){
        assert!(self.pssm_buffer[a].len() >= len);
        self.check_pssmbuff_length(b,len);
        for i in 0..len{
            self.pssm_buffer[b][i].0 = self.pssm_buffer[a][i].0.clone();
            self.pssm_buffer[b][i].1 = self.pssm_buffer[a][i].1;
        }
        self.pssm_buffer[b][len].1 = self.pssm_buffer[a][len].1;
    }
    
    pub fn copy_ali_a_to_b(&mut self,a:usize,b:usize,len:usize){
        assert!(self.alignment_buffer[a].len() >= len);
        self.check_alibufflen(b,len);

        self.alignment_buffer[b] = self.alignment_buffer[a].clone();

    }

    pub fn release_pssmbuff(&mut self,id:usize){
        assert!(self.pssmbuff_used[id]);
        self.pssmbuff_used[id] = false;
    }

    //文字列を保持しておくバッファのサイズを確認し、小さい場合は冗長性を少し加えて大きくする
    pub fn check_pssmbuff_length(&mut self,id:usize,len:usize){
        if self.pssm_buffer[id].len() < len+1{
            self.pssm_buffer[id] = vec![(vec![0.0;self.vec_size],0.0,0.0);len+6];
        }
    }
    
    pub fn format_pssmbuff(&mut self,id:usize,len:usize){
        self.check_pssmbuff_length(id, len);
        for vv in 0..=len{
            self.pssm_buffer[id][vv].0.fill(0.0);
            self.pssm_buffer[id][vv].1 = 0.0;
        }
    }

    //呼び出し元から直接 dot_product に飛ばしてもいいかも
    pub fn calc_match_score(a:&Vec<f32>,b:&Vec<f32>)->f32{
        unsafe{
            return matrix_process::dot_product(a,b);
        }
    }
    
    pub fn pssm_colget(&self,lid:usize,pos:usize)->&(Vec<f32>,f32,f32){
        return &self.pssm_buffer[lid][pos];
    }

    //別関数になるかも
    pub fn pssm_set(&mut self,lid:usize,pos:usize,val:&(Vec<f32>,f32,f32)){
        self.pssm_buffer[lid][pos].0.copy_from_slice(&val.0);
        self.pssm_buffer[lid][pos].1 = val.1;//二番目は前の残基とギャップにならず連続しているもののウェイト合計。Gap Open として考慮される
        self.pssm_buffer[lid][pos].2 = val.2;//三番目は前の残基とギャップになっているもののウェイト合計。Gap Extend として考慮される
    }

    //この辺違う関数になるはず
    /*
    pub fn pssm_add(&mut self,lid:usize,pos:usize,c:char,num:i32){
        self.pssm_buffer[lid][pos][self.charmap[c as usize]] += num;
        assert!(self.pssm_buffer[lid][pos][self.charmap[c as usize]] >= 0);
    }

    pub fn pssm_add_i(&mut self,lid:usize,pos:usize,c:usize,num:i32){
        self.pssm_buffer[lid][pos][c] += num;
        assert!(self.pssm_buffer[lid][pos][c] >= 0);
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
        
        let alid = a.pssmbuff_id as usize;
        let blid = b.pssmbuff_id as usize;
        
        let mut currentpenal = 0.0;
        
        // B 側 N 末にギャップを入れる
        for ii in 0..=aalen{
            if ii != 0{
                if ii == 1{
                    currentpenal = self.pssm_colget(alid,ii-1).1*gap_open_penalty+self.pssm_colget(alid,ii-1).2*gap_extension_penalty;
                }else{
                    currentpenal = self.pssm_colget(alid,ii-1).1*gap_extension_penalty+self.pssm_colget(alid,ii-1).2*gap_extension_penalty;
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
                    currentpenal = self.pssm_colget(blid,ii-1).1*gap_open_penalty+self.pssm_colget(blid,ii-1).2*gap_extension_penalty;
                }else{
                    currentpenal = self.pssm_colget(blid,ii-1).1*gap_extension_penalty+self.pssm_colget(blid,ii-1).2*gap_extension_penalty;
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

        for ii in 1..=aalen{
            for jj in 1..=bblen{
                let acol = self.pssm_colget(alid,ii-1);
                let bcol = self.pssm_colget(blid,jj-1);
                let sc:f32 = ScoredSeqAligner::calc_match_score(&acol.0,&bcol.0);
                

                let diag_m:f32 = self.dp_matrix[ii-1][jj-1][DIREC_UPLEFT as usize] + sc;
                let diag_l:f32 = self.dp_matrix[ii-1][jj-1][DIREC_LEFT as usize] + sc;
                let diag_u:f32 = self.dp_matrix[ii-1][jj-1][DIREC_UP as usize] + sc;


                let lef_m:f32 = self.dp_matrix[ii-1][jj][DIREC_UPLEFT as usize] + acol.1*gap_open_penalty + acol.2*gap_extension_penalty;
                let lef_l:f32 = self.dp_matrix[ii-1][jj][DIREC_LEFT as usize] + acol.1*gap_extension_penalty + acol.2*gap_extension_penalty;
                let lef_u:f32 = self.dp_matrix[ii-1][jj][DIREC_UP as usize] + acol.1*gap_open_penalty + acol.2*gap_extension_penalty;


                let up_m:f32 = self.dp_matrix[ii][jj-1][DIREC_UPLEFT as usize] + bcol.1*gap_open_penalty + bcol.2*gap_extension_penalty;
                let up_l:f32 = self.dp_matrix[ii][jj-1][DIREC_LEFT as usize] + bcol.1*gap_open_penalty + bcol.2*gap_extension_penalty;
                let up_u:f32 = self.dp_matrix[ii][jj-1][DIREC_UP as usize] + bcol.1*gap_extension_penalty + bcol.2*gap_extension_penalty;
                

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
        if a.id >= b.id{
            panic!("a.id must be smaller than b.id! {} {} ",a.id,b.id);
        }

        let anumaliseq = a.primary_ids.len();//alignment にだけ使う
        let bnumaliseq = b.primary_ids.len();
        let numallseq = anumaliseq+bnumaliseq;
        let alignment_length = alignment.len();

        let mut new_aids = vec![];
        let lid = self.get_unused_pssmbuffid();
        self.register_pssmbuff(lid,alignment_length);
        let alignment_length = alignment.len();
        
        self.format_pssmbuff(lid,alignment_length);
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


        return ScoredSequence{
            id:-1,
            num_sequences:a.num_sequences+b.num_sequences,
            primary_ids:primary_ids,
            pssmbuff_id:lid as i32,
            alibuff_idx:new_aids,
            alignment_length:alignment_length,
            weight_sum:sumweight
        };
    }
    pub fn make_msa(&mut self,mut sequences: Vec<ScoredSequence>,gap_open_penalty:f32,gap_extension_penalty:f32,profile_only:bool)
    -> (Vec<ScoredSequence>,f32){
        let mut next_id = sequences.iter().fold(0,|s,a| s.max(a.id)) as i32 + 10000;
        let mut final_score:f32 = 0.0;
        //必ず若い方の id を先に入れること
        let mut calcd:HashMap<(i32,i32),(Vec<(i32,i32)>,f32)> = HashMap::new();
        let mut calc_priority = false;
        while sequences.len() > 1{
            let numseq = sequences.len();
            let mut maxpair:(i32,i32) = (-1,-1);
            let mut maxscore:f32 = -10000.0;
            let mut prioritymap:HashMap<usize,f32> = HashMap::new();
            'outer:for ii in 0..numseq{
                let aid = sequences[ii].id;
                for jj in (ii+1)..numseq{
                    let bid = sequences[jj].id;
                    if aid < bid{
                        if !calcd.contains_key(&(aid,bid)){
                            calcd.insert((aid,bid),self.perform_dp(&sequences[ii],&sequences[jj],gap_open_penalty,gap_extension_penalty));
                        }
                        let dscore = calcd.get(&(aid,bid)).unwrap().1;
                        if maxpair.0 == -1 || dscore > maxscore{
                            maxpair = (aid,bid);
                            maxscore = dscore;
                        }
                        if !calc_priority{
                            break 'outer;
                        }else{
                            prioritymap.insert(jj,dscore);
                        }
                    }else{
                        if !calcd.contains_key(&(bid,aid)){
                            calcd.insert((bid,aid),self.perform_dp(&sequences[jj],&sequences[ii], gap_open_penalty,gap_extension_penalty));
                        }
                        let dscore = calcd.get(&(bid,aid)).unwrap().1;
                        if maxpair.0 == -1 || dscore > maxscore{
                            maxpair = (bid,aid);
                            maxscore = dscore;
                        }
                        
                        if !calc_priority{
                            break 'outer;
                        }else{
                            prioritymap.insert(jj,dscore);
                        }
                    }
                }
                break;
            }
            assert!(maxpair.0 > -1);
            if calc_priority{
                
                let mut newseq:Vec<(ScoredSequence,f32)> = vec![];
                for pp in sequences.into_iter().enumerate(){
                    if prioritymap.contains_key(&pp.0){
                        newseq.push((pp.1,*prioritymap.get(&pp.0).unwrap()));
                    }else{
                        newseq.push((pp.1,100.0));
                    }
                }
                newseq.sort_by(|a,b|b.1.partial_cmp(&a.1).unwrap());
                sequences = newseq.into_iter().map(|a|a.0).collect();
                
                calc_priority = false;
            }
            let mut aseq:Option<ScoredSequence> = None;
            let mut bseq:Option<ScoredSequence> = None;
            for ii in 0..numseq{
                if sequences[ii].id == maxpair.0{
                    aseq = Some(sequences.remove(ii));
                    if let Some(_) = bseq{
                        break;
                    }
                    if sequences[ii].id == maxpair.1{
                        bseq = Some(sequences.remove(ii));
                        break;
                    }
                }
                
                if sequences[ii].id == maxpair.1{
                    bseq = Some(sequences.remove(ii));
                    if let Some(_) = aseq{
                        break;
                    }
                    if sequences[ii].id == maxpair.0{
                        aseq = Some(sequences.remove(ii));
                        break;
                    }
                }
            }

            let dpres = calcd.remove(&maxpair).unwrap();
            final_score = maxscore;
            let mut newgroup = ScoredSeqAligner::make_alignment(self,aseq.unwrap(),bseq.unwrap(),dpres.0,profile_only);
            newgroup.set_id(next_id);
            next_id += 1;
            sequences.insert(0,newgroup);
            if calcd.len() > sequences.len()*sequences.len(){
                let hs:HashSet<i32> = sequences.iter().map(|a| a.id).collect();
                let mut rempairs:Vec<(i32,i32)> = vec![];
                for hh in calcd.iter(){
                    if !hs.contains(&(hh.0).0) || !hs.contains(&(hh.0).1){
                        rempairs.push((hh.0).clone());
                    } 
                }
                for rr in rempairs.iter(){
                    calcd.remove(rr);
                }
            }
        }
        return (sequences,final_score);
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
    pub pssmbuff_id:i32,
    pub alibuff_idx:Vec<usize>,
    pub alignment_length:usize,
    pub weight_sum:f32
}

impl ScoredSequence{
    pub fn new(a:Vec<char>,a_pssm:Vec<Vec<f32>>,aligner:&mut ScoredSeqAligner,register_id:bool)-> ScoredSequence{
        let lid = aligner.get_unused_pssmbuffid();
        let aid = aligner.get_unused_alibuffid();
        aligner.register_pssmbuff(lid,a.len());
        aligner.register_alibuff(aid,a.len());
        let alen = a.len();
        for aa in a.into_iter().zip(a_pssm.into_iter()).enumerate(){
            if aligner.charmap[(aa.1).0 as usize] != NUM_CHARTYPE{
                aligner.pssm_set(lid,aa.0,&((aa.1).1,1.0,0.0)); 
                aligner.ali_set(aid,aa.0,(aa.1).0);
            }else{
                panic!("can not process {}!",(aa.1).0);
            }
        }
        
        aligner.pssm_set(lid,alen,&(vec![0.0;aligner.vec_size],1.0,0.0));
        
        let mut idd = -1_i32;
        let mut pidd:Vec<i32> = vec![];
        if register_id{
            idd = aid as i32;
            pidd.push(aid as i32);
        }
        return ScoredSequence{
            id:idd,
            pssmbuff_id:lid as i32,
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


    /*また今度
    pub fn make_consensus_sequence(&self,aligner:&mut ScoredSeqAligner)->ScoredSequence{
        let mut counter:Vec<i32> = vec![0;NUM_CHARTYPE];　
        let mut _maxindex:i32 = -1;
        let mut maxcount:i32 = 0;
        let mut count_all:i32 = 0;
        for a_ in 0..self.alignment_length{   
            for (ee,vv) in aligner.pssm_colget(self.pssmbuff_id,_a).iter().enumerate(){
                if aligner.charmap['-' as usize] == ee{
                    continue;
                }
                counter[ee] += *vv;
                count_all +=  *vv;
                if maxcount < counter[ee]{
                    maxcount = counter[ee];
                    _maxindex = ee as i32;
                }
            }
            if _maxindex == -1{
                //外した方が良いが・・・
            }else{

            }
        }
    }
    */


    pub fn make_copy(&self,aligner:&mut ScoredSeqAligner)->ScoredSequence{
        let id = self.id;
        let alignment_length = self.get_alignment_length();
        let primary_ids = self.primary_ids.clone();
        let lid = aligner.get_unused_pssmbuffid();
        aligner.register_pssmbuff(lid, alignment_length);
        aligner.copy_pssm_a_to_b(self.pssmbuff_id as usize,lid,alignment_length);

        let nseq = self.alibuff_idx.len();
        let mut alibuff_ids:Vec<usize> = vec![];
        for nn in 0..nseq{ 
            let aid = aligner.get_unused_alibuffid();
            alibuff_ids.push(aid);
            aligner.register_alibuff(aid, alignment_length);
            aligner.copy_ali_a_to_b(self.alibuff_idx[nn],aid,alignment_length);
        }
        return ScoredSequence{
            id:id,
            primary_ids:primary_ids,
            pssmbuff_id:lid as i32,
            alibuff_idx:alibuff_ids,
            num_sequences:self.num_sequences,
            alignment_length:alignment_length,
            weight_sum:self.weight_sum
        };
    }

    pub fn check_buff(&mut self,aligner:&mut ScoredSeqAligner){
        self.set_alignment_length(self.alignment_length, aligner);
    }

    pub fn set_alignment_length(&mut self,len:usize,aligner:&mut ScoredSeqAligner){
        self.alignment_length = len;
        if self.pssmbuff_id < 0{
            let lid = aligner.get_unused_pssmbuffid();
            aligner.register_pssmbuff(lid,self.alignment_length);
            self.pssmbuff_id = lid as i32;
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
        aligner.release_pssmbuff(self.pssmbuff_id as usize);
        self.pssmbuff_id = -1;
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
    
    pub fn calc_pssm(&mut self,aligner:&mut ScoredSeqAligner){
        let lid:usize = self.pssmbuff_id as usize;
        let alilen = self.alignment_length;

        //let weights = sequence_weighting::calc_henikoff_henikoff_weight(&aligner.alignment_buffer,&self.alibuff_ids,alilen);
        for alipos in 0..alilen{
            (aligner.pssm_buffer[lid][alipos].0).fill(0.0);
            for (_eii,(aidx,ppid)) in self.alibuff_idx.iter().zip(self.primary_ids.iter()).enumerate(){
                if *ppid < 0{
                    panic!("???");
                }
                aligner.pssm_buffer[lid][alipos].0[aligner.charmap[aligner.alignment_buffer[*aidx][alipos] as usize]]
                 += aligner.weights[*ppid as usize];
            }
            let mut ungapratio = 0.0;
            let mut gapratio = 0.0;
            for (_eii,(aidx,ppid)) in self.alibuff_idx.iter().zip(self.primary_ids.iter()).enumerate(){
                // 前後で繋がっていて GAPOPEN が必要なもののウェイトを取る
                if alipos == 0{
                    if aligner.alignment_buffer[*aidx][alipos] != GAP_CHAR{
                        ungapratio +=  aligner.weights[*ppid as usize];
                    }else{
                        gapratio +=  aligner.weights[*ppid as usize];
                    }
                }else{
                    if aligner.alignment_buffer[*aidx][alipos-1] != GAP_CHAR && aligner.alignment_buffer[*aidx][alipos] != GAP_CHAR{
                        ungapratio +=  aligner.weights[*ppid as usize];
                    }else{
                        gapratio +=  aligner.weights[*ppid as usize];
                    }
                }
            }
            aligner.pssm_buffer[lid][alipos].1 = ungapratio;
            aligner.pssm_buffer[lid][alipos].2 = gapratio;
        }

        // 最後の残基以降のギャップ
        let mut ungapratio = 0.0;
        let mut gapratio = 0.0;
        for (_eii,(aidx,ppid)) in self.alibuff_idx.iter().zip(self.primary_ids.iter()).enumerate(){
            if aligner.alignment_buffer[*aidx][alilen-1] != GAP_CHAR{
                ungapratio +=  aligner.weights[*ppid as usize];
            }else{
                gapratio +=  aligner.weights[*ppid as usize];
            }
        }
        aligner.pssm_buffer[lid][alilen].1 = ungapratio;
        aligner.pssm_buffer[lid][alilen].2 = gapratio;
        self.num_sequences = self.primary_ids.len();
    }
    
}
