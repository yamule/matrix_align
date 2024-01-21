use std::collections::HashMap;
use std::collections::HashSet;
use std::io;
use regex::Regex;
//use rand::prelude::*;

use std::time::Instant;


const NUM_CHARTYPE:usize = 5;

const DIREC_UPLEFT:u8 = 0;
const DIREC_UP:u8 = 1;
const DIREC_LEFT:u8 = 2;
pub struct FreqSeqAligner{
    pub dp_matrix:Vec<Vec<f32>>,
    pub path_matrix:Vec<Vec<u8>>,
    pub charmap:Vec<usize>,
    pub alen:usize,
    pub blen:usize,
    pub alignment_buffer:Vec<Vec<char>>,
    pub letters_buffer:Vec<Vec<Vec<i32>>>,
    pub alibuff_used:Vec<bool>,
    pub letbuff_used:Vec<bool>
}
impl FreqSeqAligner {
    pub fn new(buff_len:usize,buff_seqnum:usize)->FreqSeqAligner{
        let dp_matrix:Vec<Vec<f32>> = vec![vec![0.0;buff_len+1];buff_len+1];
        let path_matrix:Vec<Vec<u8>> = vec![vec![0;buff_len+1];buff_len+1];
        let mut charmap:Vec<usize> = vec![NUM_CHARTYPE;256];
        let alignment_buffer = vec![vec!['-';buff_len];buff_seqnum];
        let letters_buffer:Vec<Vec<Vec<i32>>> = vec![vec![vec![0;NUM_CHARTYPE];buff_len];buff_seqnum];
        let buff_used:Vec<bool> = vec![false;buff_seqnum];
        let letbuff_used:Vec<bool> = vec![false;buff_seqnum];
        charmap['A' as usize] = 0;
        charmap['T' as usize] = 1;
        charmap['G' as usize] = 2;
        charmap['C' as usize] = 3;
        charmap['-' as usize] = 4;
        return FreqSeqAligner{
            dp_matrix:dp_matrix
            ,path_matrix:path_matrix
            ,charmap
            ,alen:0
            ,blen:0
            ,alignment_buffer:alignment_buffer
            ,alibuff_used:buff_used
            ,letters_buffer:letters_buffer
            ,letbuff_used:letbuff_used
        };
    }
    
    pub fn get_unused_alibuffid(&mut self)->usize{
        for vv in self.alibuff_used.iter().enumerate(){
            if !*vv.1{
                return vv.0;
            }
        }
        let newid = self.alibuff_used.len();
        self.alibuff_used.push(false);
        self.alignment_buffer.push(vec!['-';self.alignment_buffer[0].len()]);
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

    pub fn get_unused_letbuffid(&mut self)->usize{
        for vv in self.letbuff_used.iter().enumerate(){
            if !*vv.1{
                return vv.0;
            }
        }
        let newid = self.letbuff_used.len();
        self.letbuff_used.push(false);
        self.letters_buffer.push(vec![vec![0;NUM_CHARTYPE];self.letters_buffer[0].len()]);
        return newid;
    }
    
    pub fn register_letbuff(&mut self,id:usize,len:usize){
        assert!(!self.letbuff_used[id]);
        self.format_letbuff(id, len);
        self.check_letbufflen(id, len);
        self.letbuff_used[id] = true;
    }

    pub fn copy_let_a_to_b(&mut self,a:usize,b:usize,len:usize){
        assert!(self.letters_buffer[a].len() >= len);
        self.check_letbufflen(b,len);
        for i in 0..len{
            for j in 0..NUM_CHARTYPE{
                self.letters_buffer[b][i][j] = self.letters_buffer[a][i][j];
            }
        }
    }
    
    pub fn copy_ali_a_to_b(&mut self,a:usize,b:usize,len:usize){
        assert!(self.alignment_buffer[a].len() >= len);
        self.check_alibufflen(b,len);
        for i in 0..len{
            self.alignment_buffer[b][i] = self.alignment_buffer[a][i];
        }
    }

    pub fn release_letbuff(&mut self,id:usize){
        assert!(self.letbuff_used[id]);
        self.letbuff_used[id] = false;
    }

    pub fn check_letbufflen(&mut self,id:usize,len:usize){
        if self.letters_buffer[id].len() < len{
            self.letters_buffer[id] = vec![vec![0;NUM_CHARTYPE];len+5];
        }
    }
    
    pub fn format_letbuff(&mut self,id:usize,len:usize){
        self.check_letbufflen(id, len);
        for vv in 0..len{
            for vvv in self.letters_buffer[id][vv].iter_mut(){
                *vvv = 0;
            }
        }
    }

    //align された場合のスコア（マイナス）
    pub fn calc_match_score(a:&Vec<i32>,b:&Vec<i32>)->f32{
        let mut counter:Vec<i32> = vec![0;NUM_CHARTYPE];
        let mut _maxindex:usize = 0;
        let mut maxcount:i32 = 0;
        let mut count_all:i32 = 0;
        for (ee,vv) in a.iter().enumerate(){
            counter[ee] += *vv;
            count_all +=  *vv;
            if maxcount < counter[ee]{
                maxcount = counter[ee];
                _maxindex = ee;
            }
        }
        
        for (ee,vv) in b.iter().enumerate(){
            counter[ee] += *vv;
            count_all +=  *vv;
            if maxcount < counter[ee]{
                maxcount = counter[ee];
                _maxindex = ee;
            }
        }
        return (maxcount as i32 - count_all as i32 ) as f32;
    }


    pub fn calc_score(&self,a_:&Vec<Vec<i32>>,start:usize,end:usize,additional_gap:usize) -> f32{
        let mut ret:f32 = 0.0;
        let mut counter:Vec<i32> = vec![0;NUM_CHARTYPE];
        for poss in start..=end{
            let a = &a_[poss];
            for ii in 0..NUM_CHARTYPE{
                counter[ii] = 0;
            }
            counter[self.charmap['-' as usize ]] = additional_gap as i32;
            let mut _maxindex:usize = self.charmap['-' as usize ];
            let mut maxcount:i32 = additional_gap as i32;
            let mut count_all:i32 = 0;
            for (ee,vv) in a.iter().enumerate(){
                counter[ee] += *vv;
                count_all +=  *vv;
                if maxcount < counter[ee]{
                    maxcount = counter[ee];
                    _maxindex = ee;
                }
            }
            ret += (maxcount as i32 - count_all as i32)  as f32;
        }
        return ret;
    }
    
    pub fn calc_match_ex_score(&self,a:&Vec<i32>,b:&Vec<i32>)->f32{
        let mut ret:f32 = 0.0;
        for (ii,(ee,vv)) in a.iter().zip(b.iter()).enumerate(){
            if ii == self.charmap['-' as usize]{
                continue;
            }
            ret += (ee*vv) as f32;
        }
        return ret;
    }

    pub fn calc_gap_score(&self,a:&Vec<i32>,anumseq:usize,bnumseq:usize)->f32{//bseq にギャップが入った場合のスコア（マイナス）
        let numgap = a[self.charmap['-' as usize]]+bnumseq as i32;
        let numseqq = anumseq+bnumseq;
        let mut maxcount = numgap;
        for ii in 0..NUM_CHARTYPE{
            if ii == self.charmap['-' as usize]{
                continue;
            }
            maxcount = maxcount.max(a[ii]);
        }
        return (maxcount as i32 - numseqq as i32) as f32;
    }
    pub fn reconstruct_matrix(&mut self,amax:usize,bmax:usize){
        self.dp_matrix = vec![vec![0.0;bmax+1];amax+1];
        self.path_matrix  = vec![vec![0;bmax+1];amax+1];
    }

    pub fn letters_colget(&self,lid:usize,pos:usize)->&Vec<i32>{
        return &self.letters_buffer[lid][pos];
    }
    
    pub fn letters_set(&mut self,lid:usize,pos:usize,c:char,num:i32){
        self.letters_buffer[lid][pos][self.charmap[c as usize]] = num as i32;
    }
    pub fn letters_get_i(&self,lid:usize,pos:usize,c:usize)->i32{
        return self.letters_buffer[lid][pos][c];
    }
    pub fn letters_set_i(&mut self,lid:usize,pos:usize,c:usize,num:i32){
        self.letters_buffer[lid][pos][c] = num as i32;
    }
    
    pub fn letters_get(&self,lid:usize,pos:usize,c:char)->i32{
        return self.letters_buffer[lid][pos][self.charmap[c as usize]];
    }
    
    pub fn letters_add(&mut self,lid:usize,pos:usize,c:char,num:i32){
        self.letters_buffer[lid][pos][self.charmap[c as usize]] += num;
        assert!(self.letters_buffer[lid][pos][self.charmap[c as usize]] >= 0);
    }

    pub fn letters_add_i(&mut self,lid:usize,pos:usize,c:usize,num:i32){
        self.letters_buffer[lid][pos][c] += num;
        assert!(self.letters_buffer[lid][pos][c] >= 0);
    }

    pub fn ali_set(&mut self,aid:usize,pos:usize,c:char){
        self.alignment_buffer[aid][pos] = c;
    }

    pub fn ali_get(&self,aid:usize,pos:usize)->char{
        return self.alignment_buffer[aid][pos];
    }
    pub fn ali_get_i(&self,aid:usize,pos:usize)->usize{
        return self.charmap[self.alignment_buffer[aid][pos] as usize];
    }

    pub fn perform_dp(&mut self,a:&FreqSequence,b:&FreqSequence,gap_penalty:f32)->(Vec<(i32,i32)>,f32) {
        let aalen = a.alignment_length;
        let bblen = b.alignment_length;
        let recflag = self.dp_matrix.len() <= aalen || self.dp_matrix[0].len() <= bblen;
        if recflag{
            self.reconstruct_matrix(aalen+25, bblen+25);
        }
        let anumseq = a.get_num_seq();
        let bnumseq = b.get_num_seq();
        let exbias:f32 =  1.0/((anumseq*bnumseq) as f32);
        let alid = a.letbuff_id as usize;
        let blid = b.letbuff_id as usize;
        //let gap_penalty= -exbias/2.0;
        for ii in 0..=aalen{
            if ii != 0{
                self.dp_matrix[ii][0] = self.dp_matrix[ii-1][0]
                 + self.calc_gap_score(self.letters_colget(alid,ii-1)
                 ,anumseq,bnumseq) + gap_penalty*(ii as f32);
            }else{
                self.dp_matrix[ii][0] = 0.0;
            }
            self.path_matrix[ii][0] = DIREC_LEFT;
        }
        for ii in 0..=bblen{
            if ii != 0{
                self.dp_matrix[0][ii] = self.dp_matrix[0][ii-1] 
                + self.calc_gap_score(self.letters_colget(blid,ii-1),bnumseq,anumseq) + gap_penalty*(ii as f32);
            }else{
                self.dp_matrix[0][ii] = 0.0;
            }
            self.path_matrix[0][ii] = DIREC_UP;
        }
        for ii in 1..=aalen{
            for jj in 1..=bblen{
                let acol = self.letters_colget(alid,ii-1);
                let bcol = self.letters_colget(blid,jj-1);
                let sc:f32 = FreqSeqAligner::calc_match_score(acol,bcol);
                let sc2:f32 = FreqSeqAligner::calc_match_ex_score(self,acol,bcol);
                let diag:f32 = self.dp_matrix[ii-1][jj-1] + sc + sc2*exbias;
                let leff:f32 = self.dp_matrix[ii-1][jj] + self.calc_gap_score(acol,anumseq,bnumseq) + gap_penalty;
                let upp:f32 = self.dp_matrix[ii][jj-1] + self.calc_gap_score(bcol,bnumseq,anumseq) + gap_penalty;
                //println!("{} {} {}",diag,leff,upp);
                if diag >= leff && diag >= upp{
                    self.dp_matrix[ii][jj] = diag;
                    self.path_matrix[ii][jj] = DIREC_UPLEFT;
                }else{
                    if leff > upp{
                        self.dp_matrix[ii][jj] = leff;
                        self.path_matrix[ii][jj] = DIREC_LEFT;
                    }else{
                        self.dp_matrix[ii][jj] = upp;
                        self.path_matrix[ii][jj] = DIREC_UP;
                    }
                }
            }
        }
        //panic!("{:?}",self.dp_matrix);
        let mut currentx = aalen;
        let mut currenty = bblen;
        let mut aligned_tuple:Vec<(i32,i32)> = vec![];
        while currentx > 0 || currenty > 0{
            if self.path_matrix[currentx][currenty] == DIREC_UPLEFT{
                currentx -= 1;
                currenty -= 1;
                aligned_tuple.push((currentx as i32,currenty as i32));
            }else if self.path_matrix[currentx][currenty] == DIREC_UP{
                currenty -= 1;
                aligned_tuple.push((-1,currenty as i32));
            }else if self.path_matrix[currentx][currenty] == DIREC_LEFT{
                currentx -= 1;
                aligned_tuple.push((currentx as i32,-1));
            }else{
                panic!("???");
            }
        }
        aligned_tuple.reverse();
        return (aligned_tuple,self.dp_matrix[aalen][bblen]);
    }
    pub fn make_alignment(&mut self
        ,mut a:FreqSequence
        ,mut b:FreqSequence
        ,alignment:Vec<(i32,i32)>
        ,profile_only:bool // true にすると alignment の文字は a に関してのみ保持する
    )->FreqSequence{
        if a.id >= b.id{
            panic!("a.id must be smaller than b.id! {} {} ",a.id,b.id);
        }

        let anumaliseq = a.primary_ids.len();//alignment にだけ使う
        let bnumaliseq = b.primary_ids.len();
        let numallseq = anumaliseq+bnumaliseq;
        let alignment_length = alignment.len();

        let mut new_aids = vec![];
        let lid = self.get_unused_letbuffid();
        self.register_letbuff(lid,alignment_length);
        let alignment_length = alignment.len();
        
        self.format_letbuff(lid,alignment_length);
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
                for ii in 0..NUM_CHARTYPE{
                    self.letters_add_i(lid, aa, ii, self.letters_get_i(a.letbuff_id as usize,poss,ii));
                }
                for ap in 0..anumaliseq{
                    self.ali_set(new_aids[ap],aa, self.ali_get(a.alibuff_ids[ap], poss as usize));
                }
            }else{
                self.letters_add(lid, aa, '-', anumaliseq as i32);
                for ap in 0..anumaliseq{
                    self.ali_set(new_aids[ap],aa,'-');
                }
            }
            if ppos.1 > -1{
                let poss = ppos.1 as usize;
                for ii in 0..NUM_CHARTYPE{
                    self.letters_add_i(lid, aa, ii
                        , self.letters_get_i(b.letbuff_id as usize,poss,ii));
                }
                if !profile_only{
                    for ap in 0..bnumaliseq{
                        self.ali_set(new_aids[ap+anumaliseq],aa,
                            self.ali_get(b.alibuff_ids[ap], poss as usize));
                    }
                }
            }else{
                self.letters_add(lid, aa, '-', b.get_num_seq() as i32);
                if !profile_only{
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
        a.release_buffs(self);
        b.release_buffs(self);

        return FreqSequence{
            id:-1,
            num_sequences:a.num_sequences+b.num_sequences,
            primary_ids:primary_ids,
            letbuff_id:lid as i32,
            alibuff_ids:new_aids,
            alignment_length:alignment_length
        };
    }
    pub fn make_msa(&mut self,mut sequences: Vec<FreqSequence>,gap_penalty:f32,profile_only:bool)
    -> (Vec<FreqSequence>,f32){
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
                            calcd.insert((aid,bid),self.perform_dp(&sequences[ii],&sequences[jj],gap_penalty));
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
                            calcd.insert((bid,aid),self.perform_dp(&sequences[jj],&sequences[ii],gap_penalty));
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
                
                let mut newseq:Vec<(FreqSequence,f32)> = vec![];
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
            let mut aseq:Option<FreqSequence> = None;
            let mut bseq:Option<FreqSequence> = None;
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
            let mut newgroup = FreqSeqAligner::make_alignment(self,aseq.unwrap(),bseq.unwrap(),dpres.0,profile_only);
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

    

    //primid_map は、FreqSequence 内の primary_ids → distmat の index へのマップ
    //複数 Seq は 全 vs 全 の平均
    pub fn make_msa_with_distmat(&mut self, mut sequences:Vec<FreqSequence>,distmat:Vec<Vec<f32>>,primid_map:HashMap<i32,usize>,gap_penalty:f32)
    ->(Vec<FreqSequence>,f32){
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
            let dpres =  self.perform_dp(&sequences[nearest_pair_index.0 as usize],&sequences[nearest_pair_index.1 as usize],gap_penalty);


            let mut aseq:Option<FreqSequence> = None;
            let mut bseq:Option<FreqSequence> = None;
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

            let mut newgroup = FreqSeqAligner::make_alignment(self,aseq.unwrap(),bseq.unwrap(),dpres.0,false);
            newgroup.set_id(next_id);
            next_id += 1;
            sequences.insert(0,newgroup);
        }
        return (sequences,-1.0);
    }

    
    pub fn refine_subregion(&mut self,seq:&FreqSequence,start:usize,end:usize
        ,rollstart:usize,rollend:usize,duration_start:&Instant){
        let mut prevscore = self.calc_score(&self.letters_buffer[seq.letbuff_id as usize], start, end , 0);
        let mut seq_frag:Vec<FreqSequence> =vec![];
        let mut num_allgapseq:usize = 0;
        let mut frag_pid_to_orig_uid:HashMap<i32,usize> = HashMap::new();
        let mut maxid = -100000;
        for (ii,_tt) in seq.alibuff_ids.iter().enumerate(){
            let tt = &self.alignment_buffer[*_tt];
            let mut nucc:Vec<char> = vec![];
            for ttt in tt[start..=end].iter(){
                if *ttt != '-'{
                    nucc.push(*ttt);
                }
            }
            if nucc.len() == 0{
                num_allgapseq += 1;
                continue;
            }
            let mut fseq = FreqSequence::new(nucc,self);
            fseq.set_id(ii as i32);
            maxid = maxid.max(fseq.id);
            frag_pid_to_orig_uid.insert(fseq.id,ii);
            seq_frag.push(fseq);
        }
        if seq_frag.len() == 1{
            for mut ss in seq_frag.into_iter(){
                ss.release_buffs(self);
            }
            return;
        }
        maxid += 100;
        for rp in 0..seq_frag.len(){
            if rp >= rollstart && rp <= rollend{
                let num_subsample = (seq_frag.len()/2).min(10);
                let mut dummy_ids:Vec<i32> = vec![];
                let mut sequences:Vec<FreqSequence> = vec![];
                for ss_ in 0..num_subsample{
                    let mut ss = seq_frag[ss_].make_copy(self);
                    let did = maxid+1;
                    maxid = did;
                    dummy_ids.push(did);
                    ss.set_id(did);
                    ss.set_first_primary_id(did);
                    sequences.push(ss);
                }

                //プロファイルのみ作成する
                let (mut sequences,pscore) = self.make_msa(sequences, 0.0,true);
                

                let mut copied:Vec<FreqSequence> = seq_frag.iter().map(|m|m.make_copy(self)).collect();
                copied.reverse();
                for ss in copied.into_iter(){
                    sequences.push(ss);
                }

                //println!("{:?}",&self.letters_buffer[seq_frag[0].letbuff_id as usize]);
                let prevlen = end-start+1;
                
                let (mut resseq_,pscore) = self.make_msa(sequences, 0.0,false);
                let mut resseq = resseq_.pop().unwrap();
                resseq.remove_seqs(&dummy_ids,self);
                resseq.recalc_freq(self);
                let currentscore = self.calc_score(&self.letters_buffer[resseq.letbuff_id as usize],0,resseq.alignment_length-1,num_allgapseq);
                //println!("{} {} ",prevscore,currentscore);
                /*
                for ps in resseq.alibuff_ids.iter().zip(resseq.primary_ids.iter()){
                    //println!("{:?}",ps);
                    let ss = &self.alignment_buffer[*ps.0];
                    for cc in 0..resseq.alignment_length{
                        print!("{}",ss[cc]);
                    }
                    println!("");
                }
                */
                if currentscore > prevscore && prevlen >= resseq.alignment_length{
                    prevscore = currentscore;
                    for aid in seq.alibuff_ids.iter(){
                        for ii in start..=end{
                            self.ali_set(*aid,ii,'-');
                        }
                    }
                    let lid = seq.letbuff_id as usize;
                    for ii in start..=end{
                        for jj in 0..NUM_CHARTYPE{
                            self.letters_set_i(lid, ii,jj,0);
                        }
                    }
                    for (rid,frag_aid) in resseq.primary_ids.iter().zip(resseq.alibuff_ids.iter()){
                        let orig_aid = seq.alibuff_ids[*frag_pid_to_orig_uid.get(rid).unwrap()];
                        for jj in 0..resseq.get_alignment_length(){
                            let c = self.ali_get(*frag_aid,jj);
                            self.ali_set(orig_aid,start+jj, 
                            c);
                            self.letters_add(lid, start+jj,c,1);
                        }
                    }
                }
                resseq.release_buffs(self);
            }
            let p = seq_frag.remove(0);
            seq_frag.push(p);
            
            //Genocon 用
            let duration = duration_start.elapsed();
            if duration.as_millis() > 9750{
                break;
            }
        }
        for mut ss in seq_frag.into_iter(){
            ss.release_buffs(self);
        }


        //途中
    }
    

    pub fn calc_alignment_dist(&self,seq:&FreqSequence)->(Vec<Vec<f32>>,HashMap<i32,usize>){
        let mut pids = seq.primary_ids.clone();
        pids.sort();
        let primaryid_map:HashMap<i32,usize> = pids.iter().enumerate().map(|m|(*m.1,m.0)).collect();
        assert_eq!(primaryid_map.len(),pids.len());
        
        let mut dist:Vec<Vec<f32>> = vec![vec![0.0;pids.len()];pids.len()];
        for ii in 0..seq.alibuff_ids.len(){
            for jj in (ii+1)..seq.alibuff_ids.len(){
                let mut diff:i32 = 0;
                for kk in 0..seq.alibuff_ids.len(){
                    if self.ali_get(seq.alibuff_ids[ii],kk)
                        != self.ali_get(seq.alibuff_ids[jj],kk){
                        diff += 1;
                    }
                }
                dist[ii][jj] = diff as f32;
            }
        }

        return (dist,primaryid_map);
    }
    
}




pub struct FreqSequence{
    pub id:i32,
    //pub letters:Vec<Vec<usize>>,//カラムごとの文字の出現数
    //pub alignments:Vec<Vec<char>>,
    pub num_sequences:usize,
    pub primary_ids:Vec<i32>,
    pub letbuff_id:i32,
    pub alibuff_ids:Vec<usize>,
    pub alignment_length:usize
}

impl FreqSequence{
    pub fn new(a:Vec<char>,aligner:&mut FreqSeqAligner)-> FreqSequence{
        let lid = aligner.get_unused_letbuffid();
        let aid = aligner.get_unused_alibuffid();
        aligner.register_letbuff(lid,a.len());
        aligner.register_alibuff(aid,a.len());
        for aa in a.iter().enumerate(){
            if aligner.charmap[*aa.1 as usize] != NUM_CHARTYPE{
                aligner.letters_set(lid,aa.0,*aa.1,1);
                aligner.ali_set(aid,aa.0,*aa.1);
            }else{
                panic!("can not process {}!",aa.1);
            }
        }
        return FreqSequence{
            id:-1,
            letbuff_id:lid as i32,
            alibuff_ids:vec![aid],
            num_sequences:1,
            primary_ids:vec![],//id はまだ発行しない。した方が良いか？
            alignment_length:a.len()
        }
    }

    pub fn get_alignment_length(&self)->usize{
        return self.alignment_length;
    }


    /*また今度
    pub fn make_consensus_sequence(&self,aligner:&mut FreqSeqAligner)->FreqSequence{
        let mut counter:Vec<i32> = vec![0;NUM_CHARTYPE];　
        let mut _maxindex:i32 = -1;
        let mut maxcount:i32 = 0;
        let mut count_all:i32 = 0;
        for a_ in 0..self.alignment_length{   
            for (ee,vv) in aligner.letters_colget(self.letbuff_id,_a).iter().enumerate(){
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


    pub fn make_copy(&self,aligner:&mut FreqSeqAligner)->FreqSequence{
        let id = self.id;
        let alignment_length = self.get_alignment_length();
        let primary_ids = self.primary_ids.clone();
        let lid = aligner.get_unused_letbuffid();
        aligner.register_letbuff(lid, alignment_length);
        aligner.copy_let_a_to_b(self.letbuff_id as usize,lid,alignment_length);

        let nseq = self.alibuff_ids.len();
        let mut alibuff_ids:Vec<usize> = vec![];
        for nn in 0..nseq{ 
            let aid = aligner.get_unused_alibuffid();
            alibuff_ids.push(aid);
            aligner.register_alibuff(aid, alignment_length);
            aligner.copy_ali_a_to_b(self.alibuff_ids[nn],aid,alignment_length);
        }
        return FreqSequence{
            id:id,
            primary_ids:primary_ids,
            letbuff_id:lid as i32,
            alibuff_ids:alibuff_ids,
            num_sequences:self.num_sequences,
            alignment_length:alignment_length
        };
    }

    pub fn check_buff(&mut self,aligner:&mut FreqSeqAligner){
        self.set_alignment_length(self.alignment_length, aligner);
    }

    pub fn set_alignment_length(&mut self,len:usize,aligner:&mut FreqSeqAligner){
        self.alignment_length = len;
        if self.letbuff_id < 0{
            let lid = aligner.get_unused_letbuffid();
            aligner.register_letbuff(lid,self.alignment_length);
            self.letbuff_id = lid as i32;
        }
        while self.alibuff_ids.len() < self.primary_ids.len(){
            let aid = aligner.get_unused_alibuffid();
            aligner.register_alibuff(aid,self.alignment_length); 
            self.alibuff_ids.push(aid);
        }
    }

    pub fn get_num_seq(&self)-> usize{
        return self.num_sequences;
    }

    pub fn release_buffs(&mut self,aligner:&mut FreqSeqAligner){
        aligner.release_letbuff(self.letbuff_id as usize);
        self.letbuff_id = -1;
        for pp in self.alibuff_ids.iter(){
            aligner.release_alibuff(*pp);
        }
        self.alibuff_ids.clear();
    }

    pub fn set_first_primary_id(&mut self,idd:i32){
        assert!(self.primary_ids.len() == 1);
        self.primary_ids[0] = idd;
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
    pub fn split_one_with_id(&mut self,id:i32,aligner:&FreqSeqAligner)->FreqSequence{
        let mut targetidx:usize = self.get_index_of(id);
        return self.split_one(targetidx, aligner);
    }
    */

    //remids に入った primary id を持つ配列を削除する
    pub fn remove_seqs(&mut self,remids:&Vec<i32>,aligner:&mut FreqSeqAligner){
        let lid:usize = self.letbuff_id as usize;
        for ii in remids.iter(){
            let idx_:Option<usize> = self.get_index_of(*ii);
            if let Some(idx) = idx_{
                let aid = self.alibuff_ids.remove(idx);
                let _idd = self.primary_ids.remove(idx);
                let lnum = self.alignment_length;
                for jj in 0..lnum{
                    aligner.letters_add(lid,jj,aligner.ali_get(aid,jj), -1);
                }
                aligner.release_alibuff(aid);
            }
        }
    }
    
    
    
    pub fn recalc_freq(&mut self,aligner:&mut FreqSeqAligner){
        let lid:usize = self.letbuff_id as usize;
        let aidlen = self.alibuff_ids.len();
        let alilen = self.alignment_length;
        for alipos in 0..alilen{
            let mut lett:Vec<usize> = vec![0;NUM_CHARTYPE];
            for ii in 0..aidlen{
                lett[aligner.ali_get_i(self.alibuff_ids[ii],alipos)] += 1;
            }
            for jj in 0..NUM_CHARTYPE{
                aligner.letters_set_i(lid,alipos,jj,lett[jj] as i32);
            }
        }
        self.num_sequences = self.primary_ids.len();
    }
    

    pub fn calc_penalty(&self,aligner:&FreqSeqAligner)->(Vec<f32>,usize){
        let mut ret:Vec<f32> = vec![0.0;self.primary_ids.len()];
        let lnum = self.alignment_length;
        let mut topindex:usize = 0;
        let mut maxpenal:f32 = 0.0;
        for ii in 0..lnum{
            let mut maxchar = 0;
            let mut maxcount = 0;
            for cc in 0..NUM_CHARTYPE{
                let ll = aligner.letters_get_i(self.letbuff_id as usize,ii,cc);
                if maxcount < ll{
                    maxchar = cc;
                    maxcount = ll;
                }
            }
            for (cc,aa) in self.alibuff_ids.iter().enumerate(){
                if aligner.ali_get_i(*aa,ii) != maxchar{
                    ret[cc] += 1.0;
                    if maxpenal < ret[cc]{
                        topindex = cc;
                        maxpenal = ret[cc];
                    }
                }
            }
        }
        return (ret,topindex);
    }
}

pub fn read_stdins(num:usize)->Vec<String>{
    let REGEX_TAILBLANK:Regex = Regex::new(r"[\s]*$").unwrap();
    let mut ret:Vec<String> = vec![];
    for ii in 0..num{
        let mut inn:String = String::new();
        io::stdin().read_line(&mut inn).expect("Failed to read line.");
        ret.push(REGEX_TAILBLANK.replace_all(&inn,"").to_string());
    }
    return ret;
}



fn main(){

}


fn main_() {
    let start = Instant::now();
    let tss = read_stdins(1);
    let m = tss[0].parse::<usize>().unwrap();
    let mut aligner = FreqSeqAligner::new(1000,m*5);
    //let mut aligner = FreqSeqAligner::new(10,m);
    let tss = read_stdins(m);
    let mut seq_orig:Vec<FreqSequence> =vec![];
    let mut seq_chars:Vec<Vec<char>> = vec![];
    for tt in tss.into_iter(){
        seq_chars.push(tt.chars().into_iter().collect());
    }
    let mut orig_id_map:HashMap<i32,usize> = HashMap::new();
    for ii in 0..seq_chars.len(){
        seq_orig.push(FreqSequence::new(seq_chars[ii].clone(),&mut aligner));
        seq_orig[ii].set_id(ii as i32);
        orig_id_map.insert(ii as i32,ii);
    }

    let mut max_score:f32 = 0.0;
    let mut max_seq:Option<FreqSequence> = None;
    
    if true{//profile+ subalign
        let mut ids_processed:HashSet<i32> = HashSet::new();
        for ii in 0..m{
            ids_processed.insert(seq_orig[0].id);
            let start_one = start.elapsed().as_millis();
            let mut sequences = vec![];
            let slen = seq_orig.len();
            let mut maxid:i32 = (1000+slen) as i32;

            let mut children:HashMap<i32,Vec<i32>> = HashMap::new();//もとの ID からコピーされた ID へのマップ
            let mut dummy_ids:Vec<i32> = vec![];
            
            for ss in seq_orig.iter(){
                children.insert(ss.id,vec![ss.id]);
            }
            let num_subsample = m;
            //seq_orig.iter().map(|m|m.make_copy(&mut aligner)).collect();
            for ss_ in 0..num_subsample{//profile 作成用のダミー配列を追加
                if ss_ >= seq_orig.len(){
                    break;
                }
                let mut ss = seq_orig[ss_].make_copy(&mut aligner);
                let did = maxid+1;
                maxid = did;
                dummy_ids.push(did);
                children.get_mut(&ss.id).unwrap().push(did);
                ss.set_id(did);
                ss.set_first_primary_id(did);
                sequences.push(ss);
            }
            
            //profile のみ作る
            let (mut sequences,_score) =  aligner.make_msa(sequences,0.0,true);
            //panic!("{:?}",aligner.letters_buffer[sequences[0].letbuff_id as usize]);
            //alignment 用の配列を用意
            let mut copied:Vec<FreqSequence> = seq_orig.iter().map(|m|m.make_copy(&mut aligner)).collect();
            copied.reverse();
            for ss in copied.into_iter(){
                sequences.push(ss);
            }
            
            
            //テスト用ファイルは"D:\dummy\download\test_lsim_00.txt"
            //timer は STDIN の後に発動させているので、前に戻すこと
            let (mut res_,_score) =  aligner.make_msa(sequences,0.0,false);
            let mut res = res_.pop().unwrap();
            res.remove_seqs(&dummy_ids, &mut aligner);
            res.recalc_freq(&mut aligner);
            let mut prev_score = aligner.calc_score(&aligner.letters_buffer[res.letbuff_id as usize],0,res.alignment_length-1,0);
            //let (mut distmat,mut orig_id_map) = aligner.calc_alignment_dist(&res);
            
            for stepp in vec![80,40,20]{
                if stepp > res.get_alignment_length() /2 {
                    continue;
                }
                let refstep =res.get_alignment_length()/stepp;
                for zshift in vec![(0,m/3),(m/3+1,m/3*2),(m/3*2+1,m-1)]{//ゲノコン用
                    for pshift in vec![0,stepp/2]{
                        //println!("step {} {} ",stepp,pshift);
                        for pp in 0..refstep{
                            let startpos = pp*stepp+pshift;
                            let mut end = startpos+stepp;
                            if end >= res.get_alignment_length()-1-6{
                                end = res.get_alignment_length()-1;
                            }
                            
                            if end -startpos < 5{
                                break;
                            }

                            if end > res.get_alignment_length()-1{
                                end = res.get_alignment_length()-1;
                            }
                            if refstep -1 == pp{
                                end = res.get_alignment_length()-1;
                            }
                            aligner.refine_subregion(&res,startpos,end,zshift.0,zshift.1,&start);
                            
                            let duration = start.elapsed();
                            if duration.as_millis() > 9750{
                                break ;
                            }
                        }
                        let duration = start.elapsed();
                        if duration.as_millis() > 9750{
                            break ;
                        }
                    }
                    let duration = start.elapsed();
                    if duration.as_millis() > 9750{
                        break ;
                    }
                }
                let duration = start.elapsed();
                if duration.as_millis() > 9750{
                    break ;
                }
            }
            let current_score = aligner.calc_score(&aligner.letters_buffer[res.letbuff_id as usize],0,res.alignment_length-1,0);
            //println!("score: {}->{}",prev_score,current_score);
            

            if let None = max_seq{
                max_seq = Some(res);
                max_score = current_score;
            }else{
                if max_score < current_score{
                    max_seq = Some(res);
                    max_score = current_score;
                }else{
                    res.release_buffs(&mut aligner);
                }
            }
            
            let n = seq_orig.remove(0);
            seq_orig.reverse();//後ろの配列の方の方が重みとして小さく考えられているため Reverse する。
            seq_orig.push(n);
            for _ in 0..m{
                if !ids_processed.contains(&seq_orig[0].id){
                    break;
                }
                let n = seq_orig.remove(0);
                seq_orig.push(n);
            }

            let duration = start.elapsed();
            //println!("{}",duration.as_millis());
            if duration.as_millis()*2-start_one > 9750{
                break;
            }
        }
    }

    let max_seq = max_seq.unwrap_or_else(||panic!("???"));
    let mut indexx:Vec<(i32,usize)> = max_seq.primary_ids.iter().enumerate().map(|(i,e)|(*e,i)).collect();
    indexx.sort_by(|a,b|a.0.cmp(&b.0));
    for ps in indexx.iter(){
        //println!("{:?}",ps);
        let ss = &aligner.alignment_buffer[max_seq.alibuff_ids[ps.1]];
        for cc in 0..max_seq.alignment_length{
            print!("{}",ss[cc]);
        }
        println!("");
    }
    //println!("{:?}",sequences[0].letters);
}
