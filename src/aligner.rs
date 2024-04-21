use self::matrix_process::element_multiply;
use self::matrix_process::vector_add;

use self::misc::UnionFind;

use super::*;

#[allow(unused_imports)]
use std::time::Instant;


const NUM_CHARTYPE:usize = 28;
const ACCEPTABLE_CHARS:&str = "ACDEFGHIKLMNPQRSTVWXY-";
const GAP_CHAR:char = '-';
const DIREC_UPLEFT:u8 = 0;
const DIREC_UP:u8 = 1;
const DIREC_LEFT:u8 = 2;

#[derive(Clone,Debug)]
pub enum AlignmentType {
    Local,
    Global,
}


#[derive(Clone,Debug)]
pub enum ScoreType {
    DistanceZscore,
    DotProduct,
}

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

#[derive(Clone,Debug)]
pub struct DPResult{
    pub alignment:Vec<(i32,i32)>,
    pub match_scores:Vec<f32>,
    pub score:f32
}


#[derive(Clone,Debug)]
pub struct GapPenaltyAutoAdjustParam{
    pub a1:f32,
    pub a2:f32
}


#[derive(Clone,Debug)]
pub struct ProfileAligner{
    pub dp_matrix:Vec<Vec<Vec<f32>>>,
    pub path_matrix:Vec<Vec<Vec<u8>>>,
    pub charmap:Vec<usize>,
    pub vec_size:usize,
    pub alen:usize,
    pub blen:usize,
    pub gap_open_penalty:f32,
    pub gap_extension_penalty:f32,
    pub penalty_warning:bool,
    pub alignment_type:AlignmentType,
    pub score_type:ScoreType,
    pub col_norm:bool,
    pub gap_penalty_auto_adjust:bool,
    pub auto_adjust_param:GapPenaltyAutoAdjustParam,
}
impl ProfileAligner {
    pub fn new(vec_size:usize,buff_len:usize,gap_penalty_param:Option<(f32,f32)>
        ,alignment_type:AlignmentType,score_type:ScoreType,gap_penalty_auto_adjust_:Option<GapPenaltyAutoAdjustParam>)->ProfileAligner{
        
        let mut gap_open_penalty:f32 = 0.0;
        let mut gap_extension_penalty:f32 = 0.0;
        let mut autoadjustflag = false;
        if let Some(x) = gap_penalty_param{
            if let Some(_y) = gap_penalty_auto_adjust_{
                panic!("one of gap_penalty_param or gap_penalty_auto_adjust_ should be None.");
            }
            gap_open_penalty = x.0;
            gap_extension_penalty = x.1;
        }else{
            autoadjustflag = true;
        }

        
        if let None = gap_penalty_auto_adjust_{
            if let None = gap_penalty_param{
                panic!("one of gap_penalty_param or gap_penalty_auto_adjust_ should not be None.");
            }
        }
        
        let mut gap_penalty_auto_adjust = gap_penalty_auto_adjust_.unwrap_or(
            GapPenaltyAutoAdjustParam{a1:1.0,a2:1.0});
        
        let dp_matrix:Vec<Vec<Vec<f32>>> = vec![vec![vec![];1];1];
        let path_matrix:Vec<Vec<Vec<u8>>> = vec![vec![vec![];1];1];
        let mut charmap:Vec<usize> = vec![NUM_CHARTYPE;256];
        
        let cc:Vec<char> = ACCEPTABLE_CHARS.chars().into_iter().collect();
        for ee in cc.into_iter().enumerate(){
            charmap[ee.1 as usize] = ee.0;
        }
        let col_norm = false;

        let mut ret = ProfileAligner{
            dp_matrix:dp_matrix
            ,path_matrix:path_matrix
            ,charmap
            ,vec_size:vec_size
            ,alen:0
            ,blen:0
            ,gap_open_penalty
            ,gap_extension_penalty
            ,penalty_warning:false
            ,alignment_type:alignment_type
            ,score_type:score_type
            ,col_norm:col_norm
            ,gap_penalty_auto_adjust:autoadjustflag
            ,auto_adjust_param:gap_penalty_auto_adjust
        };
        ret.reconstruct_matrix(buff_len, buff_len);
        return ret;
    }
    
    pub fn reconstruct_matrix(&mut self,amax:usize,bmax:usize){
        self.dp_matrix = vec![vec![vec![0.0;3];bmax+1];amax+1];
        self.path_matrix  = vec![vec![vec![0;3];bmax+1];amax+1];
    }

    pub fn perform_dp(&mut self,a:&SequenceProfile,b:&SequenceProfile)->DPResult {
        let mut gap_extension_penalty = self.gap_extension_penalty;
        let mut gap_open_penalty = self.gap_open_penalty;

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
        
        let mut aavec:Vec<&Vec<f32>> = vec![];
        let mut aweight:Vec<f32> = vec![];

        let mut aavec_colnorm:Vec<Vec<f32>> = vec![];
        let mut bbvec_colnorm:Vec<Vec<f32>> = vec![];

        if self.col_norm{
            for ii in 0..aalen{
                    let st = matrix_process::calc_stats(&a.gmat[ii].match_vec);
                    let mut arr = a.gmat[ii].match_vec.clone();
                    unsafe{
                        matrix_process::element_add(&mut arr,st.mean*-1.0);
                        if st.var > 0.0{
                            matrix_process::element_multiply(&mut arr,1.0/st.var);
                        }
                    }
                    aavec_colnorm.push(arr);
            }

            for ii in 0..bblen{
                let st = matrix_process::calc_stats(&b.gmat[ii].match_vec);
                let mut arr = b.gmat[ii].match_vec.clone();
                unsafe{
                    matrix_process::element_add(&mut arr,st.mean*-1.0);
                    if st.var > 0.0{
                        matrix_process::element_multiply(&mut arr,1.0/st.var);
                    }
                }
                bbvec_colnorm.push(arr);
        }
        }
        for ii in 0..aalen{
            if self.col_norm{
                aavec.push(&aavec_colnorm[ii]);
            }else{
                aavec.push(&a.gmat[ii].match_vec);
            }
            aweight.push(a.gmat[ii].match_ratio);
        }
        
        let mut bbvec:Vec<&Vec<f32>> = vec![];
        let mut bweight:Vec<f32> = vec![];
        for ii in 0..bblen{
            if self.col_norm{
                bbvec.push(&bbvec_colnorm[ii]);
            }else{
                bbvec.push(&b.gmat[ii].match_vec);
            }
            bweight.push(b.gmat[ii].match_ratio);
        }
        //バッファに入れようかと思ったが、結局新しく領域を確保していたのでやめた
        //match_ratio についてもここで渡しておくこと
        let match_score:Vec<Vec<f32>> = match self.score_type{
            ScoreType::DistanceZscore => {
                gmat::calc_dist_zscore_matrix(&aavec, &bbvec,Some(&aweight),Some(&bweight))
            },
            ScoreType::DotProduct => {
                gmat::calc_dot_product_matrix(&aavec, &bbvec)
            }
        };

        if self.gap_penalty_auto_adjust{
            let mut zmax = std::f32::NEG_INFINITY;
            let mut zmin = std::f32::INFINITY;
            for rr in 0..match_score.len(){
                for cc in 0..match_score[rr].len(){
                    zmax = zmax.max(match_score[rr][cc]);
                    zmin = zmin.min(match_score[rr][cc]);
                }
            }
            gap_open_penalty  = zmax*self.auto_adjust_param.a2+zmin*self.auto_adjust_param.a2;
            gap_extension_penalty = gap_open_penalty*0.05;
            //println!("Adjusted gap penalty open:{} extend:{}",gap_open_penalty,gap_extension_penalty);
        }


        let mut currentpenal:f32;
        
        // B 側 N 末にギャップを入れる
        for ii in 0..=aalen{
            if ii != 0{
                
                match self.alignment_type{
                    AlignmentType::Global => {
                        currentpenal = b.gmat[0].connected_ratio*gap_open_penalty+b.gmat[0].gapped_ratio*gap_extension_penalty*(ii as f32 - 1.0);
                    },
                    AlignmentType::Local => {
                        currentpenal = 0.0;
                    }
                }
                self.dp_matrix[ii][0][DIREC_LEFT as usize] = self.dp_matrix[ii-1][0][DIREC_LEFT as usize] + currentpenal*(1.0-a.gmat[ii-1].del_ratio);

                self.dp_matrix[ii][0][DIREC_UPLEFT as usize] = std::f32::NEG_INFINITY;
                self.dp_matrix[ii][0][DIREC_UP as usize] = std::f32::NEG_INFINITY;
            }
            self.path_matrix[ii][0][0] = DIREC_LEFT;
            self.path_matrix[ii][0][1] = DIREC_LEFT;
            self.path_matrix[ii][0][2] = DIREC_LEFT;
        }

        // A 側 N 末にギャップを入れる
        let mut currentpenal;
        for ii in 0..=bblen{
            if ii != 0{
                match self.alignment_type{
                    AlignmentType::Global => {
                        currentpenal = a.gmat[0].connected_ratio*gap_open_penalty+a.gmat[0].gapped_ratio*gap_extension_penalty*(ii as f32 - 1.0);
                    },
                    AlignmentType::Local => {
                        currentpenal = 0.0;
                    }
                }

                self.dp_matrix[0][ii][DIREC_UP as usize] = self.dp_matrix[0][ii-1][DIREC_UP as usize]+currentpenal*(1.0-b.gmat[ii-1].del_ratio);
                self.dp_matrix[0][ii][DIREC_UPLEFT as usize] = std::f32::NEG_INFINITY;
                self.dp_matrix[0][ii][DIREC_LEFT as usize] = std::f32::NEG_INFINITY;
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
                let abweight = (acol.match_ratio*0.5+bcol.match_ratio*0.5);
                let sc:f32 = match_score[ii-1][jj-1]*abweight;
                
                let diag_m:f32 = self.dp_matrix[ii-1][jj-1][DIREC_UPLEFT as usize] + sc;
                let diag_l:f32 = self.dp_matrix[ii-1][jj-1][DIREC_LEFT as usize] + sc;
                let diag_u:f32 = self.dp_matrix[ii-1][jj-1][DIREC_UP as usize] + sc;

                
                let lef_m:f32 = self.dp_matrix[ii-1][jj][DIREC_UPLEFT as usize]
                + (bcol.connected_ratio*gap_open_penalty + bcol.gapped_ratio*gap_extension_penalty)*abweight;
                let lef_l:f32 = self.dp_matrix[ii-1][jj][DIREC_LEFT as usize]
                + (bcol.connected_ratio*gap_extension_penalty + bcol.gapped_ratio*gap_extension_penalty)*abweight;
                let lef_u:f32 = self.dp_matrix[ii-1][jj][DIREC_UP as usize]
                + (bcol.connected_ratio*gap_open_penalty + bcol.gapped_ratio*gap_extension_penalty)*abweight;

                let up_m:f32 = self.dp_matrix[ii][jj-1][DIREC_UPLEFT as usize]
                + (acol.connected_ratio*gap_open_penalty + acol.gapped_ratio*gap_extension_penalty)*abweight;
                let up_l:f32 = self.dp_matrix[ii][jj-1][DIREC_LEFT as usize]
                + (acol.connected_ratio*gap_open_penalty + acol.gapped_ratio*gap_extension_penalty)*abweight;
                let up_u:f32 = self.dp_matrix[ii][jj-1][DIREC_UP as usize]
                + (acol.connected_ratio*gap_extension_penalty + acol.gapped_ratio*gap_extension_penalty)*abweight;

                let px = vec![
                    (DIREC_UPLEFT,(diag_m,diag_l,diag_u)),
                    (DIREC_LEFT,(lef_m,lef_l,lef_u)),
                    (DIREC_UP,(up_m,up_l,up_u)),
                ];
                //println!("{} {} {}",diag,leff,upp);
                for pp in px.iter(){
                    let poss = pp.0;
                    let (mut _m,mut _l,mut _u) = pp.1;
                    if let AlignmentType::Local = self.alignment_type{
                        _m = _m.max(0.0);
                        _l = _l.max(0.0);
                        _u = _u.max(0.0);
                    }
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
        match self.alignment_type{
            AlignmentType::Global =>{
                for ii in 0..3{
                    if maxscore < self.dp_matrix[currentx][currenty][ii]{
                        currentpos = ii as u8;
                        maxscore = self.dp_matrix[currentx][currenty][ii];
                    }
                }
            },
            AlignmentType::Local => {
                for xx in 0..=aalen{
                    for yy in 0..=bblen{
                        for ii in 0..3{
                            if maxscore < self.dp_matrix[xx][yy][ii]{
                                currentpos = ii as u8;
                                currentx = xx;
                                currenty = yy;
                                maxscore = self.dp_matrix[xx][yy][ii];
                            }
                        }
                    }
                }
            }
        }
        let mut startingx = currentx;
        let mut startingy = currenty;
        let mut nexpos = self.path_matrix[currentx][currenty][currentpos as usize];//前のセルのどこから来たか
        let mut aligned_tuple:Vec<(i32,i32)> = vec![];
        let mut ret_score:Vec<f32> = vec![];
        while currentx > 0 || currenty > 0{
            if currentpos == DIREC_UPLEFT{
                aligned_tuple.push((currentx as i32 -1 ,currenty as i32 -1));
                if currentx > 0 && currenty > 0{
                    ret_score.push(match_score[currentx as usize -1][currenty as usize -1]);
                }

                currentx -= 1;
                currenty -= 1;
                
                if let AlignmentType::Local = self.alignment_type{
                    if self.dp_matrix[currentx][currenty][nexpos as usize] <= 0.0{
                        break;
                    }
                }
            }else if currentpos == DIREC_UP{
                aligned_tuple.push((-1,currenty as i32 -1));

                currenty -= 1;
                
                
                if let AlignmentType::Local = self.alignment_type{
                    if self.dp_matrix[currentx][currenty][nexpos as usize] <= 0.0{
                        break;
                    }
                }

            }else if currentpos == DIREC_LEFT{
                aligned_tuple.push((currentx as i32 -1,-1));

                currentx -= 1;

                if let AlignmentType::Local = self.alignment_type{
                    if self.dp_matrix[currentx][currenty][nexpos as usize] <= 0.0{
                        break;
                    }
                }

            }else{
                panic!("???");
            }
            currentpos = nexpos;
            nexpos = self.path_matrix[currentx][currenty][currentpos as usize];//前のセルのどこから来たかを示す
            
        }


        //println!("{:?} vs {:?}:{}",a.headers,b.headers,maxscore);
        
        if let AlignmentType::Local = self.alignment_type{
            while currentx > 0{
                aligned_tuple.push((currentx as i32 -1,-1));
                currentx -= 1;
            }
            while currenty > 0{
                aligned_tuple.push((-1, currenty as i32 -1));
                currenty -= 1;
            }
        }
        aligned_tuple.reverse();
        
        if let AlignmentType::Local = self.alignment_type{
            startingx += 1;
            while startingx <= aalen {
                aligned_tuple.push((startingx as i32 -1,-1));
                startingx += 1;
            }
            startingy += 1;
            while startingy <= bblen{
                aligned_tuple.push((-1, startingy as i32 -1));
                startingy += 1;
            }
        }

        let mut flaga = false;
        for (ii,hh) in a.headers.iter().enumerate(){
            if hh == "seq0001" || hh == "seq0002" || hh == "seq0003"{
                flaga = true;
            }
        }

        let mut flagb = false;
        for (ii,hh) in b.headers.iter().enumerate(){
            if hh == "seq0001" || hh == "seq0002" || hh == "seq0003"{
                flagb = true;
            }
        }

        if flaga && flagb{

            for (ii,hh) in a.headers.iter().enumerate(){
                if hh == "seq0001" || hh == "seq0002" || hh == "seq0003"{
                    println!("a:{}",hh);
                    println!("{}",a.get_aligned_seq(ii).into_iter().fold("".to_owned(),|s,m|s+m.to_string().as_str()));
                }
            }
    
            for (ii,hh) in b.headers.iter().enumerate(){
                if hh == "seq0001" || hh == "seq0002" || hh == "seq0003"{
                    println!("b{}",hh);
                    println!("{}",b.get_aligned_seq(ii).into_iter().fold("".to_owned(),|s,m|s+m.to_string().as_str()));
                }
            }
    
            println!("{:?}",ret_score);
        }

        return DPResult{alignment:aligned_tuple,score:maxscore,match_scores:ret_score};
    }

    //alignment のデータをもとにプロファイルを合成して返す
    pub fn make_alignment(&mut self
        ,mut a:SequenceProfile
        ,mut b:SequenceProfile
        ,alignment:Vec<(i32,i32)>
        ,profile_only:bool // true にすると alignment の文字は a に関してのみ保持する
        ,weight:Option<(f32,f32)>)->SequenceProfile{
        assert!(!profile_only);//後で追加する
        // もし、distance のデータはプロファイルと別で持っておくならば、distance を計算した際のウエイト（カラムにおける match_ratio）の合計も持っておくこと
        // 現在はカラムの Value の平均値を使っているが、配列の長さによって分母が変わるため、DISTANCE に使う値を例えばウエイト 0.5 等で単純に合計すると値がズレる
        let anumaliseq = a.get_num_seq();
        let bnumaliseq = b.get_num_seq();
        let _numallseq = anumaliseq+bnumaliseq;
        let vec_size = a.get_vec_size();
        let alignment_length = alignment.len();

        //アラインメント情報をマージ
        let mut new_alignments:Vec<Vec<char>> = vec![];
        new_alignments.append(&mut a.member_sequences);
        new_alignments.append(&mut b.member_sequences);

        let boffset = a.alignment_mapping.len(); //b に属していたアラインメントの開始インデクス

        let mut new_alignment_mapping:Vec<Vec<(i32,i32)>> = vec![];
        new_alignment_mapping.append(&mut a.alignment_mapping);
        new_alignment_mapping.append(&mut b.alignment_mapping);

        let mut new_alignment_mapping_ids:Vec<Vec<usize>> = vec![];
        new_alignment_mapping_ids.append(&mut a.alignment_mapping_ids);

        for bb in b.alignment_mapping_ids.iter_mut(){
            for bbb in bb.iter_mut(){
                *bbb += boffset;
            }
        }

        new_alignment_mapping_ids.append(&mut b.alignment_mapping_ids);

        let mut seqmap_a:Vec<(i32,i32)> = vec![];
        let mut seqmap_b:Vec<(i32,i32)> = vec![];

        for aa in 0..alignment.len(){
            let ppos = alignment[aa];
            if ppos.0 > -1{
                seqmap_a.push((ppos.0,aa as i32));
            }
            if ppos.1 > -1{
                seqmap_b.push((ppos.1,aa as i32));
            }
        }

        let alid_a: usize = new_alignment_mapping.len();//新しく追加されたアラインメント情報のインデクス
        new_alignment_mapping.push(seqmap_a);
        new_alignment_mapping.push(seqmap_b);
        for ii in 0..new_alignment_mapping_ids.len(){
            if ii < anumaliseq{
                new_alignment_mapping_ids[ii].push(alid_a);
            }else{
                new_alignment_mapping_ids[ii].push(alid_a+1);
            }
        }

        let mut gapper:Vec<Vec<char>> = vec![vec![],vec![]];//全体ギャップ計算
        for aa in 0..alignment.len(){
            let ppos = alignment[aa];
            if ppos.0 > -1{
                gapper[0].push('X');
            }else{
                gapper[0].push('-');
            }
            if ppos.1 > -1{
                gapper[1].push('X');
            }else{
                gapper[1].push('-');
            }
        }


        let mut mergestart = 0_usize;
        let mut mergeend = gapper.len()-1;

        match self.alignment_type{
            AlignmentType::Global =>{

            },
            AlignmentType::Local => {
                for ii in 0..alignment.len(){
                    if alignment[ii].0 > -1 && alignment[ii].1 > -1{
                        mergestart = ii;
                        break;
                    }
                }
                for ii_ in 0..alignment.len(){
                    let ii = alignment.len() -1 -ii_;
                    if alignment[ii].0 > -1 && alignment[ii].1 > -1{
                        mergeend = ii;
                        break;
                    }
                }
            }
        }

        let mut headers:Vec<String> = vec![];
        headers.append(&mut a.headers);
        headers.append(&mut b.headers);

        let mut ret:SequenceProfile = SequenceProfile::new(headers.into_iter().zip(new_alignments.into_iter()).collect()
        , gapper[0].len(), vec_size,None,None,None);

        ret.alignment_mapping = new_alignment_mapping;
        ret.alignment_mapping_ids = new_alignment_mapping_ids;

        let mut aweight = if let Some(x) = weight{x.0}else{a.get_weight_sum()};
        let mut bweight = if let Some(x) = weight{x.1}else{b.get_weight_sum()};
        
        if aweight == 0.0 && bweight == 0.0{
            eprintln!("WARNING: The weight of profile 1 and 2 is 0. 0.5 was assigned.");
            aweight = 0.5;
            bweight = 0.5;
        }
        if aweight == 0.0{
            eprintln!("WARNING: The weight of profile 1 is 0. 0.0001 was assigned. The weight of profile 2 was changed to 1.0.");
            bweight = 1.0;
            aweight = 0.0001;
        }
        if bweight == 0.0{
            eprintln!("WARNING: The weight of profile 2 is 0. 0.0001 was assigned. The weight of profile 1 was changed to 1.0.");
            aweight = 1.0;
            bweight = 0.0001;
        }
        
        let mut ex_weights:Vec<(f32,f32,f32,f32)> = vec![(0.0,0.0,0.0,0.0);alignment_length+1];

        for (wei,alichar,sprof) in vec![(aweight,&gapper[0],&a),(bweight,&gapper[1],&b)]{
            let mut poscount = 0_usize;
            for alipos in 0..alignment_length{
                let mut mergeflag = match self.alignment_type{
                    AlignmentType::Global => {true},
                    AlignmentType::Local =>{
                        if alipos >= mergestart && alipos <= mergeend{
                            true
                        }else{
                            false
                        }
                    }
                };
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
                        if mergeflag{
                            gapratio += wei*sprof.gmat[poscount].gapped_ratio+ wei*sprof.gmat[poscount].connected_ratio;
                        }
                    }
                }else{
                    if alichar[alipos-1] != GAP_CHAR && alichar[alipos] != GAP_CHAR {
                        ungapratio += wei*sprof.gmat[poscount].connected_ratio;
                        gapratio += wei*sprof.gmat[poscount].gapped_ratio;
                    }else{
                        if mergeflag || alichar[alipos] != GAP_CHAR{
                            gapratio += wei*sprof.gmat[poscount].gapped_ratio+ wei*sprof.gmat[poscount].connected_ratio;
                        }
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
                    if mergeflag{
                        sum_weight_del += wei;
                    }
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

            assert!((ungapratio+gapratio) > 0.0,"{} {} {:?} {:?}",alignment_length,alipos,weight,ex_weights[alipos]);

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
                //LOCAL の場合はどんな値であろうと関係ないのでここは flag とか使わず放置
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


    pub fn make_msa_with_edge(&mut self,mut sequences: Vec<SequenceProfile>,profile_only:bool,mut edges:Vec<(usize,usize)>,merge_all:bool)
    -> Vec<Vec<(SequenceProfile,f32)>>{
        let mut uff:UnionFind = UnionFind::new(sequences.len());
        let mut bags:Vec<Option<(SequenceProfile,f32)>> = sequences.into_iter().map(|m|Some((m,0.0))).collect();
        
        edges.reverse();//pop なので Reverse
        while edges.len() > 0{
            let e_ = edges.pop().unwrap();
            let a = uff.find(e_.0);
            let b = uff.find(e_.1);
            if a == b{
                continue;
            }
            bags.push(None);
            let aseq:(SequenceProfile,f32) = bags.swap_remove(a).unwrap_or_else(||panic!("???"));
            bags.push(None);
            let bseq:(SequenceProfile,f32) = bags.swap_remove(b).unwrap_or_else(||panic!("???"));
            let mut r = self.make_msa(vec![aseq.0,bseq.0], profile_only);
            assert!(r.len() == 1);
            uff.union(a,b);
            let rres = r.pop().unwrap();
            bags[a] = Some((rres.0,rres.1));
        }

        if merge_all{
            let mut allseq:Vec<SequenceProfile> = vec![];
            let mut scores:Vec<f32> = vec![];
            for bb in bags.into_iter(){
                if let Some(x) = bb{
                    allseq.push(x.0);
                    scores.push(x.1);
                }
            }
            if allseq.len()  > 1{
                let res = self.make_msa(allseq, profile_only);
                return vec![res];
            }else{
                return vec![vec![(allseq.pop().unwrap(),scores[0])]];
            }
        }
        return bags.into_iter().filter(
                |m| if let Some(x) =  m{true}else{false}
                ).map(|m|  if let Some(x) = m{vec![(x.0,x.1)]}else{panic!("???");}).collect();
    }
    pub fn make_msa(&mut self,mut sequences: Vec<SequenceProfile>,profile_only:bool)
    -> Vec<(SequenceProfile,f32)>{
        let mut final_score:f32 = 0.0;
        sequences.reverse();
        let mut center_seq = sequences.pop().unwrap();
        let mut firstrun = true;
        while sequences.len() > 0{
            let bseq = sequences.pop().unwrap();
            let dpres;
            let newgroup;
            if firstrun{
                dpres = self.perform_dp(&center_seq,&bseq);
                newgroup = self.make_alignment(center_seq,bseq,dpres.alignment,profile_only,None);
                firstrun = false;
            }else{
                dpres = self.perform_dp(&bseq,&center_seq);
                newgroup = self.make_alignment(bseq,center_seq,dpres.alignment,profile_only,None);
            };
            final_score = dpres.score;
            
            center_seq = newgroup;
        }
        // 今のところ全部マージする
        return vec![(center_seq,final_score)];
    }

    //配列のウェイトを計算する
    //全配列の MSA が作られている前提
    pub fn calc_weights(alignments:&Vec<Vec<char>>)->Vec<f32>{
        return sequence_weighting::calc_henikoff_henikoff_weight(&alignments);
    }

}

#[derive(Clone,Debug)]
pub struct SequenceProfile{
    pub gmat:Vec<GMatColumn>,
    pub member_sequences:Vec<Vec<char>>,
    pub alignment_mapping:Vec<Vec<(i32,i32)>>, // 現在の配列の場所は次のアラインメントが行われた場合どこに移動するか
    pub alignment_mapping_ids:Vec<Vec<usize>>, //member_sequences の配列は alignment_mapping のどの mapping をどの順番で使うか
    pub headers:Vec<String>,//Fasta Header 行
    pub seq_weights:Vec<f32>,
    pub is_dummy:bool
}

impl SequenceProfile{
    pub fn new(alignments_:Vec<(String,Vec<char>)>,alignment_length:usize,vec_size:usize,seq_weights_:Option<Vec<f32>>
        ,gmat_:Option<Vec<Vec<f32>>>,gap_state_:Option<Vec<(f32,f32,f32,f32)>>)-> SequenceProfile{
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
            for xx in 0..alignment_length{
                gmat.push(GMatColumn::new(vec_size,None,None));
            }
            gmat.push(GMatColumn::new(vec_size,None,None));// ギャップ情報だけあるカラム
        }

        let mut alignments:Vec<Vec<char>> = vec![];
        let mut headers:Vec<String> = vec![];
        let mut alimap_ids:Vec<Vec<usize>> = vec![];
        for aa in alignments_.into_iter(){
            alignments.push(aa.1);
            headers.push(aa.0);
            alimap_ids.push(vec![]);
        }
        return SequenceProfile{
            member_sequences:alignments,
            alignment_mapping:vec![],
            alignment_mapping_ids:alimap_ids,
            is_dummy:false,
            headers:headers,
            gmat:gmat,
            seq_weights:seq_weights_.unwrap_or_else(||vec![1.0;seqnum])
        }
    }

    pub fn gmat_str(&self)->Vec<String>{
        let mut ret:Vec<String> = vec![];
        let achar = self.get_aligned_seq(0);
        assert_eq!(achar.len(),self.gmat.len()-1);
        for ii in 0..self.gmat.len(){
            let gg = &self.gmat[ii];
            ret.push(format!("@\t{}\t{}\t{}\t{}",gg.match_ratio,gg.del_ratio,gg.connected_ratio,gg.gapped_ratio));
            let mut pret:Vec<String> = vec![];
            if ii < self.gmat.len()-1{
                pret.push(format!("{}",achar[ii]));
                for ff in self.gmat[ii].match_vec.iter(){
                    pret.push(format!("{}",ff));
                }
                ret.push(pret.join("\t"));
            }
        }
        return ret;
    }

    pub fn create_merged_dummy(&self)->SequenceProfile{
        let alignments:Vec<Vec<char>> = vec![vec!['X';self.member_sequences[0].len()]];
        let headers:Vec<String> = vec!["dummy".to_owned()];
        let gmat = self.gmat.clone();
        let seq_weights = vec![self.get_weight_sum()];
        return SequenceProfile{
            member_sequences:alignments,
            alignment_mapping:vec![],
            alignment_mapping_ids:vec![vec![]],
            headers:headers,
            gmat:gmat,
            seq_weights:seq_weights,
            is_dummy:true
        };
    }

    pub fn get_alignment_length(&self)->usize{// 使い回す方針を取るならこの辺変える
        return self.gmat.len()-1;
    }

    pub fn get_num_seq(&self)-> usize{
        return self.member_sequences.len();
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

    //アラインされた文字列をマッピング情報から作成して返す。
    pub fn get_aligned_seq(&self,idx:usize)->Vec<char>{//a3m にし易いように char で返す。
        
        if self.alignment_mapping_ids[idx].len() == 0{
            return self.member_sequences[idx].clone();
        }
        let mut mmax = self.member_sequences[idx].len() as i32-1;
        for mapp in self.alignment_mapping_ids[idx].iter(){
            for jj in self.alignment_mapping[*mapp].iter(){
                mmax = mmax.max(jj.1);
            }
        }
        mmax += 1;

        let mut current_pos:Vec<i32> = vec![-1;mmax as usize];
        let mut updated_step_current:Vec<i32> = vec![-1;mmax as usize];
        let mut next_pos:Vec<i32> = vec![-1;mmax as usize];
        let mut updated_step_next:Vec<i32> = vec![-1;mmax as usize];
        for ii in 0..self.member_sequences[idx].len(){
            current_pos[ii] = ii as i32;
        }
        //println!("{:?}",self.alignment_mapping_ids[idx]);
        for (eii,mapp) in self.alignment_mapping_ids[idx].iter().enumerate(){
            //println!("{} {:?}",mapp,self.alignment_mapping[*mapp]);    
            for jj in self.alignment_mapping[*mapp].iter(){
                if updated_step_current[jj.0 as usize] == eii as i32 -1{
                    next_pos[jj.1 as usize] = current_pos[jj.0 as usize];
                    updated_step_next[jj.1 as usize] = eii as i32;
                }else{
                    next_pos[jj.1 as usize] = -1;
                }
            }
            let t = current_pos;
            current_pos = next_pos;
            next_pos = t;

            let t = updated_step_current;
            updated_step_current = updated_step_next;
            updated_step_next = t;
        }
        let mstep = self.alignment_mapping_ids[idx].len() -1;
        let mut ret:Vec<char> = vec!['-';mmax as usize];
        for (eii,uu) in updated_step_current.into_iter().enumerate(){
            if uu == mstep as i32{
                if current_pos[eii] > -1{
                    ret[eii] = self.member_sequences[idx][current_pos[eii] as usize];
                }
            }
        }
        return ret;
    }    
}
