use flate2::bufread::GzDecoder;
use flate2::Compression;
use flate2::write::GzEncoder;
use std::collections::{HashMap, VecDeque};
use std::io::{BufWriter,Write,BufReader,BufRead};
use std::fs::File;
use std::vec;
use regex::Regex;
use super::aligner::SequenceProfile;
use rayon;
use rayon::prelude::*;

#[allow(dead_code)]
#[derive(Clone)]
pub struct SeqData{
    pub name:String,
    pub desc:String,
    pub seq:Vec<String>,
}
impl SeqData{
    pub fn new()->SeqData{
        return SeqData{name:"".to_string(),desc:"".to_string(),seq:Vec::new()};
    }
    pub fn create(n:String,d:String,s:String)->SeqData{
        return SeqData{name:n,desc:d,seq:SeqData::line_to_seq(s,true)};
    }
    pub fn line_to_seq(s:String,retainws:bool)->Vec<String>{
        //1 文字 1 アミノ酸を意味する文字列が与えられた場合
        if retainws {
            //pdb の ss.txt のパース目的
            let char_vec: Vec<char> = s.chars().collect();
            return char_vec.into_iter().filter(|m| m != &'\r' &&  m != &'\n' ).map(|m| m.to_string()).collect();
        }else{
            let char_vec: Vec<char> = s.chars().collect();
            return char_vec.into_iter().filter(|m| !m.is_whitespace() ).map(|m| m.to_string()).collect();
        }
    }

    /*
    pub fn filt_ws(&mut self){
        let mut vv:Vec<String> = Vec::new();
        vv.append(&mut self.seq);
        self.seq =  vv.into_iter().filter(|m| !m.is_whitespace() ).collect();
    }
    */
    pub fn load_fasta(filename:&str,retainws:bool)-> Vec<SeqData>{
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let mut ret:Vec<SeqData> = Vec::new();
        let mut seqbuff:Vec<String> = Vec::new();
        let mut currentname:String = "".to_string();
        let mut currentdesc:String = "".to_string();
        
        for (_lcount,line) in reader.lines().enumerate() {
            let line = line.unwrap();
            let rres = line.find(">");
            match rres{
                Some(x)=>{
                    //名前も配列もないものは無視される
                    if seqbuff.len() > 0 || currentname.len() > 0{
                        let mut vv:Vec<String> = Vec::new();
                        vv.append(&mut seqbuff);//中身だけ移動
                        ret.push(SeqData{name:currentname.clone(),desc:currentdesc.clone(),seq:vv});
                    }
                    if x > 0{
                        eprintln!("> was found at {}. This line was used as header anyway.",x);
                    }
                    let line = line.trim();
                    let mut nameflag = true;
                    let mut namebuff:Vec<char> = Vec::new();
                    let mut descbuff:Vec<char> = Vec::new();
                    for (sii,sss) in line.chars().enumerate(){
                        if nameflag{
                            if sii == 0 && sss == '>'{
                                continue;
                            }
                            if sss.is_whitespace(){
                                if namebuff.len()>0{
                                    nameflag = false;
                                }
                                continue;
                            }
                            namebuff.push(sss);
                        }else{
                            descbuff.push(sss);
                        }
                    }
                    currentname = namebuff.iter().collect();
                    currentdesc = descbuff.iter().collect();
                },
                _=>{
                    seqbuff.append(&mut SeqData::line_to_seq(line,retainws));
                }
            }

        }
        
        if currentname.len() > 0 ||  seqbuff.len() > 0{
            ret.push(SeqData{name:currentname,desc:currentdesc,seq:seqbuff});
        }

        return ret;
    }
}

pub fn save_lines(filename:&str,contents:Vec<String>,gzipped:bool){
    
    let file = File::create(filename).unwrap_or_else(|e| panic!("Save {} was failed! {:?}",filename,e));
        
    let mut writer: Box<dyn Write> = if gzipped {
        Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
    } else {
        Box::new(BufWriter::new(file))
    };
    for cc in contents.into_iter(){
        writer.write_all(cc.as_bytes()).unwrap_or_else(|e|panic!("{:?}",e));
        writer.write_all(b"\n").unwrap_or_else(|e|panic!("{:?}",e));
    }
    
}

pub fn load_lines(filename:&str,gzipped:bool)->Vec<String>{
    let mut ret:Vec<String> = vec![];
    let file = File::open(filename).unwrap_or_else(|e| panic!("Loading {} was failed! {:?}",filename,e));
        
    let reader: Box<dyn BufRead> = if gzipped {
        Box::new(BufReader::new(GzDecoder::new(BufReader::new(file))))
    } else {
        Box::new(BufReader::new(file))
    };
    
    for (lcount,line) in reader.lines().enumerate() {
        if let Ok(x) = line{
            ret.push(x);
        }else{
            panic!("File read error at line {}. {:?} ",lcount,line)
        }
    }
    return ret;

}

pub fn parse_gmat_block(lines:Vec<String>)-> (String,Vec<char>,Vec<Vec<f32>>,Option<Vec<(f32,f32,f32,f32)>>){
    let head_matcher:Regex = Regex::new(r">[\s]*([^\r\n]+)").unwrap();
    let mut seqname = "".to_owned();
    let mut vecsize = -1_i32;
    let mut ret_c:Vec<char> = vec![];
    let mut ret_f:Vec<Vec<f32>> = vec![];
    let mut ret_ex:Vec<(f32,f32,f32,f32)> = vec![];
    for line in lines.into_iter(){
        let ptt:Vec<String> = line.split_ascii_whitespace().map(|m|m.to_owned()).collect();
        if ptt.len() == 0 || ptt[0].len() == 0{
            continue;
        }
        let cc:Vec<char> = ptt[0].chars().into_iter().collect();
        if cc[0] == '@'{
            assert!(ptt.len() == 5);
            ret_ex.push(
                (
                    ptt[1].parse::<f32>().unwrap(),
                    ptt[2].parse::<f32>().unwrap(),
                    ptt[3].parse::<f32>().unwrap(),
                    ptt[4].parse::<f32>().unwrap()
                )
            );
            continue;
        }
        if line.starts_with(">"){
            let cres = head_matcher.captures(&line).unwrap();
            if let Some(xx) = cres.get(1){
                seqname = xx.as_str().to_string();
            }
            continue;
        }
        if cc.len() != 1{
            panic!("{} is not an expected string.",ptt[0]);
        }
        let mut val:Vec<f32> = vec![];
        if vecsize == -1{
            vecsize = ptt.len() as i32 -1;
        }else{
            assert_eq!(vecsize, ptt.len() as i32 -1);
        }
        for ii in 1..ptt.len(){
            val.push(
                ptt[ii].parse::<f32>().unwrap_or_else(|e| panic!("{:?}",e))
            );
        }
        ret_c.push(cc[0]);
        ret_f.push(val);
    }
    if ret_ex.len() > 0{
        return (seqname,ret_c,ret_f,Some(ret_ex));
    }
    return (seqname,ret_c,ret_f,None);
}
pub fn load_multi_gmat(filename:&str,gzipped:bool)-> Vec<(String,Vec<char>,Vec<Vec<f32>>,Option<Vec<(f32,f32,f32,f32)>>)>{
    //名前、一文字アミノ酸、何らかのベクトル、ギャップなどの補足情報
    let mut ret:Vec<(String,Vec<char>,Vec<Vec<f32>>,Option<Vec<(f32,f32,f32,f32)>>)> = vec![];
    let file = File::open(filename).unwrap_or_else(|e| panic!("Loading {} was failed! {:?}",filename,e));
        
    let reader: Box<dyn BufRead> = if gzipped {
        Box::new(BufReader::new(GzDecoder::new(BufReader::new(file))))
    } else {
        Box::new(BufReader::new(file))
    };
    
    let mut linebuff:Vec<String> = vec![];
    for (_lcount,line) in reader.lines().enumerate() {
        if let Ok(x) = line{
            if x.starts_with("#"){ //コメント行
                continue;
            }
            if x.starts_with(">"){
                if linebuff.len() > 0{
                    let pres = parse_gmat_block(linebuff);
                    ret.push(pres);
                }

                linebuff = vec![];
            }
            if x.starts_with("\r") || x.starts_with("\n") || x.len() == 0{
                continue;
            }
            linebuff.push(x);
        }
    }
    if linebuff.len() > 0{
        ret.push(parse_gmat_block(linebuff));
    }
    return ret;
}

//max_samples を超えるまで gmat を読み込んで返す。一つのファイルに複数の gmat が保存されている場合 max_samples を超えることがある。
pub fn parallel_load_multi_gmat(filenames:&mut VecDeque<String>,max_samples:usize,num_threads:usize)-> Vec<(String,Vec<char>,Vec<Vec<f32>>,Option<Vec<(f32,f32,f32,f32)>>)>{
    let mut ret:Vec<(String,Vec<char>,Vec<Vec<f32>>,Option<Vec<(f32,f32,f32,f32)>>)> = vec![];
    while filenames.len() > 0{
        let mut loader_minibatch:Vec<(usize,String)> = vec![];
        while filenames.len() >0{ 
            let ii = filenames.pop_front().unwrap();
            loader_minibatch.push((loader_minibatch.len(),ii));
            if loader_minibatch.len() >= num_threads{
                break;
            }
        }
        let mut res:Vec<(usize,Vec<(String,Vec<char>,Vec<Vec<f32>>,Option<Vec<(f32,f32,f32,f32)>>)>)> = loader_minibatch.into_par_iter().map(|v|{
            let fileindex = v.0;
            let filename = v.1;
            let pres = load_multi_gmat(&filename,filename.ends_with(".gz"));
            return (fileindex,pres);
        }).collect();
        res.sort_by(|a,b| (a.0.cmp(&b.0)));
        for mut rr in res.into_iter(){
            ret.append(&mut rr.1);
        }
        if ret.len() >= max_samples{
            break;
        }
    }
    return ret;

}

pub fn fasta_to_gmat(seq:&Vec<SeqData>,smat:&HashMap<String,Vec<f32>>)-> Vec<(String,Vec<char>,Vec<Vec<f32>>,Option<Vec<(f32,f32,f32,f32)>>)>{
    let mut gmat1_:Vec<(String,Vec<char>,Vec<Vec<f32>>,Option<Vec<(f32,f32,f32,f32)>>)> = vec![];
    for ss in seq.iter(){
        let name = ss.name.clone();
        let mut values:Vec<Vec<f32>> = vec![];
        let mut cc:Vec<char> = vec![];
        let mut cs:Vec<String> = vec![];// 単に char 化して to_string するのが面倒かと思って作っているだけ
        for sss in ss.seq.iter(){
            if sss.len() > 1{
                eprintln!("Only single char is allowed for representation. {} was changed to X.",sss);
                cc.push('X');
            }else if !smat.contains_key(sss){
                eprintln!("The value for {} is not defined. {} was changed to X.",sss,sss);
                cc.push('X');
            }else{
                let c:Vec<char> = sss.chars().into_iter().collect();
                cc.push(c[0]);
                cs.push(sss.clone());
            }
        }
        for c in cs.iter(){
            values.push(smat.get(c).unwrap().clone());
        }
        gmat1_.push((
            name,cc,values,None
        ));
    }
    return gmat1_;
}

pub fn save_gmat(filename:&str,seqs:&Vec<(usize,&SequenceProfile)>,gzipped:bool){
    let mut buff:Vec<String> = vec![];
    for ss in seqs.iter(){
        let baseindex = ss.0;
        let proff = ss.1;
        buff.push(
            ">".to_owned()+proff.headers[baseindex].as_str()
        );
        let aliseq = proff.get_aligned_seq(baseindex);
        let alilen = proff.get_alignment_length();
        for ii in 0..alilen{
            buff.push(
                format!("@\t{}\t{}\t{}\t{}",proff.gmat[ii].match_ratio,proff.gmat[ii].del_ratio,proff.gmat[ii].match_to_del,proff.gmat[ii].del_to_del)
            );
            let mut ptt:Vec<String> = vec![];
            if ii < aliseq.len(){
                ptt.push(
                    aliseq[ii].to_string()
                );
            }else{
                ptt.push("-".to_owned());
            }
            for vv in proff.gmat[ii].match_vec.iter(){
                ptt.push(format!("{:.5}",vv));
            }
            buff.push(
                ptt.join("\t")
            );
        }
        buff.push(
            format!("@\t{}\t{}\t{}\t{}",proff.gmat[alilen].match_ratio,proff.gmat[alilen].del_ratio,proff.gmat[alilen].match_to_del,proff.gmat[alilen].del_to_del)
        );
    }
    save_lines(filename,buff, gzipped);
}