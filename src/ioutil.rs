use flate2::bufread::GzDecoder;
use std::io::{Read,BufWriter,Write,BufReader,BufRead};
use std::fs::File;
use std::collections::HashMap;
use std::vec;

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

pub fn load_pssm_matrix(filename:&str,gzipped:bool)-> (Vec<char>,Vec<Vec<f32>>){
    
    let mut ret_c:Vec<char> = vec![];
    let mut ret_f:Vec<Vec<f32>> = vec![];
    
    let file = File::open(filename).unwrap_or_else(|e| panic!("Loading {} was failed! {:?}",filename,e));
        
    let reader: Box<dyn BufRead> = if gzipped {
        Box::new(BufReader::new(GzDecoder::new(BufReader::new(file))))
    } else {
        Box::new(BufReader::new(file))
    };
    
    let mut vecsize = -1_i32;
    for (lcount_,line) in reader.lines().enumerate() {
        let lcount = lcount_+1;
        if let Ok(x) = line{
            let ptt:Vec<String> = x.split_ascii_whitespace().map(|m|m.to_owned()).collect();
            if ptt.len() < 2{
                eprintln!("line:{} {:?} was skipped.",lcount,x);
                continue;
            }

            let cc:Vec<char> = ptt[0].chars().into_iter().collect();
            if cc[0] == '#'{ //コメント行
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
                    ptt[ii].parse::<f32>().unwrap_or_else(|e| panic!("line:{} {}",lcount,e))
                );
            }
            ret_c.push(cc[0]);
            ret_f.push(val);
        }
    }
    return (ret_c,ret_f);
}


#[test]
fn pssm_load_check(){
    //単にエラーが出ないかだけ
    let v = load_pssm_matrix("./example_files/esm2_650m_example_output/d1g43a_.res.gz",true);
    for vv in v.0.into_iter().zip(v.1.into_iter()){
        let mut arr:Vec<f32> = vec![];
        for ii in 0..10{
            arr.push(vv.1[ii]);
        }
        //println!("{} {} {:?}",vv.0,vv.1.len(),arr);
    }
}