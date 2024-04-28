use std::collections::{HashMap,HashSet};




pub struct SimpleArgParse{
    items:HashMap<String,String>,
    order:Vec<String>,
    mapper:HashMap<String,String>, //基準文字列もしくは省略形->基準文字列へのマップ
    mapped_rev:HashMap<String,String>,//基準文字列→省略形へのマップ。省略形が無い場合は None
    required:HashSet<String>,
    user_defined:HashSet<String>,//デフォルトでなくユーザーが設定したもの
    acceptable_values:HashMap<String,HashSet<String>>,
    descriptions:Vec<String>
}
impl SimpleArgParse{
    pub fn new(mut itemlist:Vec<(&str,Option<&str>,&str,Option<&str>,Vec<&str>,bool)>)->SimpleArgParse{

        let mut items:HashMap<String,String> = HashMap::new();
        let mut required:HashSet<String> = HashSet::new();
        let mut mapper:HashMap<String,String> = HashMap::new();
        let mut mapped_rev:HashMap<String,String> = HashMap::new();
        let mut descriptions:Vec<String> = vec![];
        let mut acceptable_values:HashMap<String,HashSet<String>> = HashMap::new();
        let mut order:Vec<String> = vec![];
        
        for ii in itemlist.iter(){
            let name = ii.0;
            if name == "--help"{
                panic!("--help is used for printing help message. So it is not acceptable.");
            }
        }
        
        itemlist.push(("--help",Some("-h")
        ,": Print help message."
        ,None,vec![],false));
        
        for ii in itemlist.iter(){
            let name = ii.0;
            let altname = ii.1;
            let descstring = ii.2;
            let defvalue = ii.3;
            let acc = &ii.4;
            let req = ii.5;

            order.push(name.to_string());

            let mut desc = "".to_owned();
            desc += name;
            if let Some(x) = altname{
                desc += ", ";
                desc += x;
            }
            desc += ":\t";
            desc += descstring;
            if let Some(x) = defvalue{
                desc += format!(" Default {}.",x).as_str();
            }
            if acc.len() > 0{
                desc += format!(" Acceptable values are {}.",acc.join(",")).as_str();
            }
            if req{
                desc += " Required."
            }
            descriptions.push(desc);
            
            if req{
                required.insert(name.to_string());
            }

            if mapper.contains_key(name){ //重複チェックも兼ねている
                panic!("Duplicated key {}.",name);
            }

            mapper.insert(name.to_string(), name.to_string());
            if let Some(x) = altname{
                if mapper.contains_key(x){
                    panic!("Duplicated key {}.",x);
                }
                mapper.insert(x.to_string(), name.to_string());
                mapped_rev.insert(name.to_string(),x.to_string());
            }
            
            if acc.len() > 0{
                let mut hs:HashSet<String> = HashSet::new();
                for a in acc.iter(){
                    hs.insert(a.to_string());
                }
                acceptable_values.insert(name.to_string(),hs);
            }

            if let Some(x) = defvalue{
                items.insert(name.to_string(),x.to_string());
            }
        }

        return SimpleArgParse{
            items:items,
            order:order,
            required:required,
            mapper:mapper,
            mapped_rev:mapped_rev,
            acceptable_values:acceptable_values,
            user_defined:HashSet::new(),
            descriptions:descriptions
        };
    }

    //単語が入っている状態の文字列の Vec を与えて key-value のペアに成型する
    fn preprocess(&self,args_:Vec<String>)->HashMap<String,String>{
        let mut args:Vec<String> = vec![];
        for aa in args_{
            if self.mapper.contains_key(&aa){
                args.push(self.mapper.get(&aa).unwrap().clone());
            }else{
                args.push(aa);
            }
        }
        assert!(args.len() > 0);
        let fst = args.remove(0);
        assert!(fst.contains("matrix_align"));
        let mut ret:HashMap<String,String> = HashMap::new();
        let mut non_key_count:usize = 0;
        let mut ii:usize = 0;
        while ii < args.len(){
            if args[ii].starts_with("--"){
                assert!(!ret.contains_key(&args[ii]),"{} has already been assigned.",args[ii]);
                if ii == args.len()-1 || args[ii+1].starts_with("--"){
                    ret.insert(args[ii].clone(),"true".to_owned());
                }else{
                    //
                    ret.insert(args[ii].clone(),args[ii+1].clone()); 
                    ii += 1;
                }
            }else{
                ret.insert(format!("--nokey_{}",non_key_count),args[ii].clone());
                non_key_count += 1;
            }
            ii += 1;
        }
        return ret;
    }

    pub fn print_help(&self){
        for vv in self.descriptions.iter(){
            println!("{}",vv);
        }
    }

    pub fn eprint_help(&self){
        for vv in self.descriptions.iter(){
            eprintln!("{}",vv);
        }
    }

    pub fn print_items(&self){
        for ss in self.order.iter(){
            if ss == "--help"{
                continue;
            }
            if let Some(x) = self.items.get(ss){
                println!("{}: {}",ss,x);
            }else{
                println!("{}: None",ss);
            }
        }
    }

    pub fn parse(&mut self,a:Vec<String>){
        let preprocessed = self.preprocess(a);
        let mut errors:Vec<String> = vec![];
        let mut processed:HashSet<String> = HashSet::new();
        for (kk,vv) in preprocessed.iter(){
            self.user_defined.insert(kk.to_string());
            if !self.mapper.contains_key(kk){
                if kk.starts_with("--nokey"){
                    errors.push(
                        format!("{} lacks tag.",vv)
                    );
                }else{
                    errors.push(
                        format!("{} is not an acceptable option.",kk)
                    );
                }
            }else{
                match self.acceptable_values.get(kk){
                    Some(x) => {
                        if !x.contains(vv){
                            panic!("{} is not expected for {}.",vv,kk);
                        }
                    },
                    _ => {} 
                }
                if processed.contains(kk){
                    if self.mapped_rev.contains_key(kk){
                        errors.push(
                            format!("{} or {} were assigned multiple times.",kk,self.mapped_rev.get(kk).unwrap())
                        );
                    }else{
                        errors.push(
                            format!("{} was assigned multiple times.",kk)
                        );
                    }
                }
                processed.insert(kk.to_string());
                self.items.insert(kk.to_string(),vv.to_string());
            }
        }
        if processed.contains("--help"){
            self.print_help();
            std::process::exit(0);
        }
        for rr in self.required.iter(){
            if !self.items.contains_key(rr){
                errors.push(format!("{} is required.",rr));
            }
        }
        if errors.len() > 0{
            self.eprint_help();
            eprintln!("------");
            for ee in errors.iter(){
                eprintln!("{}",ee);
            }
            panic!();
        }
    }

    pub fn check_key(&self,k:&str){
        if !self.mapper.contains_key(k){
            panic!("Error in code. {} is not registered as an option.",k);
        }
    }

    pub fn is_generous_false(&self,k:&str)->bool{
        self.check_key(k);
        if !self.items.contains_key(k){
            return true;
        }
        let chk = self.items.get(k).unwrap().clone().to_lowercase();
        if chk == "0" || chk == "false" || chk == "none"|| chk == "nulla"{
            return true;
        }
        return false;
    }

    pub fn user_defined(&self,k:&str)-> bool{
        self.check_key(k);
        return self.user_defined.contains(k);
    }

    pub fn get_string(&self,k:&str)->Option<String>{
        self.check_key(k);
        if let None = self.items.get(k){
            return None;
        }
        return Some(self.items.get(k).unwrap().clone());
    }
    
    pub fn get_int(&self,k:&str)->Option<i128>{
        self.check_key(k);
        if let None = self.items.get(k){
            return None;
        }
        let ret:i128 = self.items.get(k).unwrap().parse::<i128>()
        .expect(&format!("Failed to parse the value of '{}' '{}' to int.",k,self.items.get(k).unwrap()));
        return Some(ret);
    }
    
    pub fn get_float(&self,k:&str)->Option<f64>{
        self.check_key(k);
        if let None = self.items.get(k){
            return None;
        }
        let ret:f64 = self.items.get(k).unwrap().parse::<f64>()
        .expect(&format!("Failed to parse the value of '{}' '{}' to float.",k,self.items.get(k).unwrap()));
        return Some(ret);
    }
    
    pub fn get_bool(&self,k:&str)->Option<bool>{
        self.check_key(k);
        if let None = self.items.get(k){
            return None;
        }
        let chk = self.items.get(k).unwrap().clone().to_lowercase();
        if chk == "1" || chk == "true"{
            return Some(true);
        }
        if chk == "0" || chk == "false"{
            return Some(false);
        }
        panic!(
            "Failed to parse the value of '{}' '{}' as a boolean. Acceptable values are 0, 1, false, or true.",
            k,self.items.get(k).unwrap()
        );
    }
}