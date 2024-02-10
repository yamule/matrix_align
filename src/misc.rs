
//b を a 基準で a3m フォーマットの 2 エントリ目以降のエントリ形式にする
//つまり a でギャップの領域について小文字になる
pub fn subject_ali_to_a3m(a:&[char],b:&[char])->Vec<String>{
    assert_eq!(a.len(),b.len());
    let mut ret:Vec<String> = vec![];
    for (aa,bb) in a.iter().zip(b.iter()){
        if *aa == '-'{
            if *bb != '-'{
                ret.push(bb.to_ascii_lowercase().to_string());
            }
        }else{
            ret.push(bb.to_ascii_uppercase().to_string());
        }
    }
    return ret;
}

//b 側の文字が小文字である場合 a 側にギャップを入れて返す
pub fn query_a3m_to_ali(a:&[char],b:&[char])->Vec<String>{
    let mut ret:Vec<String> = vec![];
    let mut apos = 0;
    let mut bpos = 0;
    while apos < a.len(){
        if b[bpos] == '-'{
            ret.push(a[apos].to_string());
            bpos += 1;
            apos += 1;
            continue;
        }
        if b[bpos].is_ascii_lowercase(){
            ret.push("-".to_owned());
            bpos += 1;
            continue;
        }
        ret.push(a[apos].to_string());
        bpos += 1;
        apos += 1;
    }
    while bpos < b.len(){
        ret.push("-".to_owned());
        bpos += 1;
    }
    return ret;
}
