
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



pub struct UnionFind{
    node_to_group:Vec<usize>
}
impl UnionFind{
    pub fn new(num_elements:usize)->UnionFind{
        return UnionFind{
            node_to_group:(0..num_elements).collect()
        };
    }
    pub fn union(&mut self,a_:usize,b_:usize)-> usize{
        let a = self.find(a_);
        let b = self.find(b_);
        let mn = a.min(b);
        let mx = a.max(b);
        self.node_to_group[mx] = mn;
        return mn;
    }
    pub fn find(&mut self,focus_index:usize) -> usize{
        let pret = self.node_to_group[focus_index];
        if pret == focus_index{
            return pret;
        }
        let pret = self.find(pret);
        self.node_to_group[focus_index] = pret;
        return pret;
    }
}


//BinaryHeap とかに使う
//BinaryHeap は数値の大きい順ということを忘れない
//https://qiita.com/hatoo@github/items/fa14ad36a1b568d14f3e
// Partial orderなものをTotal orderにする
#[derive(PartialEq, PartialOrd,Debug)]
pub struct FloatWrap<T>(pub T);

impl<T: PartialEq> Eq for FloatWrap<T> {}

impl<T: PartialOrd> Ord for FloatWrap<T> {
    fn cmp(&self, other: &FloatWrap<T>) -> std::cmp::Ordering {
        self.0.partial_cmp(&other.0).unwrap()
    }
}
