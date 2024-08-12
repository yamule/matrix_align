use std::collections::HashMap;

pub fn load_from_lines(lines:Vec<String>)-> HashMap<String,Vec<f32>>{
    let mut ret:HashMap<String,Vec<f32>> = HashMap::new();
    for (_ii,line) in lines.into_iter().enumerate() {
        let bs = line.trim();
        if bs.starts_with("#"){
            continue;
        }
        let ptt:Vec<String> = bs.split_whitespace().map(|m| m.to_string()).collect();
        let mut values:Vec<f32> = vec![];
        for ll in 1..ptt.len(){
            let ss = match ptt[ll].parse::<f32>(){
                Ok(x)=>{
                    x
                },
                _=>{
                    eprintln!("{} can not be parsed! zero was assigned",ptt[ll]);
                    0.0
                }
            };
            values.push(ss);
        }
        assert!(!ret.contains_key(&ptt[0]));
        ret.insert(ptt[0].clone(),values);
    }
    return ret;
}

pub fn get_blosum62_matrix()->HashMap<String,Vec<f32>>{
    let mut lin:Vec<String> = Vec::new();
    //https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
    //#  Matrix made by matblas from blosum62.iij
    //#  * column uses minimum score
    //#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
    //#  Blocks Database = /data/blocks_5.0/blocks.dat
    //#  Cluster Percentage: >= 62
    //#  Entropy =   0.6979, Expected =  -0.5209

    // removed BZX* columns
    lin.push("#  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V".to_string());
    lin.push("A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0".to_string());
    lin.push("R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3".to_string());
    lin.push("N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3".to_string());
    lin.push("D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3".to_string());
    lin.push("C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1".to_string());
    lin.push("Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2".to_string());
    lin.push("E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2".to_string());
    lin.push("G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3".to_string());
    lin.push("H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3".to_string());
    lin.push("I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3".to_string());
    lin.push("L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1".to_string());
    lin.push("K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2".to_string());
    lin.push("M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1".to_string());
    lin.push("F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1".to_string());
    lin.push("P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2".to_string());
    lin.push("S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2".to_string());
    lin.push("T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0".to_string());
    lin.push("W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3".to_string());
    lin.push("Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1".to_string());
    lin.push("V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4".to_string());
    lin.push("B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3".to_string());
    lin.push("Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2".to_string());
    lin.push("X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1".to_string());

    return load_from_lines(lin);
}
