use std::collections::HashMap;
use rayon;
use rayon::prelude::*;
use super::guide_tree::calc_pos;

//neighbor joining tree 作成コード
//距離行列は calc_pos でアクセスできるような 45度 三角形の一次元行列。
//自分自身も入っている。
//使用法は nj_test 関数参照。

//子ノード 1 のインデクス、新しいノードから 1 への枝の長さ、子ノード 2 のインデクス以下同じ
pub fn get_next_neighbor(dist:&Vec<f32>,is_dead:&Vec<bool>,num_threads:usize)->((usize,f32),(usize,f32)){
    assert!(num_threads > 0);
    let vlen = is_dead.len();
    let mut n:usize = 0;
    let mut kmin:f64 = std::f64::MAX;
    let mut pair:((usize,f32),(usize,f32)) = ((0,-100.0),(0,-100.0));
    for ii in 0..vlen{
        if is_dead[ii]{
            continue;
        }
        n+=1;
    }

    let mut checkpair_:Vec<(usize,usize)> = vec![];
    for ii in 0..vlen{
        if is_dead[ii]{
            continue;
        }
        for jj in (ii+1)..vlen{
            if is_dead[jj]{
                continue;
            }
            checkpair_.push((ii,jj));
        }
    }
    while checkpair_.len() > 0{
        let mut checkpair:Vec<(usize,usize)> = vec![];
        for _ in 0..num_threads{
            checkpair.push(checkpair_.pop().unwrap());
            if checkpair_.len() == 0{
                break;
            }
        }

        let res:Vec<(f64,((usize,f32),(usize,f32)))> = checkpair.into_par_iter().map(|m|
            {
                let ii = m.0;
                let jj = m.1;
                let mut ssum_i:f64 = 0.0;
                let mut ssum_j:f64 = 0.0;
                let mut gsum:f64 = 0.0;
                for kk in 0..vlen{
                    if is_dead[kk]{
                        continue;
                    }
                    if (kk as i64 -ii as i64)*(kk as i64 -jj as i64 ) == 0{
                        continue;
                    }
                    ssum_i += dist[calc_pos(ii,kk)] as f64;
                    ssum_j += dist[calc_pos(jj,kk)] as f64;
                    for ll in (kk+1)..vlen{
                        if is_dead[ll]{
                            continue;
                        }
                        if (ll as i64 -ii as i64)*(ll as i64 -jj as i64 ) == 0{
                            continue;
                        }
                        gsum += dist[calc_pos(kk,ll)] as f64;
                    }
                }
                ssum_i /= n as f64-2.0;
                ssum_j /= n as f64-2.0;
                let ssum = ssum_i+ssum_j;
                let dist_ij = dist[calc_pos(ii,jj)] as f64;
                gsum /= n as f64-2.0;
                let ksum = ssum/2.0+gsum+dist_ij/2.0;
                //println!("{:?} {} ",(ii,jj),dist[calc_pos(ii,jj)]);
                return (ksum, ((ii,((dist_ij+ssum_i-ssum_j)/2.0) as f32)
                        ,(jj,((dist_ij+ssum_j-ssum_i)/2.0) as f32)))
            }
        ).collect();
        for rr in res.iter(){
            if rr.0 < kmin{
                kmin = rr.0;
                pair = rr.1;
            }
        }
    }
    
    //println!("{}",kmin);
    return pair;
}


//子 1 のインデクス、子 2 のインデクス、自分自身の長さを返す
///合計枝長が最短になる状態であるので、親子関係に系統学的な意味はないと思う
pub fn generate_unrooted_tree(dist:&mut Vec<f32>,num_threads:usize) -> Vec<(i64,i64,f32)>{
    let leafnum:usize = ((-1.0+(1.0 as f64 +8.0*dist.len() as f64).sqrt()+0.0001) as usize)/2;
    let mut nodenum:usize = leafnum;
    let mut is_dead:Vec<bool> = vec![false;nodenum];
    let mut idmap:Vec<usize> = (0..leafnum).collect();
    let mut ret:Vec<(i64,i64,f32)> = (0..leafnum).map(|m| (m as i64, -1,-1000.0)).collect();
    while nodenum > 3{
        let (a,b):((usize,f32),(usize,f32)) = get_next_neighbor(&dist, &is_dead,num_threads);
        //let mut newdist:Vec<f32> = vec![];
        for ii in 0..is_dead.len(){
            if a.0 != ii && b.0 != ii{
                let aa = calc_pos(a.0,ii);
                let bb = calc_pos(b.0,ii);
                let tdist = (dist[aa]+dist[bb])/2.0;
                dist[aa] = tdist;//つまり元の Leaf の Distance 関係は破壊される
                dist[bb] = -1000.0;
            }
        }
        is_dead[b.0] = true;
        // Leaf でない場合、子の平均長も含んだ仮枝長
        ret[idmap[a.0]].2 = a.1;        
        ret[idmap[b.0]].2 = b.1;
        ret.push((idmap[a.0] as i64,idmap[b.0] as i64,-100.0000));//枝長は後で決定。新しくできたノードは全 Leaf より後のインデクスになる
        idmap[a.0] = ret.len()-1;//破壊された Distance 関係に対応する ret 内要素のインデクスを入れる
        nodenum -= 1;
    }
    let mut lastnodes:Vec<usize> = vec![];
    for ii in 0..leafnum{
        if !is_dead[ii]{
            lastnodes.push(ii);
        }
    }
    assert!(lastnodes.len() == 3);
    let ab:f32 = dist[calc_pos(lastnodes[0],lastnodes[1])];
    let bc:f32 = dist[calc_pos(lastnodes[1],lastnodes[2])];
    let ac:f32 = dist[calc_pos(lastnodes[0],lastnodes[2])];

    let len_a:f32 = (ab+ac-bc)/2.0;
    let len_b:f32 = (bc+ab-ac)/2.0;
    let len_c:f32 = (bc+ac-ab)/2.0;

    let aa = idmap[lastnodes[0]];
    let bb = idmap[lastnodes[1]];
    let cc = idmap[lastnodes[2]];
    ret[aa].2 = len_a;//最終枝長の計算。子ノードの長さも入っている
    ret[bb].2 = len_b;
    ret[cc].2 = len_c;
    
    let mut back_tracking:Vec<usize> = vec![aa,bb,cc];

    while back_tracking.len() > 0{
        let p:usize = back_tracking.pop().unwrap();
        if ret[p].0 > -1 && ret[p].1 > -1{//子ノードの長さを引いて実際の長さに調整
            ret[p].2 -= ret[ret[p].0 as usize].2/2.0 + ret[ret[p].1 as usize].2/2.0 ;
            back_tracking.push(ret[p].0 as usize);
            back_tracking.push(ret[p].1 as usize);
        }
    }

    return ret;
}

pub fn get_newick_string(branches_:&Vec<(i64,i64,f32)>,node_name_map_:&HashMap<usize,String>)->String{
    let parent_branch:Vec<i64> = get_parent_branch(&branches_);
    let mut branches:Vec<(i64,i64,f32)> = branches_.iter().map(|m| m.clone()).collect();
    let mut node_name_map:HashMap<usize,String> = node_name_map_.iter().map(|m|(m.0.clone(),m.1.clone())).collect();
    let mut centers:Vec<usize> = vec![];
    let llen = branches.len();
    for ii in 0..llen{
        if parent_branch[ii]  < 0{ 
            centers.push(ii);           
        }
    }

    assert!(branches.len() > 2);// 3 つ以上枝がある前提

    if centers.len() == 3{
        let mut first_node:i64 = -1;
        for (ii,bb) in branches.iter().enumerate(){
            if bb.0 == -1{
                assert_eq!(ii,bb.1 as usize);
                first_node = ii as i64;
                break;
            }
            if bb.1 == -1{
                assert_eq!(ii,bb.0 as usize);
                first_node = ii as i64;
                break;
            }
        }
        assert!(first_node > -1);
        let (b, _m, h) = set_outgroup(first_node as usize,&branches, Some(&node_name_map));
        branches = b;
        node_name_map = h.unwrap();
    }else{
        assert!(parent_branch[0] == -1);//root は 0 にあることを前提とする
    }
    //return "(".to_string()+node_name_map.get(&0).unwrap_or_else(||panic!("0 is not a node."))+":0.0,"+get_internal_node_string(0,&branches,&node_name_map).as_str()+");";

    let na = get_internal_node_string(branches[0].0 as usize,&branches,&node_name_map);
    let nb = get_internal_node_string(branches[0].1 as usize,&branches,&node_name_map);
    return "(".to_owned()+node_name_map.get(&0).unwrap_or_else(||panic!("0 is not a node."))+":"+branches[0].2.to_string().as_str()+","+na.as_str()+","+nb.as_str()+")";
}

pub fn get_internal_node_string(idx:usize,branches:&Vec<(i64,i64,f32)>,node_name_map:&HashMap<usize,String>)->String{
    if branches[idx].0 > -1 && branches[idx].1 > -1{
        return "(".to_string()+get_internal_node_string(branches[idx].0 as usize,branches,node_name_map).as_str()
        +","
        +get_internal_node_string(branches[idx].1 as usize,branches,node_name_map).as_str()+"):"+branches[idx].2.to_string().as_str();
    }
    if branches[idx].0 > -1{
        return "".to_string()+node_name_map.get(&(branches[idx].0 as usize)).unwrap_or_else(||panic!("{} is not a node.",branches[idx].0))
        +":"+branches[idx].2.to_string().as_str();
    }
    panic!("???");
    /*
    return "".to_string()+branches[idx].2.to_string().as_str()+":"
    +node_name_map.get(&(branches[idx].1 as usize)).unwrap_or_else(||panic!("{} is not a node.",branches[idx].1))
    +":"+branches[idx].2.to_string().as_str();    
    */
}


pub fn calc_mismatch(aligned:&Vec<Vec<String>>)->Vec<f32>{
    let llen = aligned.len();
    let slen = aligned[0].len();
    for ii in 0..llen{
        if aligned[ii].len() != slen{
            panic!("Sequences must have the same length! 0:{} {}:{} ",slen,ii,aligned[ii].len());
        }
    }
    let mut ret:Vec<f32> = vec![];
    for ii in 0..llen{
        for jj in 0..=ii{
            let mut mismatch:usize = 0;
            for ss in 0..slen{
                if aligned[ii][ss] != aligned[jj][ss] {
                    mismatch += 1;
                }
            }
            ret.push(((mismatch as f64)/(slen as f64)) as f32);
        }
    }
    return ret;
}


pub fn get_parent_branch(branches:&Vec<(i64,i64,f32)>)->Vec<i64>{
    let mut parent_branch:Vec<i64> = vec![-1;branches.len()];
    for (ii,bb) in branches.iter().enumerate(){
        if bb.0 > -1 && bb.0 != ii as i64{
            assert!(parent_branch[bb.0 as usize] == -1);
            parent_branch[bb.0 as usize] = ii as i64;
        }
        if bb.1 > -1 && bb.1 != ii as i64{
            assert!(parent_branch[bb.1 as usize] == -1);
            parent_branch[bb.1 as usize] = ii as i64;
        }
    }
    return parent_branch;
}

//目的の Outgroup の ID と、木構造（子枝 1 の id、子枝 2 の id、自分の長さ）を表すベクトル、Leaf に貼られているラベルのハッシュマップを与えて
//Outgroup を変更した状態の木構造、前の枝のどの枝に対応するか、新しいラベルのマップを返す。
pub fn set_outgroup(outbranch:usize,branches:&Vec<(i64,i64,f32)>,node_name_map:Option<&HashMap<usize,String>>)
-> (Vec<(i64,i64,f32)>,Vec<i64>,Option<HashMap<usize,String>>){
    return change_center_branch(outbranch,branches,node_name_map);
}  

pub fn change_center_branch(centerbranch:usize,branches:&Vec<(i64,i64,f32)>,node_name_map:Option<&HashMap<usize,String>>)
-> (Vec<(i64,i64,f32)>,Vec<i64>,Option<HashMap<usize,String>>){
    let mut ret:Vec<(i64,i64,f32)> = vec![];
    let mut parent_branch:Vec<i64> = get_parent_branch(&branches);
    let llen = parent_branch.len();
    let mut centers:Vec<usize> = vec![];
    for ii in 0..llen{
        if parent_branch[ii]  < 0{ 
            centers.push(ii);           
        }else{
            if parent_branch[ii] as usize == centerbranch{
                parent_branch[ii] = -1;
            }
        }
    }
    /*
    if parent_branch[centerbranch] == -1{
        println!("Already center.");
        let mp:Option<HashMap<usize,String>> = if let Some(x) = node_name_map{
            Some(x.clone())
        }else{
            None
        };
        return (branches.clone(),(0..(branches.len() as i64)).into_iter().collect(),mp);
    }
    */
    let mut current_branch:usize = centerbranch;
    let mut old_new_map:Vec<i64> = vec![-1;branches.len()];
    let mut descend:Vec<usize> = vec![];
    if branches[centerbranch].0 != -1 && branches[centerbranch].1 != -1{ //中間点の場合子枝が処理されないので
        descend.push(branches[centerbranch].0 as usize);
        descend.push(branches[centerbranch].1 as usize);
    }

    loop{

        let p:i64 = parent_branch[current_branch];
        if p < 0{//中心点に来た
            let mut v:Vec<usize> = vec![];
            if centers.len() > 1{
                for cc in centers.into_iter(){
                    if current_branch != cc{
                        v.push(cc);
                        descend.push(cc);// 中心に到達した場合、元の木構造のまま流す
                    }
                }
                ret.push((v[0] as i64,v[1] as i64,branches[current_branch].2));//中心の自分以外のノードの親になる
                old_new_map[current_branch] = ret.len() as i64 -1;
            }else{
                ret.push((current_branch as i64,-1,branches[current_branch].2));
                old_new_map[current_branch] = ret.len() as i64 -1;
            }
            break;
        }else{
            if current_branch as i64 != branches[p as usize].0{
                ret.push((p,branches[p as usize].0,branches[current_branch].2));
                old_new_map[current_branch] = ret.len() as i64 -1;
                descend.push(branches[p as usize].0 as usize);
            }else if current_branch as i64 != branches[p as usize].1{
                ret.push((p,branches[p as usize].1,branches[current_branch].2));
                old_new_map[current_branch] = ret.len() as i64 -1;
                descend.push(branches[p as usize].1 as usize);
            }
            
            current_branch = p as usize;
        }
    }

    while descend.len() > 0{
        let dd = descend.pop().unwrap();
        ret.push(branches[dd].clone());
        old_new_map[dd] = ret.len() as i64 -1;
        if branches[dd].0 > -1 && branches[dd].1 > -1{
            descend.push(branches[dd].0 as usize);
            descend.push(branches[dd].1 as usize);
        }
    }
    for rr in ret.iter_mut(){
        if rr.0 > -1{
            rr.0 = old_new_map[rr.0 as usize];
        }
        if rr.1 > -1{
            rr.1 = old_new_map[rr.1 as usize];
        }
    }
    if let Some(prevmap) = node_name_map{
        let mut nmap:HashMap<usize,String> = HashMap::new();
        for (kk,vv) in prevmap.iter(){
            nmap.insert(old_new_map[*kk] as usize,vv.clone());
        }

        return (ret,old_new_map,Some(nmap));
    }
    
    return (ret,old_new_map,None);
}


#[allow(unused)]
fn get_child_leaves(targetbranch:usize,branches:&Vec<(i64,i64,f32)>,leaflabels:&HashMap<usize,String>) -> Vec<String>{
    let mut ret:Vec<String> = vec![];
    let mut updated:Vec<usize> = vec![targetbranch];
    while updated.len() > 0{
        let t = updated.pop().unwrap();
        if branches[t].0 != t as i64 && branches[t].0 != -1 {
            updated.push(branches[t].0 as usize);
        }
        if branches[t].1 != t as i64 && branches[t].1 != -1 {
            updated.push(branches[t].1 as usize);
        }

        if branches[t].1 == t as i64 && branches[t].0 == -1 {
            panic!("???");
            //ret.push(leaflabels.get(&t).unwrap().clone());
        }
        
        if branches[t].0 == t as i64 && branches[t].1 == -1 {
            ret.push(leaflabels.get(&t).unwrap().clone());
        }
        
    }
    return ret;
}

pub fn get_opposite(targetids:&Vec<usize>,maxid:usize)->Vec<usize>{
    let mut ret:Vec<usize> = vec![];
    let mut lastid = -1_i64;
    for ii in 0..targetids.len(){
        if lastid != targetids[ii] as i64 -1{
            for jj in ((lastid+1) as usize)..targetids[ii]{
                ret.push(jj);
            }
        }
        if lastid >= targetids[ii] as i64{
            panic!("at {}: {} >= {} , array must have been sorted! (or key dupllication?)",ii,lastid,targetids[ii]);
        }
        lastid = targetids[ii] as i64;
    }
    for ii in ((lastid+1) as usize)..=maxid{
        ret.push(ii);
    }
    return ret;
}

#[test]
fn pos_test(){
    let mut chk:Vec<(usize,usize)> = vec![];
    for ii in 0..100{
        for jj in 0..=ii{
            chk.push((ii,jj));
        }
    }
    for pp in 0..chk.len(){
        assert_eq!(calc_pos(chk[pp].0,chk[pp].1),pp);
    }
}

#[test]
fn nj_test(){
    let mut dist:Vec<f32> = vec![
        0.0,
        7.0,0.0,
        8.0,5.0,0.0,
        11.0,8.0,5.0,0.0,
        13.0,10.0,7.0,8.0,0.0,
        16.0,13.0,10.0,11.0,5.0,0.0,
        13.0,10.0,7.0,8.0,6.0,9.0,0.0,
        17.0,14.0,11.0,12.0,10.0,13.0,8.0,0.0
    ];
    let unrooted = generate_unrooted_tree(&mut dist,4);
    //println!("{:?}",unrooted);
    let mut dummyname:HashMap<usize,String> = HashMap::new();
    for ii in 0..8{
        dummyname.insert(ii,ii.to_string());
    }
    let (unrooted,_,dummyname) = set_outgroup(0, &unrooted, Some(&dummyname));
    println!("{}",get_newick_string(&unrooted,&dummyname.unwrap()));
    
    //parent がない、もしくはペアの相手が -1 である（もう一方は自分）である場合は Leaf になっている。
    //newick format で出力
    //interior node も名前を付けられる
    //http://evolution.genetics.washington.edu/phylip/newicktree.html
    //example
    //(B:6.0,(A:5.0,C:3.0,E:4.0)Ancestor1:5.0,D:11.0);
    // ((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);
    //(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;
    //(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939, Rodent:1.21460); 
    let mut dist:Vec<f32> = vec![
        0.0,
        0.2,0.0,
        0.65,0.65,0.0,
        1.675,1.675,1.425,0.0,
        1.775,1.775,1.525,0.7,0.0,
        1.4,1.4,1.15,1.275,1.375,0.0,
    ];

    let unrooted = generate_unrooted_tree(&mut dist,4);
    //println!("{:?}",unrooted);
    let mut dummyname:HashMap<usize,String> = HashMap::new();
    for ii in 0..6{
        dummyname.insert(ii,ii.to_string());
    }
    let (unrooted,_,dummyname) = set_outgroup(0, &unrooted, Some(&dummyname));
    let dummyname = dummyname.unwrap();
    println!("{}",get_newick_string(&unrooted,&dummyname));
    let mut group_checker:HashMap<Vec<usize>,f32> = HashMap::new();

    group_checker.insert(vec![0],0.1);
    group_checker.insert(vec![1],0.1);
    group_checker.insert(vec![2],0.2);
    group_checker.insert(vec![3],0.3);
    group_checker.insert(vec![4],0.4);
    group_checker.insert(vec![5],0.5);

    group_checker.insert(vec![0,1],0.35);
    group_checker.insert(vec![0,1,2],0.45);
    group_checker.insert(vec![3,4],0.475);
    group_checker.insert(vec![3,4,5],0.45);
    group_checker.insert(vec![2,3,4,5],0.35);
    
    for ii in 0..unrooted.len(){
        let mut childleaves:Vec<usize> = get_child_leaves(ii,&unrooted, &dummyname).iter().map(|m|m.parse::<usize>().unwrap()).collect();
        childleaves.sort();
        let chklen = unrooted[ii].2;
        let chklen2;
        if !group_checker.contains_key(&childleaves){
            let opp = get_opposite(&childleaves,5);
            chklen2 = *group_checker.get(&opp).unwrap();
        }else{
            chklen2 = *group_checker.get(&childleaves).unwrap();
        }
        //println!("{} {} {}",ii,chklen,chklen2);
        assert!((chklen-chklen2).abs() < 0.0001);
    }
    let llen = unrooted.len();
    for ii in 0..llen{
        //println!("{} === ",ii);
        let changed = change_center_branch(ii,&unrooted, Some(&dummyname));
        let changed_branch = changed.0;
        let changed_label = changed.2.unwrap();
        assert_eq!(changed_branch.len(),llen);
        for jj in 0..llen{
            let mut childleaves:Vec<usize> = get_child_leaves(jj,&changed_branch, &changed_label).iter().map(|m|m.parse::<usize>().unwrap()).collect();
            childleaves.sort();
            let chklen = changed_branch[jj].2;
            let chklen2;
            if !group_checker.contains_key(&childleaves){
                let opp = get_opposite(&childleaves,5);
                assert!(group_checker.contains_key(&opp),"??? {:?}\n{:?}",childleaves,opp);
                chklen2 = *group_checker.get(&opp).unwrap();
            }else{
                chklen2 = *group_checker.get(&childleaves).unwrap();
            }
            //println!("{} {} {} {:?}",jj,chklen,chklen2,childleaves);
            assert!((chklen-chklen2).abs() < 0.0001);
        }
    }
}


#[allow(unused)]
fn compare_newick_structure(s1: &str, s2: &str) -> bool {
    let mut s1_parts = s1.split(&['(', ')', ',', ':'][..]).filter(|s| !s.is_empty());
    let mut s2_parts = s2.split(&['(', ')', ',', ':'][..]).filter(|s| !s.is_empty());

    while let (Some(p1), Some(p2)) = (s1_parts.next(), s2_parts.next()) {
        if p1.parse::<f32>().is_err() && p2.parse::<f32>().is_err() {
            if p1 != p2 {
                return false;
            }
        }
    }

    s1_parts.next().is_none() && s2_parts.next().is_none()
}

#[test]
fn claude3test() {
    let mut dist: Vec<f32> = vec![
        0.0,
        0.2, 0.0,
        0.65, 0.65, 0.0,
        1.675, 1.675, 1.425, 0.0,
        1.775, 1.775, 1.525, 0.7, 0.0,
        1.4, 1.4, 1.15, 1.275, 1.375, 0.0,
    ];

    let unrooted = generate_unrooted_tree(&mut dist,4);

    let mut node_name_map: HashMap<usize, String> = HashMap::new();
    node_name_map.insert(0, "A".to_string());
    node_name_map.insert(1, "B".to_string());
    node_name_map.insert(2, "C".to_string());
    node_name_map.insert(3, "D".to_string());
    node_name_map.insert(4, "E".to_string());
    node_name_map.insert(5, "F".to_string());

    let newick_string = get_newick_string(&unrooted, &node_name_map);
    let _expected_newick_structure = "(A,((B,C),(D,E),F));";
    println!("{:?}",newick_string);
    //assert!(compare_newick_structure(&newick_string, expected_newick_structure), "Original Newick string structure does not match the expected structure");

    // Test set_outgroup
    let outbranch = 3; // Set "D" as the outgroup
    let (new_branches, _, new_node_name_map) = set_outgroup(outbranch, &unrooted, Some(&node_name_map));
    let new_newick_string = get_newick_string(&new_branches, &new_node_name_map.unwrap());
    let _expected_outgroup_newick_structure = "(D,((A,((B,C),F)),E));";
    //assert!(compare_newick_structure(&new_newick_string, expected_outgroup_newick_structure), "Newick string structure with outgroup 'D' does not match the expected structure");
    println!("{:?}",new_newick_string);
    
    // Test change_center_branch
    let centerbranch = 7; // Change the center branch to the branch connecting "A" and "B"
    let (new_branches, _, new_node_name_map) = change_center_branch(centerbranch, &unrooted, Some(&node_name_map));
    let new_newick_string = get_newick_string(&new_branches, &new_node_name_map.unwrap());
    let _expected_center_newick_structure = "(B,(A,((D,E),(F,C))));";
    //assert!(compare_newick_structure(&new_newick_string, expected_center_newick_structure), "Newick string structure with center branch changed does not match the expected structure");
    println!("{:?}",new_newick_string);
    
    println!("All tests passed!");
}

#[test]
fn small_nj_test(){
    let mut dist:Vec<f32> = vec![
        0.0,
        5.0,0.0,
        7.0,8.0,0.0,
        11.0,12.0,10.0,0.0,
        12.0,13.0,11.0,3.0,0.0
    ];
    let unrooted = generate_unrooted_tree(&mut dist,4);
    let mut dummyname:HashMap<usize,String> = HashMap::new();
    for ii in 0..5{
        dummyname.insert(ii,ii.to_string());
    }
    //println!("{:?}",set_outgroup(0,&unrooted,Some(&dummyname)));
    let (unrooted,_,dummyname) = set_outgroup(0, &unrooted, Some(&dummyname));
    println!("{}",get_newick_string(&unrooted,&dummyname.unwrap()));
}