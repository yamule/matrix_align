use std::collections::HashMap;

//neighbor joining tree 作成コード
//距離行列は calc_pos でアクセスできるような 45度 三角形の一次元行列。
//自分自身も入っている。
//使用法は nj_test 関数参照。
pub fn calc_pos(aa:usize,bb:usize) -> usize{
    let mut a:f64 = aa as f64;
    let mut b:f64 = bb as f64;
    if a < b{// a の方が常に大きい
        let c = b;
        b = a;
        a = c;
    }
    return (a*a/2.0+a/2.0+b +0.000001) as usize;
}

//子ノード 1 のインデクス、新しいノードから 1 への枝の長さ、子ノード 2 のインデクス以下同じ
pub fn get_next_neighbor(dist:&Vec<f32>,is_dead:&Vec<bool>)->((usize,f32),(usize,f32)){
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
    for ii in 0..vlen{
        if is_dead[ii]{
            continue;
        }
        for jj in (ii+1)..vlen{
            if is_dead[jj]{
                continue;
            }
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
            if ksum < kmin{
                kmin = ksum;
                pair = ((ii,((dist_ij+ssum_i-ssum_j)/2.0) as f32)
                       ,(jj,((dist_ij+ssum_j-ssum_i)/2.0) as f32));
            }
        }
    }
    //println!("{}",kmin);
    return pair;
}


pub fn generate_unrooted_tree(dist:&mut Vec<f32>) -> Vec<(i64,i64,f32)>{
    let leafnum:usize = ((-1.0+(1.0 as f64 +8.0*dist.len() as f64).sqrt()+0.0001) as usize)/2;
    let mut nodenum:usize = leafnum;
    let mut is_dead:Vec<bool> = vec![false;nodenum];
    let mut idmap:Vec<usize> = (0..leafnum).collect();
    let mut ret:Vec<(i64,i64,f32)> = (0..leafnum).map(|m| (m as i64, -1,-1000.0)).collect();
    while nodenum > 3{
        let (a,b):((usize,f32),(usize,f32)) = get_next_neighbor(&dist, &is_dead);
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

    //新規の枝の長さの計算//木の生成後に後ろからできるはず
    //配列中に挿入
    //3 になったらBreak
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
        let (b,h) = set_outgroup(first_node as usize,&branches, &node_name_map);
        branches = b;
        node_name_map = h;
    }

    return "(".to_string()+node_name_map.get(&0).unwrap_or_else(||panic!("0 is not a node."))+":"+branches[0].2.to_string().as_str()+","+get_internal_node_string(0,&branches,&node_name_map).as_str()+");";
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


pub fn set_outgroup(outbranch:usize,branches:&Vec<(i64,i64,f32)>,node_name_map:&HashMap<usize,String>)
-> (Vec<(i64,i64,f32)>,HashMap<usize,String>){
    let mut ret:Vec<(i64,i64,f32)> = vec![];
    let parent_branch:Vec<i64> = get_parent_branch(&branches);
    let llen = parent_branch.len();
    let mut centers:Vec<usize> = vec![];
    for ii in 0..llen{
        if parent_branch[ii]  < 0{ 
            centers.push(ii);           
        }
    }
    let mut current_branch:usize = outbranch;
    let mut old_new_map:Vec<i64> = vec![-1;branches.len()];
    let mut descend:Vec<usize> = vec![];
    loop{

        let p:i64 = parent_branch[current_branch];
        if p < 0{
            let mut v:Vec<usize> = vec![];
            if centers.len() > 1{
                for cc in centers.into_iter(){
                    if current_branch != cc{
                        v.push(cc);
                        descend.push(cc);
                    }
                }
                ret.push((v[0] as i64,v[1] as i64,branches[current_branch].2));
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
    let mut nmap:HashMap<usize,String> = HashMap::new();
    for (kk,vv) in node_name_map.iter(){
        nmap.insert(old_new_map[*kk] as usize,vv.clone());
    }
    return (ret,nmap);
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
    let unrooted = generate_unrooted_tree(&mut dist);
    println!("{:?}",unrooted);
    let mut dummyname:HashMap<usize,String> = HashMap::new();
    for ii in 0..8{
        dummyname.insert(ii,ii.to_string());
    }
    println!("{:?}",set_outgroup(0,&unrooted,&dummyname));
    println!("{}",get_newick_string(&unrooted,&dummyname));
    
    //parent がない、もしくはペアの相手が -1 である（もう一方は自分）である場合は Leaf になっている。
    //newick format で出力
    //interior node も名前を付けられる
    //http://evolution.genetics.washington.edu/phylip/newicktree.html
    //example
    //(B:6.0,(A:5.0,C:3.0,E:4.0)Ancestor1:5.0,D:11.0);
    // ((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);
    //(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;
    //(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939, Rodent:1.21460); 
}