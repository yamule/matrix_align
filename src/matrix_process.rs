use std::arch::x86_64::*;


//二つの一次元配列の内積を取る
#[cfg(all(not(target_feature = "avx2"),target_feature = "sse3"))]
pub unsafe fn dot_product(a: &[f32], b: &[f32]) -> f32 {
    assert_eq!(a.len(), b.len());
    assert!(a.len() % 4 == 0);

    //eprintln!("running sse3");
    let mut sum = 0.0;
    for i in (0..a.len()).step_by(4) {
        // SIMD型にデータをロード
        let a_chunk = _mm_loadu_ps(a[i..].as_ptr());
        let b_chunk = _mm_loadu_ps(b[i..].as_ptr());

        // 乗算と加算
        let mul = _mm_mul_ps(a_chunk, b_chunk);

        //[a[0]+a[1],a[2]+a[3],b[0]+b[1],b[2]+b[3]] を返す
        let sum_chunk = _mm_hadd_ps(mul, mul);

        // スカラーに変換して合計
        
        sum += _mm_cvtss_f32(sum_chunk) + 
        _mm_cvtss_f32(_mm_shuffle_ps(sum_chunk, sum_chunk, 0x55));//0x55 は 01010101 と展開されつまり結果の配列については a[1],a[1],b[1],b[1] となる。
    
    }
    sum
}


#[cfg(target_feature = "avx2")]
pub unsafe fn dot_product(a: &[f32], b: &[f32]) -> f32 {
    assert_eq!(a.len(), b.len());
    assert!(a.len() % 8 == 0);
    //eprintln!("running avx2");
    let mut sum = 0.0;
    let mut sum_vec = _mm256_setzero_ps();
    for i in (0..a.len()).step_by(8) {
        let a_chunk = _mm256_loadu_ps(a[i..].as_ptr());
        let b_chunk = _mm256_loadu_ps(b[i..].as_ptr());
        let mul_chunk = _mm256_mul_ps(a_chunk, b_chunk);
        sum_vec = _mm256_add_ps(sum_vec, mul_chunk);
    }
    // 水平加算
    sum_vec = _mm256_hadd_ps(sum_vec, sum_vec);
    sum_vec = _mm256_hadd_ps(sum_vec, sum_vec);
    let mut temp = [0.0; 8];
    _mm256_storeu_ps(temp.as_mut_ptr(), sum_vec);
    sum = temp[0] + temp[4];
    sum
}


#[cfg(all(not(target_feature = "sse3"),not(target_feature = "avx2")))]
pub fn dot_product(a: &[f32], b: &[f32]) -> f32 {
    assert_eq!(a.len(), b.len());
    //eprintln!("running native");
    return a.iter().zip(b.iter()).map(|(x,y)|x*y).sum();
}





type SimdOpAvx2 = unsafe fn(__m256, __m256) -> __m256;
type SimdOpSse3 = unsafe fn(__m128, __m128) -> __m128;

// 二つの一次元配列を与えて min, max, add（sum）等を行い、第一要素の配列に代入する
#[cfg(all(not(target_feature = "avx2"),target_feature = "sse3"))]
pub unsafe fn element_op(vec1:&mut [f32],val:f32,op: SimdOpSse3){
    let val_vec = _mm_set1_ps(val);
    let mut i = 0;
    while i + 3 < vec1.len() {
        let chunk1 = _mm_loadu_ps(vec1.as_ptr().add(i));
        let result = op(chunk1, val_vec);
        _mm_storeu_ps(vec1.as_mut_ptr().add(i), result);
        i += 4;
    }
}

#[cfg(all(not(target_feature = "avx2"),target_feature = "sse3"))]
pub unsafe fn element_max(vec1:&mut [f32],val:f32){
    element_op(vec1,val,_mm_max_ps);
    // ベクタのサイズが4の倍数でない場合の処理
    if vec1.len()%4 != 0{
        let i = vec1.len() - vec1.len()%4;
        for j in i..vec1.len() {
            vec1[j] = vec1[j].max(val);
        }
    }
}

#[cfg(all(not(target_feature = "avx2"),target_feature = "sse3"))]
pub unsafe fn element_min(vec1:&mut [f32],val:f32){
    element_op(vec1,val,_mm_min_ps);
    // ベクタのサイズが4の倍数でない場合の処理
    if vec1.len()%4 != 0{
        let i = vec1.len() - vec1.len()%4;
        for j in i..vec1.len() {
            vec1[j] = vec1[j].min(val);
        }
    }
}

#[cfg(all(not(target_feature = "avx2"),target_feature = "sse3"))]
pub unsafe fn element_add(vec1:&mut [f32],val:f32){
    element_op(vec1,val,_mm_add_ps);
    // ベクタのサイズが4の倍数でない場合の処理
    if vec1.len()%4 != 0{
        let i = vec1.len() - vec1.len()%4;
        for j in i..vec1.len() {
            vec1[j] = vec1[j] + val;
        }
    }
}

//一次元配列をスカラーで乗算する
#[cfg(all(not(target_feature = "avx2"),target_feature = "sse3"))]
pub unsafe fn element_multiply(vec1: &mut [f32], val: f32) {
    element_op(vec1,val,_mm_mul_ps);
    if vec1.len()%4 != 0{
        let i = vec1.len() - vec1.len()%4;
        // ベクタのサイズが4の倍数でない場合の処理
        for j in i..vec1.len() {
            vec1[j] *= val;
        }
    }
}

#[cfg(target_feature = "avx2")]
pub unsafe fn element_op(vec1:&mut [f32],val:f32,op:SimdOpAvx2){
    let val_vec = _mm256_set1_ps(val);
    let mut i = 0;
    while i + 7 < vec1.len() {
        let chunk1 = _mm256_loadu_ps(vec1.as_ptr().add(i));
        let result = op(chunk1,val_vec);
        _mm256_storeu_ps(vec1.as_mut_ptr().add(i), result);
        i += 8;
    }
}

#[cfg(target_feature = "avx2")]
pub unsafe fn element_max(vec1:&mut [f32],val:f32){
    element_op(vec1,val,_mm256_max_ps);
    if vec1.len()%8 != 0{
        let i = vec1.len() - vec1.len()%8;
        for j in i..vec1.len() {
            vec1[j] = vec1[j].max(val);
        }
    }
}

#[cfg(target_feature = "avx2")]
pub unsafe fn element_min(vec1:&mut [f32],val:f32){
    element_op(vec1,val,_mm256_min_ps);
    if vec1.len()%8 != 0{
        let i = vec1.len() - vec1.len()%8;
        for j in i..vec1.len() {
            vec1[j] = vec1[j].min(val);
        }
    }
}

#[cfg(target_feature = "avx2")]
pub unsafe fn element_add(vec1:&mut [f32],val:f32){
    element_op(vec1,val,_mm256_add_ps);
    if vec1.len()%8 != 0{
        let i = vec1.len() - vec1.len()%8;
        for j in i..vec1.len() {
            vec1[j] = vec1[j] + val;
        }
    }
}

#[cfg(target_feature = "avx2")]
pub unsafe fn element_multiply(vec1: &mut [f32], val: f32) {
    element_op(vec1,val,_mm256_mul_ps);
    if vec1.len()%8 != 0{
        let i = vec1.len() - vec1.len()%8;
        for j in i..vec1.len() {
            vec1[j] *= val;
        }
    }
}


#[cfg(all(not(target_feature = "sse3"),not(target_feature = "avx2")))]
pub fn element_min(vec1: &mut [f32], val: f32) {
    for elem in vec1 {
        *elem = elem.min(val);
    }
}
#[cfg(all(not(target_feature = "sse3"),not(target_feature = "avx2")))]
pub fn element_max(vec1: &mut [f32], val: f32) {
    for elem in vec1 {
        *elem = elem.max(val);
    }
}
#[cfg(all(not(target_feature = "sse3"),not(target_feature = "avx2")))]
pub fn element_add(vec1: &mut [f32], val: f32) {
    for elem in vec1 {
        *elem += val;
    }
}
#[cfg(all(not(target_feature = "sse3"),not(target_feature = "avx2")))]
pub fn element_multiply(vec1: &mut [f32], val: f32) {
    for elem in vec1 {
        *elem *= val;
    }
}

// 二つの一次元配列を与えて min, max, add（sum）等を行い、第一要素の配列に代入する
#[cfg(all(not(target_feature = "avx2"),target_feature = "sse3"))]
pub unsafe fn vector_op(vec1:&mut [f32],vec2:& [f32],op: SimdOpSse3){
    assert_eq!(vec1.len(),vec2.len());
    let mut i = 0;
    while i + 3 < vec1.len() {
        let chunk1 = _mm_loadu_ps(vec1.as_ptr().add(i));
        let chunk2 = _mm_loadu_ps(vec2.as_ptr().add(i));
        let result = op(chunk1, chunk2);
        _mm_storeu_ps(vec1.as_mut_ptr().add(i), result);
        i += 4;
    }
    //残りは別で処理
}


#[cfg(all(not(target_feature = "avx2"),target_feature = "sse3"))]
pub unsafe fn vector_max(vec1:&mut [f32],vec2:& [f32]){
    vector_op(vec1,vec2,_mm_max_ps);
    // ベクタのサイズが4の倍数でない場合の処理
    if vec1.len()%4 != 0{
        let i = vec1.len() - vec1.len()%4;
        for j in i..vec1.len() {
            vec1[j] = vec1[j].max(vec2[j]);
        }
    }
}


#[cfg(all(not(target_feature = "avx2"),target_feature = "sse3"))]
pub unsafe fn vector_min(vec1:&mut [f32],vec2:& [f32]){
    vector_op(vec1,vec2,_mm_min_ps);
    // ベクタのサイズが4の倍数でない場合の処理
    if vec1.len()%4 != 0{
        let i = vec1.len() - vec1.len()%4;
        for j in i..vec1.len() {
            vec1[j] = vec1[j].min(vec2[j]);
        }
    }
}

#[cfg(all(not(target_feature = "avx2"),target_feature = "sse3"))]
pub unsafe fn vector_add(vec1:&mut [f32],vec2:& [f32]){
    vector_op(vec1,vec2,_mm_add_ps);
    // ベクタのサイズが4の倍数でない場合の処理
    if vec1.len()%4 != 0{
        let i = vec1.len() - vec1.len()%4;
        for j in i..vec1.len() {
            vec1[j] = vec1[j] + vec2[j];
        }
    }
}

#[cfg(all(not(target_feature = "avx2"),target_feature = "sse3"))]
pub unsafe fn vector_multiply(vec1:&mut [f32],vec2:& [f32]){
    vector_op(vec1,vec2,_mm_mul_ps);
    // ベクタのサイズが4の倍数でない場合の処理
    if vec1.len()%4 != 0{
        let i = vec1.len() - vec1.len()%4;
        for j in i..vec1.len() {
            vec1[j] = vec1[j]*vec2[j];
        }
    }
}


#[cfg(target_feature = "avx2")]
pub unsafe fn vector_op(vec1:&mut [f32],vec2:&[f32],op:SimdOpAvx2){
    assert_eq!(vec1.len(),vec2.len());
    let mut i = 0;
    while i + 7 < vec1.len() {
        let chunk1 = _mm256_loadu_ps(vec1.as_ptr().add(i));
        let chunk2 = _mm256_loadu_ps(vec2.as_ptr().add(i));
        let result = op(chunk1,chunk2);
        _mm256_storeu_ps(vec1.as_mut_ptr().add(i), result);
        i += 8;
    }
}


#[cfg(target_feature = "avx2")]
pub unsafe fn vector_max(vec1:&mut [f32],vec2:&[f32]){
    vector_op(vec1,vec2,_mm256_max_ps);
    
    // ベクタのサイズが8の倍数でない場合の処理
    if vec1.len()%8 != 0{
        let i = vec1.len() - vec1.len()%8;
        for j in i..vec1.len() {
            vec1[j] = vec1[j].max(vec2[j]);
        }
    }
}

#[cfg(target_feature = "avx2")]
pub unsafe fn vector_min(vec1:&mut [f32],vec2:&[f32]){
    vector_op(vec1,vec2,_mm256_min_ps);
    
    // ベクタのサイズが8の倍数でない場合の処理
    if vec1.len()%8 != 0{
        let i = vec1.len() - vec1.len()%8;
        for j in i..vec1.len() {
            vec1[j] = vec1[j].min(vec2[j]);
        }
    }
}

#[cfg(target_feature = "avx2")]
pub unsafe fn vector_add(vec1:&mut [f32],vec2:&[f32]){
    vector_op(vec1,vec2,_mm256_add_ps);
    
    // ベクタのサイズが8の倍数でない場合の処理
    if vec1.len()%8 != 0{
        let i = vec1.len() - vec1.len()%8;
        for j in i..vec1.len() {
            vec1[j] = vec1[j]+vec2[j];
        }
    }
}


#[cfg(target_feature = "avx2")]
pub unsafe fn vector_multiply(vec1:&mut [f32],vec2:&[f32]){
    vector_op(vec1,vec2,_mm256_mul_ps);
    
    // ベクタのサイズが8の倍数でない場合の処理
    if vec1.len()%8 != 0{
        let i = vec1.len() - vec1.len()%8;
        for j in i..vec1.len() {
            vec1[j] = vec1[j]*vec2[j];
        }
    }
}


#[cfg(all(not(target_feature = "sse3"),not(target_feature = "avx2")))]
pub fn vector_max(vec1:&mut [f32],vec2:&[f32]){
    for (elem1,elem2) in vec1.iter_mut().zip(vec2.iter()) {
        *elem1 = (*elem1).max(*elem2);
    }
}

#[cfg(all(not(target_feature = "sse3"),not(target_feature = "avx2")))]
pub fn vector_min(vec1:&mut [f32],vec2:&[f32]){
    for (elem1,elem2) in vec1.iter_mut().zip(vec2.iter()) {
        *elem1 = (*elem1).min(*elem2);
    }
}

#[cfg(all(not(target_feature = "sse3"),not(target_feature = "avx2")))]
pub fn vector_add(vec1:&mut [f32],vec2:&[f32]){
    for (elem1,elem2) in vec1.iter_mut().zip(vec2.iter()) {
        *elem1 += *elem2;
    }
}

#[cfg(all(not(target_feature = "sse3"),not(target_feature = "avx2")))]
pub fn vector_multiply(vec1:&mut [f32],vec2:&[f32]){
    for (elem1,elem2) in vec1.iter_mut().zip(vec2.iter()) {
        *elem1 *= *elem2;
    }
}

//一次元配列の各要素を二乗する
#[cfg(all(not(target_feature = "avx2"),target_feature = "sse3"))]
pub unsafe fn vector_square(vec: &mut [f32]){
    let mut i = 0;
    while i + 3 < vec.len() {
        let chunk = _mm_loadu_ps(vec.as_ptr().add(i));
        let result = _mm_mul_ps(chunk, chunk);
        _mm_storeu_ps(vec.as_mut_ptr().add(i), result);
        i += 4;
    }
    // ベクタのサイズが4の倍数でない場合の処理
    for j in i..vec.len() {
        vec[j] = vec[j]*vec[j] ;
    }
}


#[cfg(target_feature = "avx2")]
pub unsafe fn vector_square(vec: &mut [f32]){
    let mut i = 0;
    while i + 7 < vec.len() {
        let chunk = _mm256_loadu_ps(vec.as_ptr().add(i));
        let result = _mm256_mul_ps(chunk, chunk);
        _mm256_storeu_ps(vec.as_mut_ptr().add(i), result);
        i += 8;
    }
    // ベクタのサイズが8の倍数でない場合の処理
    for j in i..vec.len() {
        vec[j] = vec[j]*vec[j] ;
    }
}

#[cfg(all(not(target_feature = "sse3"),not(target_feature = "avx2")))]
pub fn vector_square(vec: &mut[f32]){
    for j in 0..vec.len() {
        vec[j] = vec[j]*vec[j] ;
    }
}



//一次元配列の各要素のルートを取る
#[cfg(all(not(target_feature = "avx2"),target_feature = "sse3"))]
pub unsafe fn vector_sqrt(vec: &mut [f32]){
    let mut i = 0;
    while i + 3 < vec.len() {
        let chunk = _mm_loadu_ps(vec.as_ptr().add(i));
        let result = _mm_sqrt_ps(chunk);
        _mm_storeu_ps(vec.as_mut_ptr().add(i), result);
        i += 4;
    }
    // ベクタのサイズが4の倍数でない場合の処理
    for j in i..vec.len() {
        vec[j] = vec[j].sqrt();
    }
}


#[cfg(target_feature = "avx2")]
pub unsafe fn vector_sqrt(vec: &mut [f32]){
    let mut i = 0;
    while i + 7 < vec.len() {
        let chunk = _mm256_loadu_ps(vec.as_ptr().add(i));
        let result = _mm256_sqrt_ps(chunk);
        _mm256_storeu_ps(vec.as_mut_ptr().add(i), result);
        i += 8;
    }
    // ベクタのサイズが8の倍数でない場合の処理
    for j in i..vec.len() {
        vec[j] = vec[j].sqrt() ;
    }
}

#[cfg(all(not(target_feature = "sse3"),not(target_feature = "avx2")))]
pub fn vector_sqrt(vec: &mut[f32]){
    for j in 0..vec.len() {
        vec[j] = vec[j].sqrt();
    }
}


//一次元配列を 2 つ与えてユークリッド距離を計算する
pub fn calc_euclid_dist(vec1: &Vec<f32>,vec2: &Vec<f32>)->f32{
    let mut mvec1:Vec<f32> = vec1.clone();//破壊するので Clone
    unsafe{
        element_multiply(&mut mvec1,-1.0);
        vector_add(&mut mvec1,&vec2);
        vector_square(&mut mvec1);
    }
    let ret:f32 = mvec1.into_iter().sum(); //まあ多分ベクトル化してくれるのでは・・・
    return ret.sqrt();
}
#[derive(Debug)]
pub struct VectorStats{
    pub mean:f32,
    pub var:f32,
    pub count:usize
}

pub fn calc_stats(vec1:&Vec<f32>)->VectorStats{
    assert!(vec1.len() != 0);
    let count:usize = vec1.len();
    let mmean:f32 = vec1.iter().sum::<f32>()/(count as f32);
    let mut mvec1 = vec1.clone();
    unsafe{
        element_add(&mut mvec1,mmean*-1.0);
        vector_square(&mut mvec1);
    }
    let vvar:f32 = mvec1.into_iter().sum::<f32>()/(count as f32);
    return VectorStats{
        mean:mmean,var:vvar,count:count
    };
}

#[test]
fn matrix_test(){
    let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let b = vec![8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
    unsafe{
        let result = dot_product(&a, &b);
        println!("Dot Product: {}", result);
        println!("Euclid Distance: {}", calc_euclid_dist(&a,&b));
    }
    
}
