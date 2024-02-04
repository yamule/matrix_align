

/*
use std::simd::f32x4;

fn dot_product(a: &[f32], b: &[f32]) -> f32 {
    assert_eq!(a.len(), b.len());
    assert!(a.len() % 4 == 0);

    let mut sum = f32x4::splat(0.0);
    for i in (0..a.len()).step_by(4) {
        let a_chunk = f32x4::from_slice(&a[i..i + 4]);
        let b_chunk = f32x4::from_slice(&b[i..i + 4]);
        sum += a_chunk * b_chunk;
    }

    sum.reduce_sum()
}
*/


use std::arch::x86_64::*;
#[cfg(not(target_feature = "avx2"))]
#[cfg(target_feature = "sse3")]
pub fn dot_product(a: &[f32], b: &[f32]) -> f32 {
    assert_eq!(a.len(), b.len());
    assert!(a.len() % 4 == 0);

    eprintln!("running sse3");
    let mut sum = 0.0;
    for i in (0..a.len()).step_by(4) {
        unsafe {
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
    }
    sum
}


#[cfg(target_feature = "avx2")]
pub fn dot_product(a: &[f32], b: &[f32]) -> f32 {
    assert_eq!(a.len(), b.len());
    assert!(a.len() % 8 == 0);

    eprintln!("running avx2");
    let mut sum = 0.0;
    unsafe {
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
    }
    sum
}

#[cfg(not(target_feature = "sse3"))]
#[cfg(not(target_feature = "avx2"))]
pub fn dot_product(a: &[f32], b: &[f32]) -> f32 {
    assert_eq!(a.len(), b.len());
    //eprintln!("running native");
    return a.iter().zip(b.iter()).map(|(x,y)|x*y).sum();
}


#[cfg(not(target_feature = "avx2"))]
#[cfg(target_feature = "sse3")]
unsafe fn multiply_elements(vec: &mut [f32], factor: f32) {
    let factor_vec = _mm_set1_ps(factor);

    let mut i = 0;
    while i + 3 < vec.len() {
        let chunk = _mm_loadu_ps(vec.as_ptr().add(i));
        let result = _mm_mul_ps(chunk, factor_vec);
        _mm_storeu_ps(vec.as_mut_ptr().add(i), result);
        i += 4;
    }

    // ベクタのサイズが4の倍数でない場合の処理
    for j in i..vec.len() {
        vec[j] *= factor;
    }
}


#[cfg(target_feature = "avx2")]
unsafe fn multiply_elements(vec: &mut [f32], factor: f32) {
    let factor_vec = _mm256_set1_ps(factor);

    let mut i = 0;
    while i + 7 < vec.len() {
        let chunk = _mm256_loadu_ps(vec.as_ptr().add(i));
        let result = _mm256_mul_ps(chunk, factor_vec);
        _mm256_storeu_ps(vec.as_mut_ptr().add(i), result);
        i += 8;
    }

    // ベクタのサイズが8の倍数でない場合の処理
    for j in i..vec.len() {
        vec[j] *= factor;
    }
}


#[cfg(not(target_feature = "sse3"))]
#[cfg(not(target_feature = "avx2"))]
pub fn multiply_elements(vec: &mut [f32], factor: f32) {
    for elem in vec {
        *elem *= factor;
    }
}

#[test]
fn test(){
    let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let b = vec![8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
    let result = dot_product(&a, &b);
    println!("Dot Product: {}", result);
}
