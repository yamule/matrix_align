

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
            // 出力を a[1],a[1],b[1],b[1] とし、最初の要素を取り出す
            sum += _mm_cvtss_f32(sum_chunk) + 
            _mm_cvtss_f32(_mm_shuffle_ps(sum_chunk, sum_chunk, 0x55));
        }
    }
    sum
}


#[cfg(target_feature = "avx2")]
fn dot_product(a: &[f32], b: &[f32]) -> f32 {
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
fn dot_product(a: &[f32], b: &[f32]) -> f32 {
    assert_eq!(a.len(), b.len());
    eprintln!("running native");
    return a.iter().zip(b.iter()).map(|(x,y)|x*y).sum();
}


#[test]
fn test(){
    let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let b = vec![8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
    let result = dot_product(&a, &b);
    println!("Dot Product: {}", result);
}
