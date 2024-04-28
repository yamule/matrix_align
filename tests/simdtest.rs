#[cfg(test)]
mod simdtest{
    use std::{collections::HashMap};
    use matrix_align::matrix_process::*;
    use rand::{Rng, SeedableRng};
    use rand::rngs::StdRng;
    const ERROR_TOLERANCE:f32 = 0.00001;
    pub fn is_close_vec(vec1:&Vec<f32>,vec2:&Vec<f32>,tol:f32){
        assert_eq!(vec1.len(),vec2.len());
        for ii in 0..vec1.len(){
            assert!(
                (vec1[ii]-vec2[ii]).abs() <  tol, "Large difference in index {}. {} vs {}\n{:?}\n{:?}",ii,vec1[ii],vec2[ii],vec1,vec2
            );
        }
    }
    
    pub fn is_close(val1:f32,val2:f32,tol:f32){
        assert!((val1-val2).abs() < tol,"Large difference was foun {} vs {}.",val1,val2);
    }

    #[allow(unused_mut)]
    #[test]
    pub fn test_simd_all(){
        unsafe {

            let mut rng = StdRng::seed_from_u64(123);
                    
            check_simd();
            
            let mut vec1_orig:Vec<f32> = vec![1.0,2.4,3.2,2.4,3.2,2.4,3.2,2.4,3.2];
            let mut vec2_orig:Vec<f32> = vec![5.3,1.2,1.1,2.4,3.2,2.4,3.2,2.4,3.2];
            let mut val = 32.0;

            let mut vec1 = vec1_orig.clone();
            let vec2 = vec2_orig.clone();
            let mut vec1_native = vec1_orig.clone();
            let vec2_native = vec2_orig.clone();
            
            let a = dot_product(&vec1,&vec2);
            let b = dot_product_native(&vec1_native,&vec2_native);
            is_close(a,b,0.00001);
            
            let mut vec1 = vec1_orig.clone();
            let mut vec1_native = vec1_orig.clone();
            element_min(&mut vec1, val);
            element_min_native(&mut vec1_native, val);
            is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);

            let mut vec1 = vec1_orig.clone();
            let mut vec1_native = vec1_orig.clone();
            element_max(&mut vec1, val);
            element_max_native(&mut vec1_native, val);
            is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
            
            let mut vec1 = vec1_orig.clone();
            let mut vec1_native = vec1_orig.clone();
            element_add(&mut vec1, val);
            element_add_native(&mut vec1_native, val);
            is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
            
            let mut vec1 = vec1_orig.clone();
            let mut vec1_native = vec1_orig.clone();
            element_multiply(&mut vec1,val);
            element_multiply_native(&mut vec1_native,val);
            is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);

            let mut vec1 = vec1_orig.clone();
            let vec2 = vec2_orig.clone();
            let mut vec1_native = vec1_orig.clone();
            let vec2_native = vec2_orig.clone();
            vector_max(&mut vec1,&vec2);
            vector_max_native(&mut vec1_native,&vec2_native);
            is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
            
            let mut vec1 = vec1_orig.clone();
            let vec2 = vec2_orig.clone();
            let mut vec1_native = vec1_orig.clone();
            let vec2_native = vec2_orig.clone();
            vector_min(&mut vec1,&vec2);
            vector_min_native(&mut vec1_native,&vec2_native);
            is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);

            let mut vec1 = vec1_orig.clone();
            let vec2 = vec2_orig.clone();
            let mut vec1_native = vec1_orig.clone();
            let vec2_native = vec2_orig.clone();
            vector_add(&mut vec1,&vec2);
            vector_add_native(&mut vec1_native,&vec2_native);
            is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
            
            let mut vec1 = vec1_orig.clone();
            let vec2 = vec2_orig.clone();
            let mut vec1_native = vec1_orig.clone();
            let vec2_native = vec2_orig.clone();
            vector_multiply(&mut vec1,&vec2);
            vector_multiply_native(&mut vec1_native,&vec2_native);
            is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
            
            let mut vec1 = vec1_orig.clone();
            let vec2 = vec2_orig.clone();
            let mut vec1_native = vec1_orig.clone();
            let vec2_native = vec2_orig.clone();
            vector_square(&mut vec1);
            vector_square_native(&mut vec1_native);
            is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
            
            // プラスでないと拙いのでコピーしない
            vector_sqrt(&mut vec1);
            vector_sqrt_native(&mut vec1_native);
            is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);

            for ii in 1..270{
;

                let mut vec1_orig:Vec<f32> = vec![];
                let mut vec2_orig:Vec<f32> = vec![];
                for jj in 0..ii{
                    let a:f32 = rng.gen_range(-2.0..2.0);
                    let b:f32 = rng.gen_range(-2.0..2.0);
                    vec1_orig.push(a);
                    vec2_orig.push(b);
                }
                //println!("{:?}",vec1_orig);
                //println!("{:?}",vec2_orig);
                let mut val = 32.0;
    
                let mut vec1 = vec1_orig.clone();
                let vec2 = vec2_orig.clone();
                let mut vec1_native = vec1_orig.clone();
                let vec2_native = vec2_orig.clone();
                
                let a = dot_product(&vec1,&vec2);
                let b = dot_product_native(&vec1_native,&vec2_native);
                is_close(a,b,0.00001*(ii as f32));
                
                let mut vec1 = vec1_orig.clone();
                let mut vec1_native = vec1_orig.clone();
                element_min(&mut vec1, val);
                element_min_native(&mut vec1_native, val);
                is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
    
                let mut vec1 = vec1_orig.clone();
                let mut vec1_native = vec1_orig.clone();
                element_max(&mut vec1, val);
                element_max_native(&mut vec1_native, val);
                is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
                
                let mut vec1 = vec1_orig.clone();
                let mut vec1_native = vec1_orig.clone();
                element_add(&mut vec1, val);
                element_add_native(&mut vec1_native, val);
                is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
                
                let mut vec1 = vec1_orig.clone();
                let mut vec1_native = vec1_orig.clone();
                element_multiply(&mut vec1,val);
                element_multiply_native(&mut vec1_native,val);
                is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
    
                let mut vec1 = vec1_orig.clone();
                let vec2 = vec2_orig.clone();
                let mut vec1_native = vec1_orig.clone();
                let vec2_native = vec2_orig.clone();
                vector_max(&mut vec1,&vec2);
                vector_max_native(&mut vec1_native,&vec2_native);
                is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
                
                let mut vec1 = vec1_orig.clone();
                let vec2 = vec2_orig.clone();
                let mut vec1_native = vec1_orig.clone();
                let vec2_native = vec2_orig.clone();
                vector_min(&mut vec1,&vec2);
                vector_min_native(&mut vec1_native,&vec2_native);
                is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
    
                let mut vec1 = vec1_orig.clone();
                let vec2 = vec2_orig.clone();
                let mut vec1_native = vec1_orig.clone();
                let vec2_native = vec2_orig.clone();
                vector_add(&mut vec1,&vec2);
                vector_add_native(&mut vec1_native,&vec2_native);
                is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
                
                let mut vec1 = vec1_orig.clone();
                let vec2 = vec2_orig.clone();
                let mut vec1_native = vec1_orig.clone();
                let vec2_native = vec2_orig.clone();
                vector_multiply(&mut vec1,&vec2);
                vector_multiply_native(&mut vec1_native,&vec2_native);
                is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
                
                let mut vec1 = vec1_orig.clone();
                let vec2 = vec2_orig.clone();
                let mut vec1_native = vec1_orig.clone();
                let vec2_native = vec2_orig.clone();
                vector_square(&mut vec1);
                vector_square_native(&mut vec1_native);
                is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);
                
                // プラスでないと拙いのでコピーしない
                vector_sqrt(&mut vec1);
                vector_sqrt_native(&mut vec1_native);
                is_close_vec(&vec1, &vec1_native,ERROR_TOLERANCE);

            }
        }       
    }   
}