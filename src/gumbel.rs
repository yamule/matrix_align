use rand::distributions::Open01;
use rand::Rng;
use std::f64::consts::PI;

// translated from hmmer source code by ChatGPT
pub mod gumbel {
    use super::*;

    /// Returns the probability density at `x`, `P(S=x)`.
    pub fn pdf(x: f64, mu: f64, lambda: f64) -> f64 {
        let y = lambda * (x - mu);
        if y < -6.5 || y > 710.0 {
            0.0
        } else {
            lambda * ((-y - (-y).exp()).exp())
        }
    }

    /// Returns the log of the pdf at `x`, `log P(S=x)`.
    pub fn logpdf(x: f64, mu: f64, lambda: f64) -> f64 {
        let y = lambda * (x - mu);
        if y < -708.0 || y > 708.0 {
            f64::NEG_INFINITY
        } else {
            lambda.ln() - y - (-y).exp()
        }
    }

    /// Returns the cumulative distribution at `x`, `P(S ≤ x)`.
    pub fn cdf(x: f64, mu: f64, lambda: f64) -> f64 {
        let y = lambda * (x - mu);
        if y < -6.5 {
            0.0
        } else if y > 36.0 {
            1.0
        } else {
            (-(-y).exp()).exp()
        }
    }

    /// Returns the log of the cdf at `x`, `log P(S ≤ x)`.
    pub fn logcdf(x: f64, mu: f64, lambda: f64) -> f64 {
        let y = lambda * (x - mu);
        if y < -708.0 {
            f64::NEG_INFINITY
        } else if y > 708.0 {
            0.0
        } else {
            -(-y).exp()
        }
    }

    /// Returns right tail mass above `x`, `P(S > x)`.
    pub fn surv(x: f64, mu: f64, lambda: f64) -> f64 {
        let y = lambda * (x - mu);
        let ey = -(-y).exp();

        if ey.abs() < 1e-15 {
            -ey
        } else {
            1.0 - ey.exp()
        }
    }

    /// Returns log survival at `x`, `log P(S > x)`.
    pub fn logsurv(x: f64, mu: f64, lambda: f64) -> f64 {
        let y = lambda * (x - mu);
        let ey = -(-y).exp();

        if ey.abs() < 1e-15 {
            -y
        } else if ey.exp().abs() < 1e-15 {
            -ey.exp()
        } else {
            (1.0 - ey.exp()).ln()
        }
    }

    /// Calculates the inverse CDF for a Gumbel distribution.
    pub fn invcdf(p: f64, mu: f64, lambda: f64) -> f64 {
        mu - ((-p.ln()).ln()) / lambda
    }

    /// Calculates the score at which the right tail's mass is `p`.
    pub fn invsurv(p: f64, mu: f64, lambda: f64) -> f64 {
        let log_part = if p < 1e-15 {
            (p.powf(p) - 1.0) / p
        } else {
            (- (1.0 - p).ln()).ln()
        };

        mu - log_part / lambda
    }

    /// Returns a Gumbel-distributed random sample `x`.
    pub fn sample<R: Rng + ?Sized>(rng: &mut R, mu: f64, lambda: f64) -> f64 {
        let Open01(p) = rng.sample(Open01);
        invcdf(p, mu, lambda)
    }

    /// Estimates `mu` and `lambda` from complete data.
    pub fn fit_complete(x: &[f64]) -> Result<(f64, f64), String> {
        let n = x.len();
        if n <= 1 {
            return Err("n must be greater than 1".to_string());
        }
        let n_f64 = n as f64;

        // Compute mean and variance
        let mean = x.iter().sum::<f64>() / n_f64;
        let variance = x
            .iter()
            .map(|xi| (xi - mean).powi(2))
            .sum::<f64>() / n_f64;

        // Initial guess for lambda
        let mut lambda = PI / (6.0 * variance).sqrt();

        let tol = 1e-5;
        let mut fx;
        let mut dfx;

        // Newton-Raphson iterations
        let mut i = 0;
        for _ in 0..100 {
            let (f, df) = lawless416(x, lambda);
            fx = f;
            dfx = df;
            if fx.abs() < tol {
                break;
            }
            lambda = lambda - fx / dfx;
            if lambda <= 0.0 {
                lambda = 0.001;
            }
            i += 1;
        }

        // If Newton-Raphson failed, use bisection method
        if i == 100 {
            // Bisection method
            let mut left = 0.0;
            let mut right = PI / (6.0 * variance).sqrt();
            let (mut fx, _) = lawless416(x, lambda);
            while fx > 0.0 {
                right *= 2.0;
                if right > 1000.0 {
                    return Err("Failed to bracket root in fit_complete".to_string());
                }
                let (f, _) = lawless416(x, right);
                fx = f;
            }
            for _ in 0..100 {
                let mid = (left + right) / 2.0;
                let (f, _) = lawless416(x, mid);
                fx = f;
                if fx.abs() < tol {
                    lambda = mid;
                    break;
                }
                if fx > 0.0 {
                    left = mid;
                } else {
                    right = mid;
                }
            }
            if fx.abs() >= tol {
                return Err("Even bisection search failed in fit_complete".to_string());
            }
        }

        // Compute mu
        let esum = x.iter().map(|&xi| (-lambda * xi).exp()).sum::<f64>();
        let mu = - (esum / n_f64).ln() / lambda;

        Ok((mu, lambda))
    }

    /// Estimates `mu` from complete data, given `lambda`.
    pub fn fit_complete_loc(x: &[f64], lambda: f64) -> Result<f64, String> {
        let n = x.len();
        if n <= 1 {
            return Err("n must be greater than 1".to_string());
        }
        let n_f64 = n as f64;

        let esum = x.iter().map(|&xi| (-lambda * xi).exp()).sum::<f64>();
        let mu = - (esum / n_f64).ln() / lambda;

        Ok(mu)
    }

    // Helper function for lawless416
    fn lawless416(x: &[f64], lambda: f64) -> (f64, f64) {
        let n = x.len() as f64;
        let mut esum = 0.0;
        let mut xesum = 0.0;
        let mut xxesum = 0.0;
        let mut xsum = 0.0;

        for &xi in x.iter() {
            let exp_term = (-lambda * xi).exp();
            xsum += xi;
            xesum += xi * exp_term;
            xxesum += xi * xi * exp_term;
            esum += exp_term;
        }

        let f = (1.0 / lambda) - (xsum / n) + (xesum / esum);
        let df = ((xesum / esum).powi(2)) - (xxesum / esum) - (1.0 / (lambda * lambda));
        (f, df)
    }
}

fn main() {
    // Example usage
    let mut rng = rand::thread_rng();
    let mu = -20.0;
    let lambda = 0.4;

    // Generate samples
    let samples: Vec<f64> = (0..10000)
        .map(|_| gumbel::sample(&mut rng, mu, lambda))
        .collect();

    // Fit parameters
    match gumbel::fit_complete(&samples) {
        Ok((est_mu, est_lambda)) => {
            println!("Estimated mu: {}, Estimated lambda: {}", est_mu, est_lambda);
        }
        Err(e) => println!("Error in fitting: {}", e),
    }
}
