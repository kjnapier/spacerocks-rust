
pub fn calc_E_from_M(e: f64, M: f64) -> f64 {

    let mut E;

    if e == 0.0 {
        E = M;
        return E;
    }

    if M == 0.0 {
        E = 0.0;
        return E;
    }
    
    if e < 1.0 {

        let k = 0.85;
        

        // Define initial estimate
        if M.sin() < 0.0 {
            E = M - k * e;
        }
        else {
            E = M + k * e;
        }

        // Perform Newton-Raphson estimate
        for _ in 0..10 {

            // Compute f(E), f'(E), f''(E) and f'''(E), avoiding recomputation of sine and cosine.
            let esinE = e * E.sin();
            let ecosE = e * E.cos();
            let f_E = E - esinE - M;
            let fP_E = 1.0 - ecosE;
            let fPP_E = esinE;
            let fPPP_E = ecosE;

            let delta_i1 = -f_E / fP_E;
            let delta_i2 = -f_E / (fP_E + 0.5 * delta_i1 * fPP_E);
            let delta_i3 = -f_E / (fP_E + 0.5 * delta_i2 * fPP_E + 1.0/6.0 * fPPP_E * delta_i2 * delta_i2);
            
            // Update E
            E += delta_i3;
            if delta_i3.abs() < 1.0e-16 {
                break;
        }
      }
    }

    else {

        E = M / M.abs() * (2.0 * M.abs() / e + 1.8).ln();
  		let mut F = E - e * E.sinh() + M;
  		for _ in 1..100 {

  			E = E - F / (1.0 - e * E.cosh());
  			F = E - e * E.sinh() + M;
  			if F.abs() < 1.0e-16 {
  				break;
        }
      }
    }
    return E;
}