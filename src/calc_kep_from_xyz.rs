use crate::constants::MU_BARY;
use crate::statevector::StateVector;
use crate::keplerorbit::KeplerOrbit;
use nalgebra::Vector3;



pub fn calc_kep_from_xyz(state: StateVector) -> KeplerOrbit {

    let EMIN = 1e-10;
    let IMIN = 1e-10;

    let (a, e, inc, mut arg, mut node, mut f);

    let r = state.position.norm();
    let vsq = state.velocity.dot(&state.velocity);

    let hvec = state.position.cross(&state.velocity);
    let evec = state.velocity.cross(&hvec) / MU_BARY - state.position / r;  

    let nvec = Vector3::new(-hvec.y, hvec.x, 0.0);
    let n = nvec.norm();

    a = 1.0 / (2.0 / r - vsq / MU_BARY);
    e = evec.norm();
    inc = (hvec.z / hvec.norm()).acos();

    if inc == 0.0 {
      node = 0.0;
    }
    else {
      node = (nvec.x / n).acos();
    }
    if nvec.y < 0.0 {
      node = 2.0 * std::f64::consts::PI - node;
    }
  
    // Compute the argument of pericenter (arg)
    if e < EMIN {
      arg = 0.0;
    }
    else if inc < IMIN || inc > std::f64::consts::PI - IMIN {
      // let mut theta = (state.position.x / r).acos();
      // if state.position.y < 0.0 {
      //   theta = 2.0 * std::f64::consts::PI - theta;
      // }
      let mut varpi = (evec.x / e).acos();
      if evec.y < 0.0 {
        varpi = 2.0 * std::f64::consts::PI - varpi;
      }
      if inc < std::f64::consts::PI / 2.0 {
        arg = varpi - node;
      }
      else {
        arg = node - varpi;
      }
    }
    else {
      arg = (nvec.dot(&evec) / (n * e)).acos();
      if evec.z < 0.0 {
        arg = 2.0 * std::f64::consts::PI - arg;
      }
    }
      
    // Compute the true anomaly (f)
    if e < 1.0 {
      if inc < IMIN || inc > std::f64::consts::PI - IMIN {
        // Handling the near-planar case
        if e > EMIN {
          // Near-planar, elliptical
          let theta = (state.position.x / r).acos();
          let varpi = (evec.x / e).acos();
          if inc < std::f64::consts::PI/2.0 {
            f = theta - varpi;
          }
          else {
            f = varpi - theta;
          }
        }
        else {
          // Near-planar, near-circular
          f = (state.position.x / r).acos();
          if state.velocity.x > 0.0 {
            f = 2.0 * std::f64::consts::PI - f;
          }
        }
      } 
      else {
        // Handling the non-planar case
        if e > EMIN {
          // Non-planar, elliptical
          let edotr = evec.dot(&state.position);
          let rdotv = state.position.dot(&state.velocity);
          f = (edotr / (e * r)).acos();
          if rdotv < 0.0 {
            f = 2.0 * std::f64::consts::PI - f;
          }
        }
        else {
          // Non-planar, circular
          f = (nvec.dot(&state.position) / (n * r)).acos();
          if state.position.z < 0.0 {
            f = 2.0 * std::f64::consts::PI - f;
          }
        }
      }
    }
    else {
      let mut argument = (1.0 - r / a) / e;
      if (argument - 1.0).abs() < 1e-10 {
        argument = 1.0;
      }
      let mut E = argument.acosh();
      let rdotv = state.position.dot(&state.velocity);
      if rdotv < 0.0 {
       E *= -1.0;
      }
      f = 2.0 * ((e + 1.0).sqrt() * (E / 2.0).tanh()).atan2((e - 1.0).sqrt());
    }
  
    let kep = KeplerOrbit::new(a, e, inc, arg, node, f);
    return kep;
      
}