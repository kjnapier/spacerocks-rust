use crate::constants::MU_BARY;
use crate::statevector::StateVector;
use crate::calc_kep_from_xyz;

pub struct KeplerOrbit {
    pub a: f64,
    pub e: f64,
    pub inc: f64,
    pub arg: f64,
    pub node: f64,
    pub f: f64,
}

impl KeplerOrbit {
    pub fn from_xyz(state: StateVector) -> Self {
        let kep = calc_kep_from_xyz::calc_kep_from_xyz(state);
        KeplerOrbit {
            a: kep.a,
            e: kep.e,
            inc: kep.inc,
            arg: kep.arg,
            node: kep.node,
            f: kep.f,
        }
    }

    pub fn new(a: f64, e: f64, inc: f64, arg: f64, node: f64, f: f64) -> Self {
        KeplerOrbit {
            a: a,
            e: e,
            inc: inc,
            arg: arg,
            node: node,
            f: f,
        }
    }
}
