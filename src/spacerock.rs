// access the constants from constants.rs in this directory
use crate::constants::*;
use crate::statevector::StateVector;
use crate::observatory::Observatory;
use crate::correct_for_ltt::correct_for_ltt;

use nalgebra::Vector3;

pub struct SpaceRock {
    pub name: String,
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub epoch: f64,
    pub frame: String,
    pub origin: String,
    pub mass: Option<f64>,
    // pub radius: Option<f64>,
    // pub H: Option<f64>,
    // pub G: Option<f64>,
    // pub mag: Option<f64>,
    // pub orbit: Option<KeplerOrbit>,
    // pub a: Option<f64>,
    // pub e: Option<f64>,
    // pub inc: Option<f64>,
    // pub arg: Option<f64>,
    // pub node: Option<f64>,
    // pub M: Option<f64>,
    // pub f: Option<f64>,
    // pub E: Option<f64>,
    // pub n: Option<f64>,
    // pub P: Option<f64>,
    // pub T: Option<f64>,
    // pub q: Option<f64>,
}

#[allow(dead_code)]
impl SpaceRock {

    // Instantiation Methods

    pub fn from_spice(name: &str, epoch: f64) -> Self {
        let et = spice::str2et(&format!("JD{epoch} UTC", epoch=epoch));
        let (state, _) = spice::spkezr(name, et, "J2000", "NONE", "SSB");
        SpaceRock {
            name: name.to_string(), 
            position: Vector3::new(state[0], state[1], state[2]) * KM_TO_AU,
            velocity: Vector3::new(state[3], state[4], state[5]) * KM_TO_AU * SECONDS_PER_DAY,
            epoch: epoch,
            frame: "J2000".to_string(),
            origin: "SSB".to_string(),
            mass: None
        }
    }

    pub fn from_xyz(name: &str, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64, epoch: f64) -> Self {
        SpaceRock {
            name: name.to_string(),
            position: Vector3::new(x, y, z),
            velocity: Vector3::new(vx, vy, vz),
            epoch: epoch,
            frame: "J2000".to_string(),
            origin: "SSB".to_string(),
            mass: None
        }
    }

    pub fn from_state(name: &str, state: StateVector, epoch: f64) -> Self {
        SpaceRock {
            name: name.to_string(),
            position: state.position,
            velocity: state.velocity,
            epoch: epoch,
            frame: "J2000".to_string(),
            origin: "SSB".to_string(),
            mass: None
        }
    }

    pub fn observe(&mut self, observer: &SpaceRock) -> [f64; 2] {
        self.change_frame("J2000");
        let corrected_rock = correct_for_ltt(&self, observer);
        let ra = corrected_rock.position.y.atan2(corrected_rock.position.x);
        let dec = (corrected_rock.position.z / corrected_rock.position.norm()).asin();
        return [ra, dec];
    }

    pub fn change_frame(&mut self, frame: &str) {
        if frame != self.frame {
            let inv = ROTATION_MATRICES[&self.frame].try_inverse().unwrap();
            let rot = ROTATION_MATRICES[frame] * inv;
            self.position = rot * self.position;
            self.velocity = rot * self.velocity;
            self.frame = frame.to_string();
        }
    }

    fn r_squared(&self) -> f64 {
        self.position.dot(&self.position)
    }

    pub fn r(&self) -> f64 {
        self.position.norm()
    }

    fn v_squared(&self) -> f64 {
        self.velocity.dot(&self.velocity)
    }

    fn v(&self) -> f64 {
        self.velocity.norm()
    }

}

#[allow(dead_code)]
fn separation(body1: &SpaceRock, body2: &SpaceRock) -> f64 {
    let d_pos = body1.position - body2.position;
    return d_pos.norm();
}
