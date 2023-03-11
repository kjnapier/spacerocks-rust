// access the constants from constants.rs in this directory
use crate::constants::*;
use crate::statevector::StateVector;
use crate::observatory::Observatory;

use nalgebra::Vector3;

pub struct SpaceRock {
    pub name: String,
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub epoch: f64,
    pub frame: String,
    pub origin: String,
    pub mass: Option<f64>
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

    // pub fn observe(observatory: &Observatory) -> Detection {
    //     return None
    // }

    pub fn change_frame(&mut self, frame: &str) {
        if frame != self.frame {
            let inv = ROTATION_MATRICES[&self.frame].try_inverse().unwrap();
            let rot = ROTATION_MATRICES[frame] * inv;
            self.position = rot * self.position;
            self.velocity = rot * self.velocity;
            self.frame = frame.to_string();
        }
    }

    // fn change_frame(&mut self, frame: &str) {
    //     let (state, _) = spice::spkezr(&self.name, self.epoch, frame, "NONE", "SSB");
    //     self.position = Vector3::new(state[0] * KM_TO_AU, state[1] * KM_TO_AU, state[2] * KM_TO_AU);
    //     self.velocity = Vector3::new(state[3] * KM_TO_AU * SECONDS_PER_DAY, state[4] * KM_TO_AU * SECONDS_PER_DAY, state[5] * KM_TO_AU * SECONDS_PER_DAY);
    //     self.frame = frame.to_string();
    // }

    // fn change_origin(&mut self, origin: &str) {
    //     let (state, _) = spice::spkezr(&self.name, self.epoch, "J2000", "NONE", origin);
    //     self.position = Vector3::new(state[0] * KM_TO_AU, state[1] * KM_TO_AU, state[2] * KM_TO_AU);
    //     self.velocity = Vector3::new(state[3] * KM_TO_AU * SECONDS_PER_DAY, state[4] * KM_TO_AU * SECONDS_PER_DAY, state[5] * KM_TO_AU * SECONDS_PER_DAY);
    //     self.origin = origin.to_string();
    // }

    fn r_squared(&self) -> f64 {
        self.position.dot(&self.position)
    }

    fn r(&self) -> f64 {
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
