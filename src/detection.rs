use nalgebra::Vector3;

#[allow(dead_code)]
#[derive(Debug)]
pub struct Detection {
    pub ra: f64,
    pub dec: f64,
    pub ra_rate: f64,
    pub dec_rate: f64,
    pub epoch: f64,
    pub objid: String,
    pub obscode: String,
    pub observer: Vector3<f64>,
    pub pointing_vector: Vector3<f64>
}

impl Detection {

    pub fn new(ra: f64, dec: f64, ra_rate: f64, dec_rate: f64, epoch: f64, objid: String, obscode: String, observer: Vector3<f64>) -> Self {

        let pointing_vector = compute_pointing_vector(ra, dec);
        Detection {
            ra: ra,
            dec: dec,
            ra_rate: ra_rate,
            dec_rate: dec_rate,
            epoch: epoch,
            objid: objid,
            obscode: obscode,
            observer: observer,
            pointing_vector: pointing_vector
        }
    }

    pub fn generate_orbit(&self, r: f64, r_rate: f64) -> SpaceRock {
        let rho = self.calculate_rho(r);
        let rho_rate = self.calculate_rho_rate(r, r_rate);

        let x = self.observer.position.x + rho * self.pointing_vector.x;
        let y = self.observer.position.y + rho * self.pointing_vector.y;
        let z = self.observer.position.z + rho * self.pointing_vector.z;
        let vx = self.observer.velocity.x + rho_rate * self.pointing_vector.x + rho * self.pointing_vector_rate.x;
        let vy = self.observer.velocity.y + rho_rate * self.pointing_vector.y + rho * self.pointing_vector_rate.y;
        let vz = self.observer.velocity.z + rho_rate * self.pointing_vector.z + rho * self.pointing_vector_rate.z;

        let ltt = rho / speed_of_light.value;

        let rock = SpaceRock(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, epoch=self.epoch.utc.jd - ltt, frame='J2000', origin='ssb');
        return rock;
    }

    fn calculate_rho(&self, r: f64) -> f64 {
        let A = self.observer.x * self.pointing_vector.x.cos();
        let B = (self.observer.x.powi(2) * (self.pointing_vector.x.cos().powi(2) - 1.0) + r.powi(2)).sqrt();
        return A + B;
    }

    fn calculate_rho_rate(&self, r: f64, r_rate: f64) -> f64 {
        let A = self.observer.velocity.x * self.pointing_vector.x.cos();
        let B = self.observer.x * self.pointing_vector.x.sin() * self.pointing_vector_rate.x;
        let C = r * r_rate - self.observer.x * (self.pointing_vector.x.sin().powi(2) * self.observer.velocity.x + self.pointing_vector.x.cos() * self.observer.x * self.pointing_vector.x.sin() * self.pointing_vector_rate.x);
        let D = (r.powi(2) - self.observer.x.powi(2) * self.pointing_vector.x.sin().powi(2)).sqrt();
        return A - B + (C / D);
    }

    fn solar_elongation(&self) -> f64 {
        let dot = self.observer.dot(&self.pointing_vector);
        let mag = self.observer.magnitude() * self.pointing_vector.magnitude();
        return dot / mag;
    }

    fn solar_elongation_rate_times_sin_solar_elongation(&self) -> f64 {
        let observer_r = self.observer.position.norm();
        let observer_rdot = self.observer.velocity.dot(&self.observer.position);
        let r_hat_dot = (observer_r * self.observer.velocity - self.observer.position * observer_rdot) / observer_r.powi(2);
        let elongation_rate = self.pointing_vector.dot(&r_hat_dot) + self.observer.position.unit().dot(&self.pointing_vector_rate);
        return elongation_rate;
    }

}


fn compute_pointing_vector(ra: f64, dec: f64) -> Vector3<f64> {
    let ra = ra.to_radians();
    let dec = dec.to_radians();
    let x = ra.cos() * dec.cos();
    let y = ra.sin() * dec.cos();
    let z = dec.sin();
    return Vector3::new(x, y, z);
}

