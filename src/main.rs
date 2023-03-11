use spice;
use spacerocks::spacerock::SpaceRock;
use spacerocks::calc_E_from_M::calc_E_from_M;
use spacerocks::constants::*;
use std::fs::File;
use std::io::Write;


fn main() {

    spice::furnsh("/Users/kjnapier/Desktop/research/spacerocks/spacerocks/data/spice/latest_leapseconds.tls");
    spice::furnsh("/Users/kjnapier/Desktop/research/spacerocks/spacerocks/data/spice/de423.bsp");
    spice::furnsh("/Users/kjnapier/Desktop/research/spacerocks/spacerocks/data/spice/gm_Horizons.pck");

    let mut epochs: Vec<f64> = Vec::new();
    let mut epoch = 2450000.38373;
    while epoch < 2450000.38373 + 200.0 * 365.25 {
        epochs.push(epoch);
        epoch += 10.0;
    }

    let mut file = File::create("/Users/kjnapier/Desktop/positions.csv").unwrap();
    file.write_all(b"objid,x,y,z\n").unwrap();
    let objids = &vec!["Jupiter Barycenter", "Earth", "Sun", "Neptune Barycenter", "Mars Barycenter", 
                      "Venus Barycenter", "Mercury Barycenter", "Saturn Barycenter", "Uranus Barycenter", "Pluto Barycenter", "Moon"];
    for epoch in epochs {
        for objid in objids {
            let mut body = SpaceRock::from_spice(objid, epoch);
            body.change_frame("ECLIPJ2000");
            file.write_all(format!("{objid}, {x}, {y}, {z}\n", 
                                   objid=body.name, 
                                   x=body.position.x, 
                                   y=body.position.y, 
                                   z=body.position.z).as_bytes()).unwrap();
        }
    }

}