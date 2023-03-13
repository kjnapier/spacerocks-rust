use spice;

use spacerocks::spacerock::SpaceRock;
use spacerocks::observatory::Observatory;
use spacerocks::constants::*;

use std::fs::File;
use std::io::Write;


fn main() {

    spice::furnsh("/home/kevin/Desktop/spice/latest_leapseconds.tls");
    spice::furnsh("/home/kevin/Desktop/spice/de440s.bsp");

    let mut epochs: Vec<f64> = Vec::new();
    let mut epoch = 2450000.38373;
    while epoch < 2450000.38373 + 28.0 {
        epochs.push(epoch);
        epoch += 0.01;
    }

    let mut file = File::create("/home/kevin/Desktop/positions.csv").unwrap();
    file.write_all(b"objid,x,y,z\n").unwrap();
    let objids = &vec!["Jupiter Barycenter", "Earth", "Sun", "Neptune Barycenter", "Mars Barycenter", 
                      "Venus Barycenter", "Mercury Barycenter", "Saturn Barycenter", "Uranus Barycenter", "Pluto Barycenter", "Moon"];
    for epoch in &epochs {
        for objid in objids {
            let mut body = SpaceRock::from_spice(objid, *epoch);
            body.change_frame("ECLIPJ2000");
            file.write_all(format!("{objid}, {x}, {y}, {z}\n", 
                                   objid=body.name, 
                                   x=body.position.x, 
                                   y=body.position.y, 
                                   z=body.position.z).as_bytes()).unwrap();
        }
    }

    let w84 = Observatory::from_coordinates(-30.00293494202556, -70.80642, 2207.0);
    let mut file = File::create("/home/kevin/Desktop/moon-sky.csv").unwrap();
    file.write_all(b"objid,epoch,ra,dec\n").unwrap();
    for epoch in &epochs {
        let mut body = SpaceRock::from_spice("moon", *epoch);
        println!("Epoch: {}", *epoch);
        let observer = w84.at(*epoch);
        let [ra, dec] = body.observe(&observer);
        file.write_all(format!("{objid}, {epoch}, {ra}, {dec}\n", 
                               objid=body.name, 
                               epoch=body.epoch, 
                               ra=ra.to_degrees(), 
                               dec=dec.to_degrees()).as_bytes()).unwrap();
        }

}