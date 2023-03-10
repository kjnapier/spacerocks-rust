pub const KM_TO_AU: f64 = 1.0 / 149_597_870.700;
pub const SECONDS_PER_DAY: f64 = 86_400.0;

pub const EQUAT_RAD: f64 = 6378137.0;
pub const FLATTEN: f64 = 1.0 / 298.257223563;
pub const O_M_FLATTEN: f64 = 1.0 - FLATTEN;
pub const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;

pub const MU_BARY: f64 = 0.00029630927493457475;
pub const SPEED_OF_LIGHT: f64 = 173.14463268466926; // speed of light in au/day