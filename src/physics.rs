use color_eyre::Result;

use self::euler1disot::Euler1DIsot;

mod euler1disot;

pub enum PhysicsMode {
    Euler1DIsot,
}

pub trait Physics {}

pub fn init_physics() -> Result<Box<dyn Physics>> {
    return Ok(Box::new(Euler1DIsot {}));
}
