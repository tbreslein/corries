use ndarray::ArrayD;

pub struct Variables {
    /// Primitive variables
    pub p: ArrayD<f64>,

    /// Conservative variables
    pub c: ArrayD<f64>,

    /// Speed of sound
    pub c_sound: ArrayD<f64>,

    /// Temperature
    pub temperature: ArrayD<f64>,

    /// Possible equation index for density
    pub jdensity: Option<usize>,

    /// Possible equation index for velocity in the xi direction
    pub jxivelocity: Option<usize>,

    /// Possible equation index for momentum in the xi direction
    pub jximomentum: Option<usize>,

    /// Number of equations
    pub n_equations: usize,
}
