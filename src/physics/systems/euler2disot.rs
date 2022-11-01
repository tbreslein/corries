use super::super::Physics;

impl<const S: usize, const EQ: usize> Physics<S, EQ> {
    /// Updates conservative density, xi momentum, and eta momentum according to the adiabatic Euler equations
    #[inline(always)]
    pub fn update_cons_euler2d_isot(&mut self) {
        self.update_cons_euler1d_isot();
        self.calc_cons_linear_momentum_euler(2);
    }

    /// Updates primitive density, xi velocity, and eta velocity according to the adiabatic Euler equations
    #[inline(always)]
    pub fn update_prim_euler2d_isot(&mut self) {
        self.update_prim_euler1d_isot();
        self.calc_prim_linear_velocity_euler(2);
    }
}
