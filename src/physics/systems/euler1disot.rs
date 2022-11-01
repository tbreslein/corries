use super::super::Physics;

impl<const S: usize, const EQ: usize> Physics<S, EQ> {
    /// Updates conservative density and xi momentum according to the isothermal Euler equations
    #[inline(always)]
    pub fn update_cons_euler1d_isot(&mut self) {
        self.cons.row_mut(self.jdensity).assign(&self.prim.row(self.jdensity));
        self.calc_cons_linear_momentum_euler(self.jximomentum);
    }

    /// Updates primitive density and xi velocity according to the isothermal Euler equations
    #[inline(always)]
    pub fn update_prim_euler1d_isot(&mut self) {
        self.prim.row_mut(self.jdensity).assign(&self.cons.row(self.jdensity));
        self.calc_prim_linear_velocity_euler(self.jxivelocity);
    }

    /// Calculates and updates the linear momentum in [self.cons] at row `j`.
    #[inline(always)]
    pub fn calc_cons_linear_momentum_euler(&mut self, j: usize) {
        self.cons
            .row_mut(j)
            .assign(&(&self.prim.row(j) * &self.prim.row(self.jdensity)));
    }

    /// Calculates and updates the linear velocity in [self.prim] at row `j`.
    #[inline(always)]
    pub fn calc_prim_linear_velocity_euler(&mut self, j: usize) {
        self.prim
            .row_mut(j)
            .assign(&(&self.cons.row(j) / &self.cons.row(self.jdensity)));
    }

    /// Updates the speed of sound vector; no-op for isothermal physics
    #[inline(always)]
    pub fn update_c_sound_isot(&mut self) {
        self.c_sound_helper
            .assign(&(self.adiabatic_index * &self.prim.row(self.jpressure) / &self.prim.row(self.jdensity)));
        self.c_sound_helper.mapv_inplace(f64::sqrt);
        self.c_sound.assign(&self.c_sound_helper);
    }
}
