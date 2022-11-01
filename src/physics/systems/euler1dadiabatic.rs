use super::super::Physics;

impl<const S: usize, const EQ: usize> Physics<S, EQ> {
    /// Updates conservative density, xi momentum, and energy according to the adiabatic Euler equations
    #[inline(always)]
    pub fn update_cons_euler1d_adiabatic(&mut self) {
        self.update_cons_euler1d_isot();
        self.calc_energy_euler1d();
    }

    /// Updates primitive density, xi velocity, and pressure according to the adiabatic Euler equations
    #[inline(always)]
    pub fn update_prim_euler1d_adiabatic(&mut self) {
        self.update_prim_euler1d_isot();
        self.calc_pressure_euler1d();
    }

    /// Calculates and updates the energy in [self.cons].
    #[inline(always)]
    pub fn calc_energy_euler1d(&mut self) {
        self.cons.row_mut(2).assign(
            &(&self.prim.row(2) * self.adiabatic_index.recip()
                + 0.5 * &self.prim.row(0) * &self.prim.row(1) * &self.prim.row(1)),
        );
    }

    /// Calculates and updates the pressure in [self.prim].
    #[inline(always)]
    pub fn calc_pressure_euler1d(&mut self) {
        self.prim.row_mut(2).assign(
            &(self.adiabatic_index
                * (&self.cons.row(2) - 0.5 / &self.cons.row(0) * &self.cons.row(1) * &self.cons.row(1))),
        );
    }
}
