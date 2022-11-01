use super::super::Physics;

impl<const S: usize, const EQ: usize> Physics<S, EQ> {
    /// Updates conservative density and xi momentum according to the isothermal Euler equations
    #[inline(always)]
    pub fn update_cons_euler1d_isot(&mut self) {
        self.cons.row_mut(0).assign(&self.prim.row(0));
        self.calc_cons_linear_momentum_euler(1);
    }

    /// Updates primitive density and xi velocity according to the isothermal Euler equations
    #[inline(always)]
    pub fn update_prim_euler1d_isot(&mut self) {
        self.prim.row_mut(0).assign(&self.cons.row(0));
        self.calc_prim_linear_velocity_euler(1);
    }

    /// Calculates and updates the linear momentum in [self.cons] at row `j`.
    #[inline(always)]
    pub fn calc_cons_linear_momentum_euler(&mut self, j: usize) {
        self.cons.row_mut(j).assign(&(&self.prim.row(j) * &self.prim.row(0)));
    }

    /// Calculates and updates the linear velocity in [self.prim] at row `j`.
    #[inline(always)]
    pub fn calc_prim_linear_velocity_euler(&mut self, j: usize) {
        self.prim.row_mut(j).assign(&(&self.cons.row(j) / &self.cons.row(0)));
    }
}
