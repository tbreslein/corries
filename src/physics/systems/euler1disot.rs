use ndarray::{azip, Array2};

use super::super::Physics;

impl<const S: usize, const EQ: usize> Physics<S, EQ> {
    /// Updates conservative density and xi momentum according to the isothermal Euler equations
    #[inline(always)]
    pub fn update_cons_euler1d_isot(&mut self) {
        self.cons.row_mut(self.jdensity).assign(&self.prim.row(self.jdensity));
        self.update_cons_linear_momentum_euler(self.jximomentum);
    }

    /// Updates primitive density and xi velocity according to the isothermal Euler equations
    #[inline(always)]
    pub fn update_prim_euler1d_isot(&mut self) {
        self.prim.row_mut(self.jdensity).assign(&self.cons.row(self.jdensity));
        self.update_prim_linear_velocity_euler(self.jxivelocity);
    }

    /// Calculates and updates the linear momentum in [self.cons] at row `j`.
    #[inline(always)]
    pub fn update_cons_linear_momentum_euler(&mut self, j: usize) {
        self.cons
            .row_mut(j)
            .assign(&(&self.prim.row(j) * &self.prim.row(self.jdensity)));
    }

    /// Calculates and updates the linear velocity in [self.prim] at row `j`.
    #[inline(always)]
    pub fn update_prim_linear_velocity_euler(&mut self, j: usize) {
        self.prim
            .row_mut(j)
            .assign(&(&self.cons.row(j) / &self.cons.row(self.jdensity)));
    }

    /// Updates the speed of sound vector; no-op for isothermal physics
    #[inline(always)]
    pub fn update_c_sound_isot(&mut self) {
        return;
    }

    /// Calculates the physical flux for isothermal 1D Euler
    #[inline(always)]
    pub fn calc_physical_flux_euler1d_isot(&self, flux: &mut Array2<f64>) {
        self.calc_density_flux_euler1d_isot(flux);
        self.calc_xi_momentum_flux_euler1d_isot(flux);
    }

    /// Calculates the density flux for isothermal 1D Euler
    #[inline(always)]
    pub fn calc_density_flux_euler1d_isot(&self, flux: &mut Array2<f64>) {
        azip!(
            (
                dens_flux in flux.row_mut(self.jdensity),
                &ximom in self.cons.row(self.jximomentum)
            )
            *dens_flux = ximom);
    }

    /// Calculates the xi momentum flux for isothermal 1D Euler
    #[inline(always)]
    pub fn calc_xi_momentum_flux_euler1d_isot(&self, flux: &mut Array2<f64>) {
        azip!(
            (
                ximom_flux in flux.row_mut(self.jximomentum),
                &dens in self.cons.row(self.jdensity),
                &xivel in self.prim.row(self.jxivelocity),
                &cs in &self.c_sound
            )
            *ximom_flux = dens * (xivel * xivel + cs * cs));
    }
}
