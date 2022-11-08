use ndarray::{Array2, Zip};

use super::super::Physics;

impl<const S: usize, const EQ: usize> Physics<S, EQ> {
    /// Updates conservative density, xi momentum, and eta momentum according to the adiabatic Euler equations
    #[inline(always)]
    pub fn update_cons_euler2d_isot(&mut self) {
        self.update_cons_euler1d_isot();
        self.update_cons_linear_momentum_euler(self.jetamomentum);
    }

    /// Updates primitive density, xi velocity, and eta velocity according to the adiabatic Euler equations
    #[inline(always)]
    pub fn update_prim_euler2d_isot(&mut self) {
        self.update_prim_euler1d_isot();
        self.update_prim_linear_velocity_euler(self.jetavelocity);
    }

    /// Calculates the physical flux for isothermal 2D Euler
    #[inline(always)]
    pub fn calc_physical_flux_euler2d_isot(&self, flux: &mut Array2<f64>) {
        self.calc_physical_flux_euler1d_isot(flux);
        self.calc_eta_momentum_flux_euler2d_isot(flux);
    }

    /// Calculates the eta momentum flux for isothermal 2D Euler
    #[inline(always)]
    pub fn calc_eta_momentum_flux_euler2d_isot(&self, flux: &mut Array2<f64>) {
        Zip::from(flux.row_mut(self.jetamomentum))
            .and(self.prim.row(self.jxivelocity))
            .and(self.cons.row(self.jetamomentum))
            .for_each(|eta_momentum_flux, &xi_velocity, &eta_momentum| *eta_momentum_flux = xi_velocity * eta_momentum);
    }
}
