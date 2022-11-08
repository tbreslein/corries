use ndarray::{Array2, Zip};

use super::super::Physics;

impl<const S: usize, const EQ: usize> Physics<S, EQ> {
    /// Updates conservative density and xi momentum according to the isothermal Euler equations
    #[inline(always)]
    pub fn update_cons_euler1d_isot(&mut self) {
        Zip::from(self.cons.row_mut(self.jdensity))
            .and(self.prim.row(self.jdensity))
            .for_each(|cons_dens, &prim_dens| *cons_dens = prim_dens);

        self.update_cons_linear_momentum_euler(self.jximomentum);
    }

    /// Updates primitive density and xi velocity according to the isothermal Euler equations
    #[inline(always)]
    pub fn update_prim_euler1d_isot(&mut self) {
        Zip::from(self.prim.row_mut(self.jdensity))
            .and(self.cons.row(self.jdensity))
            .for_each(|prim_dens, &cons_dens| *prim_dens = cons_dens);

        self.update_prim_linear_velocity_euler(self.jxivelocity);
    }

    /// Calculates and updates the linear momentum in [self.cons] at row `j`.
    #[inline(always)]
    pub fn update_cons_linear_momentum_euler(&mut self, j: usize) {
        Zip::from(self.cons.row_mut(j))
            .and(self.prim.row(j))
            .and(self.prim.row(self.jdensity))
            .for_each(|momentum, &velocity, &density| *momentum = velocity * density);
    }

    /// Calculates and updates the linear velocity in [self.prim] at row `j`.
    #[inline(always)]
    pub fn update_prim_linear_velocity_euler(&mut self, j: usize) {
        Zip::from(self.prim.row_mut(j))
            .and(self.cons.row(j))
            .and(self.cons.row(self.jdensity))
            .for_each(|velocity, &momentum, &density| *velocity = momentum / density);
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
        Zip::from(flux.row_mut(self.jdensity))
            .and(self.cons.row(self.jximomentum))
            .for_each(|dens_flux, &xi_momentum| *dens_flux = xi_momentum);
    }

    /// Calculates the xi momentum flux for isothermal 1D Euler
    #[inline(always)]
    pub fn calc_xi_momentum_flux_euler1d_isot(&self, flux: &mut Array2<f64>) {
        Zip::from(flux.row_mut(self.jximomentum))
            .and(self.cons.row(self.jdensity))
            .and(self.prim.row(self.jxivelocity))
            .and(&self.c_sound)
            .for_each(|xi_momentum_flux, &density, &xi_velocity, &cs| {
                *xi_momentum_flux = density * (xi_velocity * xi_velocity + cs * cs)
            });
    }
}
