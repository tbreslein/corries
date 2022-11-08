use ndarray::{Array2, Zip};

use super::super::Physics;

impl<const S: usize, const EQ: usize> Physics<S, EQ> {
    /// Updates conservative density, xi momentum, and energy according to the adiabatic Euler equations
    #[inline(always)]
    pub fn update_cons_euler1d_adiabatic(&mut self) {
        self.update_cons_euler1d_isot();
        self.update_energy_euler1d();
    }

    /// Updates primitive density, xi velocity, and pressure according to the adiabatic Euler equations
    #[inline(always)]
    pub fn update_prim_euler1d_adiabatic(&mut self) {
        self.update_prim_euler1d_isot();
        self.update_pressure_euler1d();
    }

    /// Calculates and updates the energy in [self.cons].
    #[inline(always)]
    pub fn update_energy_euler1d(&mut self) {
        Zip::from(self.cons.row_mut(self.jenergy))
            .and(self.prim.row(self.jpressure))
            .and(self.prim.row(self.jdensity))
            .and(self.prim.row(self.jxivelocity))
            .for_each(|energy, &pressure, &density, &velocity| {
                *energy = pressure * (self.adiabatic_index - 1.0).recip() + 0.5 * density * velocity * velocity
            });
    }

    /// Calculates and updates the pressure in [self.prim].
    #[inline(always)]
    pub fn update_pressure_euler1d(&mut self) {
        Zip::from(self.prim.row_mut(self.jpressure))
            .and(self.cons.row(self.jenergy))
            .and(self.cons.row(self.jdensity))
            .and(self.cons.row(self.jximomentum))
            .for_each(|pressure, &energy, &density, &momentum| {
                *pressure = (self.adiabatic_index - 1.0) * (energy - 0.5 / density * momentum * momentum)
            });
    }

    /// Updates the speed of sound vector; no-op for isothermal physics
    #[inline(always)]
    pub fn update_c_sound_adiabatic(&mut self) {
        self.c_sound_helper
            .assign(&(self.adiabatic_index * &self.prim.row(self.jpressure) / &self.prim.row(self.jdensity)));
        self.c_sound_helper.mapv_inplace(f64::sqrt);
        self.c_sound.assign(&self.c_sound_helper);
    }

    /// Calculates the physical flux for adiabatic 1D Euler
    #[inline(always)]
    pub fn calc_physical_flux_euler1d_adiabatic(&self, flux: &mut Array2<f64>) {
        self.calc_density_flux_euler1d_isot(flux);
        self.calc_xi_momentum_flux_euler1d_isot(flux);
        self.calc_energy_flux_euler1d_isot(flux);
    }

    /// Calculates the xi momentum flux for adiabatic 1D Euler
    #[inline(always)]
    pub fn calc_xi_momentum_flux_euler1d_adiabatic(&self, flux: &mut Array2<f64>) {
        Zip::from(flux.row_mut(self.jximomentum))
            .and(self.prim.row(self.jxivelocity))
            .and(self.cons.row(self.jximomentum))
            .and(self.prim.row(self.jpressure))
            .for_each(|xi_momentum_flux, &xi_velocity, &xi_momentum, &pressure| {
                *xi_momentum_flux = xi_velocity * xi_momentum + pressure
            });
    }

    /// Calculates the energy flux for adiabatic 1D Euler
    #[inline(always)]
    pub fn calc_energy_flux_euler1d_isot(&self, flux: &mut Array2<f64>) {
        Zip::from(flux.row_mut(self.jenergy))
            .and(self.cons.row(self.jenergy))
            .and(self.prim.row(self.jpressure))
            .and(self.prim.row(self.jxivelocity))
            .for_each(|energy_flux, &energy, &pressure, &xi_velocity| *energy_flux = (energy + pressure) * xi_velocity);
    }
}
