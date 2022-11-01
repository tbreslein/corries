use ndarray::{azip, Array2};

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
        self.cons.row_mut(self.jenergy).assign(
            &(&self.prim.row(self.jpressure) * self.adiabatic_index.recip()
                + 0.5
                    * &self.prim.row(self.jdensity)
                    * &self.prim.row(self.jxivelocity)
                    * &self.prim.row(self.jxivelocity)),
        );
    }

    /// Calculates and updates the pressure in [self.prim].
    #[inline(always)]
    pub fn update_pressure_euler1d(&mut self) {
        self.prim.row_mut(self.jpressure).assign(
            &(self.adiabatic_index
                * (&self.cons.row(self.jenergy)
                    - 0.5 / &self.cons.row(self.jdensity)
                        * &self.cons.row(self.jximomentum)
                        * &self.cons.row(self.jximomentum))),
        );
    }

    /// Updates the speed of sound vector; no-op for isothermal physics
    #[inline(always)]
    pub fn update_c_sound_adiabatic(&mut self) {
        return;
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
        azip!(
            (
                ximom_flux in flux.row_mut(self.jximomentum),
                &xivel in self.prim.row(self.jxivelocity),
                &ximom in self.cons.row(self.jximomentum),
                &p in self.prim.row(self.jpressure)
            )
            *ximom_flux = xivel * ximom + p);
    }

    /// Calculates the energy flux for adiabatic 1D Euler
    #[inline(always)]
    pub fn calc_energy_flux_euler1d_isot(&self, flux: &mut Array2<f64>) {
        azip!(
            (
                energy_flux in flux.row_mut(self.jenergy),
                &e in self.cons.row(self.jenergy),
                &p in self.prim.row(self.jpressure),
                &xivel in self.prim.row(self.jxivelocity)
            )
            *energy_flux = (e + p) * xivel);
    }
}
