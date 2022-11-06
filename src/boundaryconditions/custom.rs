// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [CustomBoundaryConditions] which applies a set of custom boundary conditions to a
//! `Physics` object.

use crate::{
    config::{CustomBoundaryMode, PhysicsVariable},
    mesh::Mesh,
    physics::Physics,
};

use super::{BoundaryCondition, Direction};

pub struct CustomBoundaryConditions {
    direction: Direction,
    modes: Vec<(usize, CustomBoundaryMode)>,
}

impl<const S: usize, const EQ: usize> BoundaryCondition<S, EQ> for CustomBoundaryConditions {
    fn apply(&mut self, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
        match self.direction {
            Direction::West => self.modes.iter().for_each(|(j, mode)| match mode {
                CustomBoundaryMode::Extrapolate => extrapolate_west(*j, u, mesh),
                CustomBoundaryMode::ExtrapolateDensityKepler => extrapolate_density_kepler_west(*j, u, mesh),
                CustomBoundaryMode::ExtrapolateEtaVelocityKepler => extrapolate_etavel_kepler_west(*j, u, mesh),
                CustomBoundaryMode::NearZero => near_zero_west(*j, u),
                CustomBoundaryMode::NoGradients => no_gradients_west(*j, u),
                CustomBoundaryMode::OutFlowExtrapolate => outflow_extrapolate_west(*j, u, mesh),
                CustomBoundaryMode::OutflowNoGradients => outflow_no_gradients_west(*j, u),
                CustomBoundaryMode::Reflecting => reflecting_west(*j, u),
            }),
            Direction::East => self.modes.iter().for_each(|(j, mode)| match mode {
                CustomBoundaryMode::Extrapolate => extrapolate_east(*j, u, mesh),
                CustomBoundaryMode::ExtrapolateDensityKepler => extrapolate_density_kepler_east(*j, u, mesh),
                CustomBoundaryMode::ExtrapolateEtaVelocityKepler => extrapolate_etavel_kepler_east(*j, u, mesh),
                CustomBoundaryMode::NearZero => near_zero_east(*j, u),
                CustomBoundaryMode::NoGradients => no_gradients_east(*j, u),
                CustomBoundaryMode::OutFlowExtrapolate => outflow_extrapolate_east(*j, u, mesh),
                CustomBoundaryMode::OutflowNoGradients => outflow_no_gradients_east(*j, u),
                CustomBoundaryMode::Reflecting => reflecting_east(*j, u),
            }),
        }
    }
}

impl CustomBoundaryConditions {
    pub fn new<const S: usize, const EQ: usize>(
        direction: Direction,
        modes: &Vec<(PhysicsVariable, CustomBoundaryMode)>,
        u: &Physics<S, EQ>,
    ) -> Self {
        return Self {
            direction,
            modes: modes
                .iter()
                .map(|(var, mode)| (u.convert_physics_variable_to_index(*var).unwrap(), *mode))
                .collect(),
        };
    }
}

fn extrapolate_west<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
    if mesh.is_logarithmic {
        if !(1.0 / u.prim[[j, 2]]).is_finite() || !(1.0 / u.prim[[j, 3]]).is_finite() {
            for i in 0..=1 {
                u.prim[[j, i]] = 0.0;
            }
        } else {
            for i in 1..=2 {
                u.prim[[j, 2 - i]] = u.prim[[j, 2]] * (u.prim[[j, 2]] / u.prim[[j, 3]]).abs().powi(i as i32);
            }
        }
    } else {
        for i in 1..=2 {
            u.prim[[j, 2 - i]] = (i + 1) as f64 * u.prim[[j, 2]] - i as f64 * u.prim[[j, 3]];
        }
    }
}

fn extrapolate_east<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
    if mesh.is_logarithmic {
        if !(1.0 / u.prim[[j, S - 3]]).is_finite() || !(1.0 / u.prim[[j, S - 4]]).is_finite() {
            for i in S - 2..=S - 1 {
                u.prim[[j, i]] = 0.0;
            }
        } else {
            for i in 1..=2 {
                u.prim[[j, S - 3 + i]] =
                    u.prim[[j, S - 3]] * (u.prim[[j, S - 3]] / u.prim[[j, S - 4]]).abs().powi(i as i32);
            }
        }
    } else {
        for i in 1..=2 {
            u.prim[[j, S - 3 + i]] = (i + 1) as f64 * u.prim[[j, S - 3]] - i as f64 * u.prim[[j, S - 4]];
        }
    }
}

fn extrapolate_density_kepler_west<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
    for i in (0..=1).rev() {
        u.prim[[j, i]] = u.prim[[j, i + 1]] / (mesh.xi_cent[i] * mesh.xi_cent_inv[i + 1]).sqrt();
    }
}

fn extrapolate_density_kepler_east<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
    for i in S - 2..=S - 1 {
        u.prim[[j, i]] = u.prim[[j, i - 1]] / (mesh.xi_cent[i] * mesh.xi_cent_inv[i - 1]).sqrt();
    }
}

fn extrapolate_etavel_kepler_west<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
    for i in (0..=1).rev() {
        u.prim[[j, i]] = u.prim[[j, i + 1]] / (mesh.xi_cent[i] * mesh.xi_cent_inv[i + 1]).sqrt();
    }
}

fn extrapolate_etavel_kepler_east<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
    for i in S - 2..=S - 1 {
        u.prim[[j, i]] = u.prim[[j, i - 1]] / (mesh.xi_cent[i] * mesh.xi_cent_inv[i - 1]).sqrt();
    }
}

fn near_zero_west<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>) {
    for i in 0..=1 {
        u.prim.row_mut(j)[i] = 1.0e-10 * u.prim.row(j)[2];
    }
}

fn near_zero_east<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>) {
    for i in S - 2..=S - 1 {
        u.prim.row_mut(j)[i] = 1.0e-10 * u.prim.row(j)[S - 3];
    }
}

fn no_gradients_west<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>) {
    for i in 0..=1 {
        u.prim.row_mut(j)[i] = u.prim.row(j)[2];
    }
}

fn no_gradients_east<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>) {
    for i in S - 2..=S - 1 {
        u.prim.row_mut(j)[i] = u.prim.row(j)[S - 3];
    }
}

fn outflow_no_gradients_west<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>) {
    if u.prim[[j, 2]].is_sign_positive() {
        reflecting_west(j, u);
    } else {
        no_gradients_west(j, u);
    }
}

fn outflow_no_gradients_east<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>) {
    if u.prim[[j, S - 3]].is_sign_negative() {
        reflecting_east(j, u);
    } else {
        no_gradients_east(j, u);
    }
}

fn outflow_extrapolate_west<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
    if u.prim[[j, 2]].is_sign_positive() {
        reflecting_west(j, u);
    } else {
        extrapolate_west(j, u, mesh);
    }
}

fn outflow_extrapolate_east<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
    if u.prim[[j, S - 3]].is_sign_negative() {
        reflecting_east(j, u);
    } else {
        extrapolate_east(j, u, mesh);
    }
}

fn reflecting_west<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>) {
    for i in 1..=2 {
        u.prim[[j, 2 - i]] = -1.0 * u.prim[[j, i + 1]]
    }
}

fn reflecting_east<const S: usize, const EQ: usize>(j: usize, u: &mut Physics<S, EQ>) {
    for i in 1..=2 {
        u.prim[[j, S - 3 + i]] = -1.0 * u.prim[[j, S - 2 - i]]
    }
}


#[cfg(test)]
mod tests {
    use crate::{config::{meshconfig::{MeshConfig, MeshMode}, physicsconfig::{PhysicsConfig, PhysicsMode}}, units::UnitsMode};

    use super::*;
    use approx::assert_relative_eq;
    use ndarray::Array2;
    const S: usize = 6;
    const EQ: usize = 3;
    const MESHCONFIG: MeshConfig = MeshConfig {
        mode: MeshMode::Cartesian,
        xi_in: 2.0,
        xi_out: 3.0,
        ratio_disk: 1.0,
    };
    const PHYSICSCONFIG: PhysicsConfig = PhysicsConfig {mode: PhysicsMode::Euler1DAdiabatic, units_mode: UnitsMode::SI, adiabatic_index: 1.4, c_sound_0: 1.0};

    #[test]
    fn extrapolation_test() {
        // Also test this on log meshes
        let mesh: Mesh<S> = Mesh::new(&MESHCONFIG).unwrap();
        let mut u: Physics<S, EQ> = Physics::new(&PHYSICSCONFIG);
        u.prim.assign(&Array2::from_shape_vec((EQ, S), vec![
            0.0, 0.0, 0.25, 0.25, 0.0, 0.0,
            0.0, 0.0, -2.0, -1.5, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]).unwrap());
        let u_prim_expect = Array2::from_shape_vec((EQ, S), vec![
            0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
            -3.0, -2.5, -2.0, -1.5, -1.0, -0.5,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]).unwrap();
        for j in 0..EQ {
            extrapolate_west(j, &mut u, &mesh);
            extrapolate_east(j, &mut u, &mesh);
        }
        assert_relative_eq!(u.prim, u_prim_expect, max_relative = 1.0e-12);
    }

    #[test]
    fn near_zero_test() {
        let mut u: Physics<S, EQ> = Physics::new(&PHYSICSCONFIG);
        u.prim.assign(&Array2::from_shape_vec((EQ, S), vec![
            0.0, 0.0, 0.25, 0.25, 0.0, 0.0,
            0.0, 0.0, -2.0, -1.5, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]).unwrap());
        let u_prim_expect = Array2::from_shape_vec((EQ, S), vec![
            0.25e-10, 0.25e-10, 0.25, 0.25, 0.25e-10, 0.25e-10,
            -2.0e-10, -2.0e-10, -2.0, -1.5, -1.5e-10, -1.5e-10,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]).unwrap();
        for j in 0..EQ {
            near_zero_west(j, &mut u);
            near_zero_east(j, &mut u);
        }
        assert_relative_eq!(u.prim, u_prim_expect, max_relative = 1.0e-12);
    }

    #[test]
    fn no_gradients_test() {
        let mut u: Physics<S, EQ> = Physics::new(&PHYSICSCONFIG);
        u.prim.assign(&Array2::from_shape_vec((EQ, S), vec![
            0.0, 0.0, 0.25, 0.25, 0.0, 0.0,
            0.0, 0.0, -2.0, -1.5, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]).unwrap());

        let u_prim_expect = Array2::from_shape_vec((EQ, S), vec![
            0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
            -2.0, -2.0, -2.0, -1.5, -1.5, -1.5,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]).unwrap();

        for j in 0..EQ {
            no_gradients_west(j, &mut u);
            no_gradients_east(j, &mut u);
        }

        assert_relative_eq!(u.prim, u_prim_expect, max_relative = 1.0e-12);
    }
    
    #[test]
    fn outflow_no_gradients_test() {
        let mut u: Physics<S, EQ> = Physics::new(&PHYSICSCONFIG);
        u.prim.assign(&Array2::from_shape_vec((EQ, S), vec![
            0.0, 0.0, 0.25, 0.25, 0.0, 0.0,
            0.0, 0.0, -2.0, -1.5, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]).unwrap());
        let u_prim_expect = Array2::from_shape_vec((EQ, S), vec![
            -0.25, -0.25, 0.25, 0.25, 0.25, 0.25,
            -2.0, -2.0, -2.0, -1.5, 1.5, 2.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]).unwrap();
        for j in 0..EQ {
            outflow_no_gradients_west(j, &mut u);
            outflow_no_gradients_east(j, &mut u);
        }
        assert_relative_eq!(u.prim, u_prim_expect, max_relative = 1.0e-12);
    }

    #[test]
    fn outflow_extrapolate_test() {
        // Also test this on log meshes
        let mesh: Mesh<S> = Mesh::new(&MESHCONFIG).unwrap();
        let mut u: Physics<S, EQ> = Physics::new(&PHYSICSCONFIG);
        u.prim.assign(&Array2::from_shape_vec((EQ, S), vec![
            0.0, 0.0, 0.25, 0.25, 0.0, 0.0,
            0.0, 0.0, -2.0, -1.5, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]).unwrap());
        let u_prim_expect = Array2::from_shape_vec((EQ, S), vec![
            -0.25, -0.25, 0.25, 0.25, 0.25, 0.25,
            -3.0, -2.5, -2.0, -1.5, 1.5, 2.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]).unwrap();
        for j in 0..EQ {
            outflow_extrapolate_west(j, &mut u, &mesh);
            outflow_extrapolate_east(j, &mut u, &mesh);
        }
        assert_relative_eq!(u.prim, u_prim_expect, max_relative = 1.0e-12);
    }

    #[test]
    fn reflecting_test() {
        // Also test this on log meshes
        let mut u: Physics<S, EQ> = Physics::new(&PHYSICSCONFIG);
        u.prim.assign(&Array2::from_shape_vec((EQ, S), vec![
            0.0, 0.0, 0.25, 0.25, 0.0, 0.0,
            0.0, 0.0, -2.0, -1.5, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]).unwrap());
        let u_prim_expect = Array2::from_shape_vec((EQ, S), vec![
            -0.25, -0.25, 0.25, 0.25, -0.25, -0.25,
            1.5, 2.0, -2.0, -1.5, 1.5, 2.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ]).unwrap();
        for j in 0..EQ {
            reflecting_west(j, &mut u);
            reflecting_east(j, &mut u);
        }
        assert_relative_eq!(u.prim, u_prim_expect, max_relative = 1.0e-12);
    }
}
