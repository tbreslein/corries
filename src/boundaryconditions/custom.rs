// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [CustomBoundaryConditions] which applies a set of custom boundary conditions to a
//! [Variables] object.

use crate::{config::CustomBoundaryMode, mesh::Mesh, variables::Variables};

use super::{BoundaryCondition, Direction};

/// This is a set of boundary conditions, where the exact condition is set for each equation
/// separately. This is useful, where you use generic boundary conditions and want to swap them out
/// quickly, whithout defining a new boundary condition set each time.
///
/// Check [CustomBoundaryMode] for the types of boundary conditions available.
pub struct CustomBoundaryConditions {
    direction: Direction,
    modes: Vec<(usize, CustomBoundaryMode)>,
}

impl<const E: usize, const S: usize> BoundaryCondition<E, S> for CustomBoundaryConditions {
    fn apply(&mut self, vars: &mut Variables<E, S>, mesh: &Mesh<S>) {
        match self.direction {
            Direction::West => self.modes.iter().for_each(|(j, mode)| match mode {
                CustomBoundaryMode::Extrapolate => extrapolate_west(*j, vars, mesh),
                CustomBoundaryMode::ExtrapolateDensityKepler => extrapolate_density_kepler_west(*j, vars, mesh),
                CustomBoundaryMode::ExtrapolateEtaVelocityKepler => extrapolate_etavel_kepler_west(*j, vars, mesh),
                CustomBoundaryMode::NearZero => near_zero_west(*j, vars),
                CustomBoundaryMode::NoGradients => no_gradients_west(*j, vars),
                CustomBoundaryMode::OutFlowExtrapolate => outflow_extrapolate_west(*j, vars, mesh),
                CustomBoundaryMode::OutflowNoGradients => outflow_no_gradients_west(*j, vars),
                CustomBoundaryMode::Reflecting => reflecting_west(*j, vars),
            }),
            Direction::East => self.modes.iter().for_each(|(j, mode)| match mode {
                CustomBoundaryMode::Extrapolate => extrapolate_east(*j, vars, mesh),
                CustomBoundaryMode::ExtrapolateDensityKepler => extrapolate_density_kepler_east(*j, vars, mesh),
                CustomBoundaryMode::ExtrapolateEtaVelocityKepler => extrapolate_etavel_kepler_east(*j, vars, mesh),
                CustomBoundaryMode::NearZero => near_zero_east(*j, vars),
                CustomBoundaryMode::NoGradients => no_gradients_east(*j, vars),
                CustomBoundaryMode::OutFlowExtrapolate => outflow_extrapolate_east(*j, vars, mesh),
                CustomBoundaryMode::OutflowNoGradients => outflow_no_gradients_east(*j, vars),
                CustomBoundaryMode::Reflecting => reflecting_east(*j, vars),
            }),
        }
    }
}

impl CustomBoundaryConditions {
    /// Constructs a new [CustomBoundaryConditions] object
    pub fn new(direction: Direction, modes: &[(usize, CustomBoundaryMode)]) -> Self {
        return Self {
            direction,
            modes: modes.to_vec(),
        };
    }
}

fn extrapolate_west<const E: usize, const S: usize>(j: usize, vars: &mut Variables<E, S>, mesh: &Mesh<S>) {
    if mesh.is_logarithmic {
        if !(1.0 / vars.prim[[j, 2]]).is_finite() || !(1.0 / vars.prim[[j, 3]]).is_finite() {
            for i in 0..=1 {
                vars.prim[[j, i]] = 0.0;
            }
        } else {
            for i in 1..=2 {
                vars.prim[[j, 2 - i]] =
                    vars.prim[[j, 2]] * (vars.prim[[j, 2]] / vars.prim[[j, 3]]).abs().powi(i as i32);
            }
        }
    } else {
        for i in 1..=2 {
            vars.prim[[j, 2 - i]] = (i + 1) as f64 * vars.prim[[j, 2]] - i as f64 * vars.prim[[j, 3]];
        }
    }
}

fn extrapolate_east<const E: usize, const S: usize>(j: usize, vars: &mut Variables<E, S>, mesh: &Mesh<S>) {
    if mesh.is_logarithmic {
        if !(1.0 / vars.prim[[j, S - 3]]).is_finite() || !(1.0 / vars.prim[[j, S - 4]]).is_finite() {
            for i in S - 2..=S - 1 {
                vars.prim[[j, i]] = 0.0;
            }
        } else {
            for i in 1..=2 {
                vars.prim[[j, S - 3 + i]] =
                    vars.prim[[j, S - 3]] * (vars.prim[[j, S - 3]] / vars.prim[[j, S - 4]]).abs().powi(i as i32);
            }
        }
    } else {
        for i in 1..=2 {
            vars.prim[[j, S - 3 + i]] = (i + 1) as f64 * vars.prim[[j, S - 3]] - i as f64 * vars.prim[[j, S - 4]];
        }
    }
}

fn extrapolate_density_kepler_west<const E: usize, const S: usize>(
    j: usize,
    vars: &mut Variables<E, S>,
    mesh: &Mesh<S>,
) {
    for i in (0..=1).rev() {
        vars.prim[[j, i]] = vars.prim[[j, i + 1]] / (mesh.xi_cent[i] * mesh.xi_cent_inv[i + 1]).sqrt();
    }
}

fn extrapolate_density_kepler_east<const E: usize, const S: usize>(
    j: usize,
    vars: &mut Variables<E, S>,
    mesh: &Mesh<S>,
) {
    for i in S - 2..=S - 1 {
        vars.prim[[j, i]] = vars.prim[[j, i - 1]] / (mesh.xi_cent[i] * mesh.xi_cent_inv[i - 1]).sqrt();
    }
}

fn extrapolate_etavel_kepler_west<const E: usize, const S: usize>(
    j: usize,
    vars: &mut Variables<E, S>,
    mesh: &Mesh<S>,
) {
    for i in (0..=1).rev() {
        vars.prim[[j, i]] = vars.prim[[j, i + 1]] / (mesh.xi_cent[i] * mesh.xi_cent_inv[i + 1]).sqrt();
    }
}

fn extrapolate_etavel_kepler_east<const E: usize, const S: usize>(
    j: usize,
    vars: &mut Variables<E, S>,
    mesh: &Mesh<S>,
) {
    for i in S - 2..=S - 1 {
        vars.prim[[j, i]] = vars.prim[[j, i - 1]] / (mesh.xi_cent[i] * mesh.xi_cent_inv[i - 1]).sqrt();
    }
}

fn near_zero_west<const E: usize, const S: usize>(j: usize, vars: &mut Variables<E, S>) {
    for i in 0..=1 {
        vars.prim[[j, i]] = 1.0e-10 * vars.prim[[j, 2]];
    }
}

fn near_zero_east<const E: usize, const S: usize>(j: usize, vars: &mut Variables<E, S>) {
    for i in S - 2..=S - 1 {
        vars.prim[[j, i]] = 1.0e-10 * vars.prim[[j, S - 3]];
    }
}

fn no_gradients_west<const E: usize, const S: usize>(j: usize, vars: &mut Variables<E, S>) {
    for i in 0..=1 {
        vars.prim[[j, i]] = vars.prim[[j, 2]];
    }
}

fn no_gradients_east<const E: usize, const S: usize>(j: usize, vars: &mut Variables<E, S>) {
    for i in S - 2..=S - 1 {
        vars.prim[[j, i]] = vars.prim[[j, S - 3]];
    }
}

fn outflow_no_gradients_west<const E: usize, const S: usize>(j: usize, vars: &mut Variables<E, S>) {
    if vars.prim[[j, 2]].is_sign_positive() {
        reflecting_west(j, vars);
    } else {
        no_gradients_west(j, vars);
    }
}

fn outflow_no_gradients_east<const E: usize, const S: usize>(j: usize, vars: &mut Variables<E, S>) {
    if vars.prim[[j, S - 3]].is_sign_negative() {
        reflecting_east(j, vars);
    } else {
        no_gradients_east(j, vars);
    }
}

fn outflow_extrapolate_west<const E: usize, const S: usize>(j: usize, vars: &mut Variables<E, S>, mesh: &Mesh<S>) {
    if vars.prim[[j, 2]].is_sign_positive() {
        reflecting_west(j, vars);
    } else {
        extrapolate_west(j, vars, mesh);
    }
}

fn outflow_extrapolate_east<const E: usize, const S: usize>(j: usize, vars: &mut Variables<E, S>, mesh: &Mesh<S>) {
    if vars.prim[[j, S - 3]].is_sign_negative() {
        reflecting_east(j, vars);
    } else {
        extrapolate_east(j, vars, mesh);
    }
}

fn reflecting_west<const E: usize, const S: usize>(j: usize, vars: &mut Variables<E, S>) {
    for i in 1..=2 {
        vars.prim[[j, 2 - i]] = -1.0 * vars.prim[[j, i + 1]];
    }
}

fn reflecting_east<const E: usize, const S: usize>(j: usize, vars: &mut Variables<E, S>) {
    for i in 1..=2 {
        vars.prim[[j, S - 3 + i]] = -1.0 * vars.prim[[j, S - 2 - i]];
    }
}

#[cfg(test)]
mod tests {
    use crate::config::{
        meshconfig::{MeshConfig, MeshMode},
        physicsconfig::PhysicsConfig,
    };

    use super::*;
    use approx::assert_relative_eq;
    use ndarray::Array2;
    const E: usize = 3;
    const S: usize = 6;
    const MESHCONFIG: MeshConfig = MeshConfig {
        mode: MeshMode::Cartesian,
        xi_in: 2.0,
        xi_out: 3.0,
        ratio_disk: 1.0,
    };
    const PHYSICSCONFIG: PhysicsConfig = PhysicsConfig {
        adiabatic_index: 1.4,
        units_mode: crate::UnitsMode::SI,
    };

    #[test]
    fn extrapolation_test() {
        // Also test this on log meshes
        let mesh: Mesh<S> = Mesh::new(&MESHCONFIG).unwrap();
        let mut vars: Variables<E, S> = Variables::new(&PHYSICSCONFIG);
        vars.prim.assign(
            &Array2::from_shape_vec(
                (E, S),
                vec![
                    0.0, 0.0, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0, -2.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            )
            .unwrap()
            .view(),
        );
        let prim_expect = Array2::from_shape_vec(
            (E, S),
            vec![
                0.25, 0.25, 0.25, 0.25, 0.25, 0.25, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            ],
        )
        .unwrap();
        for j in 0..E {
            extrapolate_west(j, &mut vars, &mesh);
            extrapolate_east(j, &mut vars, &mesh);
        }
        assert_relative_eq!(vars.prim, prim_expect, max_relative = 1.0e-12);
    }

    #[test]
    fn near_zero_test() {
        let mut vars: Variables<E, S> = Variables::new(&PHYSICSCONFIG);
        vars.prim.assign(
            &Array2::from_shape_vec(
                (E, S),
                vec![
                    0.0, 0.0, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0, -2.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            )
            .unwrap()
            .view(),
        );
        let prim_expect = Array2::from_shape_vec(
            (E, S),
            vec![
                0.25e-10, 0.25e-10, 0.25, 0.25, 0.25e-10, 0.25e-10, -2.0e-10, -2.0e-10, -2.0, -1.5, -1.5e-10, -1.5e-10,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            ],
        )
        .unwrap();
        for j in 0..E {
            near_zero_west(j, &mut vars);
            near_zero_east(j, &mut vars);
        }
        assert_relative_eq!(vars.prim, prim_expect, max_relative = 1.0e-12);
    }

    #[test]
    fn no_gradients_test() {
        let mut vars: Variables<E, S> = Variables::new(&PHYSICSCONFIG);
        vars.prim.assign(
            &Array2::from_shape_vec(
                (E, S),
                vec![
                    0.0, 0.0, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0, -2.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            )
            .unwrap()
            .view(),
        );

        let prim_expect = Array2::from_shape_vec(
            (E, S),
            vec![
                0.25, 0.25, 0.25, 0.25, 0.25, 0.25, -2.0, -2.0, -2.0, -1.5, -1.5, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            ],
        )
        .unwrap();

        for j in 0..E {
            no_gradients_west(j, &mut vars);
            no_gradients_east(j, &mut vars);
        }

        assert_relative_eq!(vars.prim, prim_expect, max_relative = 1.0e-12);
    }

    #[test]
    fn outflow_no_gradients_test() {
        let mut vars: Variables<E, S> = Variables::new(&PHYSICSCONFIG);
        vars.prim.assign(
            &Array2::from_shape_vec(
                (E, S),
                vec![
                    0.0, 0.0, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0, -2.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            )
            .unwrap()
            .view(),
        );
        let prim_expect = Array2::from_shape_vec(
            (E, S),
            vec![
                -0.25, -0.25, 0.25, 0.25, 0.25, 0.25, -2.0, -2.0, -2.0, -1.5, 1.5, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            ],
        )
        .unwrap();
        for j in 0..E {
            outflow_no_gradients_west(j, &mut vars);
            outflow_no_gradients_east(j, &mut vars);
        }
        assert_relative_eq!(vars.prim, prim_expect, max_relative = 1.0e-12);
    }

    #[test]
    fn outflow_extrapolate_test() {
        // Also test this on log meshes
        let mesh: Mesh<S> = Mesh::new(&MESHCONFIG).unwrap();
        let mut vars: Variables<E, S> = Variables::new(&PHYSICSCONFIG);
        vars.prim.assign(
            &Array2::from_shape_vec(
                (E, S),
                vec![
                    0.0, 0.0, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0, -2.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            )
            .unwrap()
            .view(),
        );
        let prim_expect = Array2::from_shape_vec(
            (E, S),
            vec![
                -0.25, -0.25, 0.25, 0.25, 0.25, 0.25, -3.0, -2.5, -2.0, -1.5, 1.5, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            ],
        )
        .unwrap();
        for j in 0..E {
            outflow_extrapolate_west(j, &mut vars, &mesh);
            outflow_extrapolate_east(j, &mut vars, &mesh);
        }
        assert_relative_eq!(vars.prim, prim_expect, max_relative = 1.0e-12);
    }

    #[test]
    fn reflecting_test() {
        // Also test this on log meshes
        let mut vars: Variables<E, S> = Variables::new(&PHYSICSCONFIG);
        vars.prim.assign(
            &Array2::from_shape_vec(
                (E, S),
                vec![
                    0.0, 0.0, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0, -2.0, -1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            )
            .unwrap()
            .view(),
        );
        let prim_expect = Array2::from_shape_vec(
            (E, S),
            vec![
                -0.25, -0.25, 0.25, 0.25, -0.25, -0.25, 1.5, 2.0, -2.0, -1.5, 1.5, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            ],
        )
        .unwrap();
        for j in 0..E {
            reflecting_west(j, &mut vars);
            reflecting_east(j, &mut vars);
        }
        assert_relative_eq!(vars.prim, prim_expect, max_relative = 1.0e-12);
    }
}
