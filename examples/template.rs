// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::Context, Result};
use corries::components::Runner;
use corries::prelude::*;

// =============
// CONFIGURATION
// =============

const S: usize = 10;
set_Physics_and_E!(Euler1DIsot);
type N = Hll<E, S>;
type T = RungeKuttaFehlberg<P, E, S>;

fn main() -> Result<()> {
    CorriesConfig {
        print_banner: true,
        mesh_config: MeshConfig {
            mode: MeshMode::Cartesian,
            xi_in: 1.0,
            xi_out: 2.0,
        },
        physics_config: PhysicsConfig {
            adiabatic_index: 1.4,
            units_mode: UnitsMode::SI,
        },
        boundary_condition_west: BoundaryMode::Custom(vec![
            (0, CustomBoundaryMode::NoGradients),
            (1, CustomBoundaryMode::NoGradients),
        ]),
        boundary_condition_east: BoundaryMode::Custom(vec![
            (0, CustomBoundaryMode::NoGradients),
            (1, CustomBoundaryMode::NoGradients),
        ]),
        numerics_config: NumericsConfig {
            numflux_config: NumFluxConfig::Kt {
                limiter_mode: LimiterMode::Monocent(1.2),
            },
            time_integration_config: TimeIntegrationConfig::Rkf(RkfConfig {
                rkf_mode: RKFMode::RK4,
                asc: false,
                asc_relative_tolerance: 0.001,
                asc_absolute_tolerance: 0.1,
                asc_timestep_friction: 0.08,
            }),
            iter_max: usize::MAX - 2,
            t0: 0.0,
            t_end: 10.0,
            dt_min: 1.0e-12,
            dt_max: f64::MAX,
            dt_cfl_param: 0.1,
        },
        output_counter_max: 10,
        writer_config: vec![
            OutputConfig {
                stream_mode: StreamMode::Stdout,
                formatting_mode: FormattingMode::TSV,
                string_conversion_mode: ToStringConversionMode::Scalar,
                folder_name: "".to_string(),
                should_clear_out_folder: false,
                file_name: "".to_string(),
                precision: 3,
                should_print_ghostcells: false,
                should_print_metadata: true,
                data_names: vec![DataName::Iter, DataName::T, DataName::Dt, DataName::DtKind],
            },
            OutputConfig {
                stream_mode: StreamMode::File,
                formatting_mode: FormattingMode::TSV,
                string_conversion_mode: ToStringConversionMode::Vector,
                folder_name: "results/accretiondisk".to_string(),
                should_clear_out_folder: true,
                file_name: "accretiondisk".to_string(),
                precision: 7,
                should_print_ghostcells: true,
                should_print_metadata: false,
                data_names: vec![DataName::XiCent],
            },
        ],
    }
    .init_corries::<P, N, T, E, S>(|u, _, _| {
        let breakpoint_index = (S as f64 * 0.5) as usize;
        for i in 0..breakpoint_index {
            u.cent.prim[[0, i]] = 1.0;
            u.cent.prim[[E - 1, i]] = 1.0;
        }
        for i in breakpoint_index..S {
            u.cent.prim[[0, i]] = 0.125;
            u.cent.prim[[E - 1, i]] = 0.1;
        }
        Ok(())
    })
    .context("While calling CorriesConfig::init_corries")?
    .run_corries()
}
