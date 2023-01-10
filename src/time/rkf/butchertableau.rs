// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [ButcherTableau] struct

use ndarray::{Array1, Array2};

use crate::config::numericsconfig::{RKFMode, RkfConfig};

/// Carries the butcher tableau needed for implement different Runge-Kutta-Fehlberg methods
pub struct ButcherTableau {
    /// Order of the method
    pub order: usize,

    /// Coefficient matrix a
    pub a: Array2<f64>,

    /// Coefficient vector b for high-order solutions
    pub b_high: Array1<f64>,

    /// Coefficient vector b for low-order solutions
    pub b_low: Array1<f64>,

    /// Coefficient vector c
    pub c: Array1<f64>,

    /// Whether to use automated step control
    pub asc: bool,
}

#[rustfmt::skip]
impl ButcherTableau {
    /// Constructs a new [ButcherTableau] object.
    ///
    /// # Argument
    ///
    /// * `rkfconfig` - RKF specific configuration
    pub fn new(rkfconfig: &RkfConfig) -> Self {
        let order = get_order(rkfconfig.rkf_mode);
        let asc = rkfconfig.asc;
        return match rkfconfig.rkf_mode {
            RKFMode::RK1 => Self {
                order,
                a: Array2::from_shape_vec(
                    (order, order),
                    vec![0.0]
                ).unwrap(),
                b_high: Array1::from_shape_vec(
                    order,
                    vec![1.0]
                ).unwrap(),
                b_low: Array1::from_shape_vec(
                    order,
                    vec![0.0]
                ).unwrap(),
                c: Array1::from_shape_vec(
                    order,
                    vec![0.0]
                ).unwrap(),
                asc: if asc {
                    println!("WARNING: You turned on automated step control (asc) for {:?}, but this RKF method does not support asc! This option is ignored and asc will be set to false.", rkfconfig.rkf_mode);
                    false
                } else {
                    asc
                }
            },
            RKFMode::RK2 => Self {
                order,
                a: Array2::from_shape_vec(
                    (order, order),
                    vec![0.0, 0.0,
                         0.5, 0.0]
                ).unwrap(),
                b_high: Array1::from_shape_vec(
                    order,
                    vec![0.0, 1.0]
                ).unwrap(),
                b_low: Array1::from_shape_vec(
                    order,
                    vec![0.0, 0.0]
                ).unwrap(),
                c: Array1::from_shape_vec(
                    order,
                    vec![0.0, 0.5]
                ).unwrap(),
                asc: if asc {
                    println!("WARNING: You turned on automated step control (asc) for {:?}, but this RKF method does not support asc! This option is ignored and asc will be set to false.", rkfconfig.rkf_mode);
                    false
                } else {
                    asc
                }
            },
            RKFMode::RK3 => Self {
                order,
                a: Array2::from_shape_vec(
                    (order, order),
                    vec![0.0, 0.0, 0.0,
                         0.5, 0.0, 0.0,
                         -1.0, 2.0, 0.0]
                ).unwrap(),
                b_high: Array1::from_shape_vec(
                    order,
                    vec![1.0/6.0, 2.0/3.0, 1.0/6.0]
                ).unwrap(),
                b_low: Array1::from_shape_vec(
                    order,
                    vec![0.0, 1.0, 0.0]
                ).unwrap(),
                c: Array1::from_shape_vec(
                    order,
                    vec![0.0, 0.5, 1.0]
                ).unwrap(),
                asc
            },
            RKFMode::RK4 => Self {
                order,
                a: Array2::from_shape_vec(
                    (order, order),
                    vec![0.0, 0.0, 0.0, 0.0,
                         0.5, 0.0, 0.0, 0.0,
                         0.0, 0.5, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0]
                ).unwrap(),
                b_high: Array1::from_shape_vec(
                    order,
                    vec![1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0]
                ).unwrap(),
                b_low: Array1::from_shape_vec(
                    order,
                    vec![0.0, 0.0, 0.0, 0.0]
                ).unwrap(),
                c: Array1::from_shape_vec(
                    order,
                    vec![0.0, 0.5, 0.5, 1.0]
                ).unwrap(),
                asc: if asc {
                    println!("WARNING: You turned on automated step control (asc) for {:?}, but this RKF method does not support asc! This option is ignored and asc will be set to false.", rkfconfig.rkf_mode);
                    false
                } else {
                    asc
                }
            },
            RKFMode::Heun2 => Self {
                order,
                a: Array2::from_shape_vec(
                    (order, order),
                    vec![0.0, 0.0,
                         1.0, 0.0]
                ).unwrap(),
                b_high: Array1::from_shape_vec(
                    order,
                    vec![0.5, 0.5]
                ).unwrap(),
                b_low: Array1::from_shape_vec(
                    order,
                    vec![1.0, 0.0]
                ).unwrap(),
                c: Array1::from_shape_vec(
                    order,
                    vec![0.0, 1.0]
                ).unwrap(),
                asc,
            },
            RKFMode::RKF12 => Self {
                order,
                a: Array2::from_shape_vec(
                    (order, order),
                    vec![0.0,       0.0,         0.0,
                         0.5,       0.0,         0.0,
                         1.0/256.0, 255.0/256.0, 0.0]
                ).unwrap(),
                b_high: Array1::from_shape_vec(
                    order,
                    vec![1.0/512.0, 255.0/256.0, 1.0/512.0]
                ).unwrap(),
                b_low: Array1::from_shape_vec(
                    order,
                    vec![1.0/256.0, 255.0/256.0, 0.0]
                ).unwrap(),
                c: Array1::from_shape_vec(
                    order,
                    vec![0.0, 0.5, 1.0]
                ).unwrap(),
                asc
            },
            RKFMode::RKF45 => Self {
                order,
                a: Array2::from_shape_vec(
                    (order, order),
                    vec![0.0,            0.0,           0.0,           0.0,            0.0,       0.0,
                         0.25,           0.0,           0.0,           0.0,            0.0,       0.0,
                         3.0/32.0,       9.0/32.0,      0.0,           0.0,            0.0,       0.0,
                         1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0,            0.0,       0.0,
                         439.0/216.0,   -8.0,           3680.0/513.0, -845.0/4104.0,   0.0,       0.0,
                         -8.0/27.0,      2.0,          -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0, 0.0]
                ).unwrap(),
                b_high: Array1::from_shape_vec(
                    order,
                    vec![16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0]
                ).unwrap(),
                b_low: Array1::from_shape_vec(
                    order,
                    vec![25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -0.2, 0.0]
                ).unwrap(),
                c: Array1::from_shape_vec(
                    order,
                    vec![0.0, 0.25, 0.375, 12.0/13.0, 1.0, 0.5]
                ).unwrap(),
                asc
            },
            RKFMode::SSPRK3 => Self {
                order,
                a: Array2::from_shape_vec(
                    (order, order),
                    vec![0.0,  0.0,  0.0,
                         1.0,  0.0,  0.0,
                         0.25, 0.25, 0.0]
                ).unwrap(),
                b_high: Array1::from_shape_vec(
                    order,
                    vec![1.0/6.0, 1.0/6.0, 2.0/3.0]
                ).unwrap(),
                b_low: Array1::from_shape_vec(
                    order,
                    vec![0.5, 0.5, 0.0]
                ).unwrap(),
                c: Array1::from_shape_vec(
                    order,
                    vec![0.0,  1.0, 0.5]
                ).unwrap(),
                asc
            },
            RKFMode::SSPRK5 => Self {
                order,
                a: Array2::from_shape_vec(
                    (order, order),
                    vec![0.0,     0.0,     0.0,     0.0,     0.0,
                         0.36717, 0.0,     0.0,     0.0,     0.0,
                         0.26802, 0.31720, 0.0,     0.0,     0.0,
                         0.11606, 0.13735, 0.18816, 0.0,     0.0,
                         0.11212, 0.13269, 0.18178, 0.41980, 0.0,
                    ]
                ).unwrap(),
                b_high: Array1::from_shape_vec(
                    order,
                    vec![0.17279, 0.094505, 0.12947, 0.29899, 0.30424]
                ).unwrap(),
                b_low: Array1::from_shape_vec(
                    order,
                    vec![0.12293, 0.31981, -0.15316, 0.31887, 0.39155]
                ).unwrap(),
                c: Array1::from_shape_vec(
                    order,
                    vec![0.0, 0.36717, 0.58522, 0.44156, 0.8464]
                ).unwrap(),
                asc
            },
        };
    }
    
}

/// Takes an [RKFMode] and returns the order of the method.
fn get_order(mode: RKFMode) -> usize {
    return match mode {
        RKFMode::RK1 => 1,
        RKFMode::RK2 => 2,
        RKFMode::RK3 => 3,
        RKFMode::RK4 => 4,
        RKFMode::Heun2 => 2,
        RKFMode::RKF12 => 3,
        RKFMode::RKF45 => 6,
        RKFMode::SSPRK3 => 3,
        RKFMode::SSPRK5 => 5,
    };
}
