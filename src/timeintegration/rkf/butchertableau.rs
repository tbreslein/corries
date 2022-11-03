// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use ndarray::{ArrayD, IxDyn};

use crate::config::numericsconfig::{NumericsConfig, RKFMode};

pub struct ButcherTableau {
    pub order: usize,
    pub a: ArrayD<f64>,
    pub b_high: ArrayD<f64>,
    pub b_low: ArrayD<f64>,
    pub c: ArrayD<f64>,
    pub asc: bool,
}

#[rustfmt::skip]
impl ButcherTableau {
    pub fn new(numericsconfig: &NumericsConfig) -> Self {
        let order = get_order(numericsconfig.rkf_mode);
        return match numericsconfig.rkf_mode {
            RKFMode::RK1 => Self {
                order,
                a: ArrayD::from_shape_vec(
                    IxDyn(&[order, order]),
                    vec![0.0]
                ).unwrap(),
                b_high: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![1.0]
                ).unwrap(),
                b_low: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0]
                ).unwrap(),
                c: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0]
                ).unwrap(),
                asc: numericsconfig.asc,
            },
            RKFMode::RK2 => Self {
                order,
                a: ArrayD::from_shape_vec(
                    IxDyn(&[order, order]),
                    vec![0.0, 0.0,
                         0.5, 0.0]
                ).unwrap(),
                b_high: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0, 1.0]
                ).unwrap(),
                b_low: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0, 0.0]
                ).unwrap(),
                c: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0, 0.5]
                ).unwrap(),
                asc: numericsconfig.asc,
            },
            RKFMode::RK3 => Self {
                order,
                a: ArrayD::from_shape_vec(
                    IxDyn(&[order, order]),
                    vec![0.0, 0.0, 0.0,
                         0.5, 0.0, 0.0,
                         -1.0, 2.0, 0.0]
                ).unwrap(),
                b_high: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![1.0/6.0, 2.0/3.0, 1.0/6.0]
                ).unwrap(),
                b_low: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0, 1.0, 0.0]
                ).unwrap(),
                c: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0, 0.5, 1.0]
                ).unwrap(),
                asc: numericsconfig.asc,
            },
            RKFMode::RK4 => Self {
                order,
                a: ArrayD::from_shape_vec(
                    IxDyn(&[order, order]),
                    vec![0.0, 0.0, 0.0, 0.0,
                         0.5, 0.0, 0.0, 0.0,
                         0.0, 0.5, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0]
                ).unwrap(),
                b_high: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0]
                ).unwrap(),
                b_low: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0, 0.0, 0.0, 0.0]
                ).unwrap(),
                c: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0, 0.5, 0.5, 1.0]
                ).unwrap(),
                asc: numericsconfig.asc,
            },
            RKFMode::Heun2 => Self {
                order,
                a: ArrayD::from_shape_vec(
                    IxDyn(&[order, order]),
                    vec![0.0, 0.0,
                         1.0, 0.0]
                ).unwrap(),
                b_high: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.5, 0.5]
                ).unwrap(),
                b_low: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![1.0, 0.0]
                ).unwrap(),
                c: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0, 1.0]
                ).unwrap(),
                asc: numericsconfig.asc,
            },
            RKFMode::RKF12 => Self {
                order,
                a: ArrayD::from_shape_vec(
                    IxDyn(&[order, order]),
                    vec![0.0,       0.0,         0.0,
                         0.5,       0.0,         0.0,
                         1.0/256.0, 255.0/256.0, 0.0]
                ).unwrap(),
                b_high: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![1.0/512.0, 255.0/256.0, 1.0/512.0]
                ).unwrap(),
                b_low: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![1.0/256.0, 255.0/256.0, 0.0]
                ).unwrap(),
                c: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0, 0.5, 1.0]
                ).unwrap(),
                asc: numericsconfig.asc,
            },
            RKFMode::RKF45 => Self {
                order,
                a: ArrayD::from_shape_vec(
                    IxDyn(&[order, order]),
                    vec![0.0,            0.0,           0.0,           0.0,            0.0,       0.0,
                         0.25,           0.0,           0.0,           0.0,            0.0,       0.0,
                         3.0/32.0,       9.0/32.0,      0.0,           0.0,            0.0,       0.0,
                         1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0,            0.0,       0.0,
                         439.0/216.0,   -8.0,           3680.0/513.0, -845.0/4104.0,   0.0,       0.0,
                         -8.0/27.0,      2.0,          -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0, 0.0]
                ).unwrap(),
                b_high: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0]
                ).unwrap(),
                b_low: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -0.2, 0.0]
                ).unwrap(),
                c: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0, 0.25, 0.375, 12.0/13.0, 1.0, 0.5]
                ).unwrap(),
                asc: numericsconfig.asc,
            },
            RKFMode::SSPRK3 => Self {
                order,
                a: ArrayD::from_shape_vec(
                    IxDyn(&[order, order]),
                    vec![0.0,  0.0,  0.0,
                         1.0,  0.0,  0.0,
                         0.25, 0.25, 0.0]
                ).unwrap(),
                b_high: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![1.0/6.0, 1.0/6.0, 2.0/3.0]
                ).unwrap(),
                b_low: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.5, 0.5, 0.0]
                ).unwrap(),
                c: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0,  1.0, 0.5]
                ).unwrap(),
                asc: numericsconfig.asc,
            },
            RKFMode::SSPRK5 => Self {
                order,
                a: ArrayD::from_shape_vec(
                    IxDyn(&[order, order]),
                    vec![0.0,     0.0,     0.0,     0.0,     0.0,
                         0.36717, 0.0,     0.0,     0.0,     0.0,
                         0.26802, 0.31720, 0.0,     0.0,     0.0,
                         0.11606, 0.13735, 0.18816, 0.0,     0.0,
                         0.11212, 0.13269, 0.18178, 0.41980, 0.0,
                    ]
                ).unwrap(),
                b_high: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.17279, 0.094505, 0.12947, 0.29899, 0.30424]
                ).unwrap(),
                b_low: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.12293, 0.31981, -0.15316, 0.31887, 0.39155]
                ).unwrap(),
                c: ArrayD::from_shape_vec(
                    IxDyn(&[order]),
                    vec![0.0, 0.36717, 0.58522, 0.44156, 0.8464]
                ).unwrap(),
                asc: numericsconfig.asc,
            },
        };
    }
    
}

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
