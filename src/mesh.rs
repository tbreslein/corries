// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use crate::config::meshconfig::{MeshConfig, MeshMode};
use color_eyre::Result;
use ndarray::{ArrayD, IxDyn};

/// The Mesh the simulation runs on.
///
/// This struct contains general information about the current mesh, like the number of cells,
/// and special indexes like the first non-ghost-cell cell index, but also vectors containing the
/// distance between each cell and the origin, or the geometric scales per cell.
///
/// Note that we use generalised coordinates xi, Phi, eta for the coordinate names. For example, in
/// cartesian coordinates, you would then usually name these xi = x, Phi = y, eta = z, or in
/// cylindrical coordinates you would name them xi = r (the radius), Phi = phi (the polar angle),
/// eta = z.
pub struct Mesh {
    /// The kind of mesh this is, i.e. cartesian, log-cylindrical, etc.
    pub mode: MeshMode,

    /// Whether this mesh is logarithmic
    pub is_logarithmic: bool,

    /// Whether this mesh is linear
    pub is_linear: bool,

    /// The number of cells in the computational area
    pub n_comp: usize,

    /// The number of ghost cells per mesh edge
    pub n_gc: usize,

    /// The number of all cells summed up
    pub n_all: usize,

    /// The first cell index in the mesh
    pub imin: usize,

    /// The last cell index in the mesh
    pub imax: usize,

    /// The first cell index after the western ghost cells
    pub ixi_in: usize,

    /// The last cell index before the eastern ghost cells
    pub ixi_out: usize,

    /// The radial coordinate at ixi_in
    pub xi_in: f64,

    /// The radial coordinate at ixi_out
    pub xi_out: f64,

    /// The ratio of radial distance to the origin between the massive disk and and computational
    /// area
    pub ratio_disk: f64,

    /// The number of cells in the massive disk
    pub n_disk: usize,

    /// The cell index of the break point between the massive disk and the rest and the
    /// computational area
    pub idisk_end: usize,

    /// The differential along xi
    pub dxi: f64,

    /// The differential along eta
    pub deta: f64,

    /// The differential along Phi
    pub dphi: f64,

    /// xi coordinates at the center of each cell
    pub xi_cent: ArrayD<f64>,

    /// xi coordinates at the west border of each cell
    pub xi_west: ArrayD<f64>,

    /// xi coordinates at the east border of each cell
    pub xi_east: ArrayD<f64>,

    /// Inverse of xi_cent
    pub xi_cent_inv: ArrayD<f64>,

    /// Inverse of xi_west
    pub xi_west_inv: ArrayD<f64>,

    /// Inverse of xi_east
    pub xi_east_inv: ArrayD<f64>,

    /// Geometric scale along xi at the centre of each cell
    pub h_xi_cent: ArrayD<f64>,

    /// Geometric scale along xi at the west border of each cell
    pub h_xi_west: ArrayD<f64>,

    /// Geometric scale along xi at the east border of each cell
    pub h_xi_east: ArrayD<f64>,

    /// Geometric scale along eta at the centre of each cell
    pub h_eta_cent: ArrayD<f64>,

    /// Geometric scale along eta at the west border of each cell
    pub h_eta_west: ArrayD<f64>,

    /// Geometric scale along eta at the east border of each cell
    pub h_eta_east: ArrayD<f64>,

    /// Geometric scale along Phi at the centre of each cell
    pub h_phi_cent: ArrayD<f64>,

    /// Geometric scale along Phi at the west border of each cell
    pub h_phi_west: ArrayD<f64>,

    /// Geometric scale along Phi at the east border of each cell
    pub h_phi_east: ArrayD<f64>,

    /// Square root of the scalar metric
    pub sqrt_g: ArrayD<f64>,

    /// Line element along xi
    pub line_xi: ArrayD<f64>,

    /// Inverse of line_xi
    pub line_xi_inv: ArrayD<f64>,

    /// Shorthand: area_west / (deta / dphi)
    pub d_area_xi_deta_dphi_west: ArrayD<f64>,

    /// Shorthand: area_east / (deta / dphi)
    pub d_area_xi_deta_dphi_east: ArrayD<f64>,

    /// Surface area of each cell normal to xi x Phi
    pub area_cell: ArrayD<f64>,

    /// Surface area of each cell's west border normal to Phi x eta
    pub area_west: ArrayD<f64>,

    /// Surface area of each cell's east border normal to Phi x eta
    pub area_east: ArrayD<f64>,

    /// Volume of each cell
    pub volume: ArrayD<f64>,

    /// Shorthand: deta * dPhi / dVolume
    pub deta_dphi_d_volume: ArrayD<f64>,

    /// Distance between west and east border of each cell
    pub cell_width: ArrayD<f64>,

    /// Inverse of cell_width
    pub cell_width_inv: ArrayD<f64>,

    /// Commutator coefficient for eta and xi
    pub cexe: ArrayD<f64>,

    /// Commutator coefficient for Phi and xi
    pub cpxp: ArrayD<f64>,

    /// Commutator coefficient for xi and eta
    pub cxex: ArrayD<f64>,

    /// Commutator coefficient for Phi and eta
    pub cpep: ArrayD<f64>,

    /// Commutator coefficient for xi and Phi
    pub cxpx: ArrayD<f64>,

    /// Commutator coefficient for eta and Phi
    pub cepe: ArrayD<f64>,

    /// Shorthand: -1.0 * cexe
    pub minus_cexe: ArrayD<f64>,

    /// Shorthand: -1.0 * cpxp
    pub minus_cpxp: ArrayD<f64>,

    /// Shorthand: -1.0 * cxex
    pub minus_cxex: ArrayD<f64>,

    /// Shorthand: -1.0 * cpep
    pub minus_cpep: ArrayD<f64>,

    /// Shorthand: -1.0 * cxpx
    pub minus_cxpx: ArrayD<f64>,

    /// Shorthand: -1.0 * cepe
    pub minus_cepe: ArrayD<f64>,
}

impl Mesh {
    pub fn new(meshconf: &MeshConfig) -> Result<Mesh> {
        let mode = meshconf.mode;
        let is_logarithmic = match mode {
            MeshMode::Cartesian => false,
        };
        let is_linear = !is_logarithmic;

        let n_comp = meshconf.n_comp;
        let n_gc = meshconf.n_gc;
        let n_all = n_comp + 2 * n_gc;

        let imin: usize = 0;
        let imax = imin + n_all - 1;
        let ixi_in = imin + n_gc;
        let ixi_out = imax - n_gc;
        let xi_in = meshconf.xi_in;
        let xi_out = meshconf.xi_out;

        let ratio_disk = meshconf.ratio_disk;
        let n_disk = (n_comp as f64 * ratio_disk) as usize;
        let idisk_end = n_gc + n_disk;

        let dxi = match mode {
            MeshMode::Cartesian => (xi_out - xi_in) / n_comp as f64,
        };
        let dphi = match mode {
            MeshMode::Cartesian => 1.0,
        };
        let deta = match mode {
            MeshMode::Cartesian => 1.0,
        };

        let xi_cent = {
            let mut v: Vec<f64> = Vec::with_capacity(n_all);
            let mut f = |x: f64| {
                for i in imin..=imax {
                    v.push(x + (0.5 * (2.0 * (i as f64 - 1.0) - 1.0) * dxi));
                }
            };
            match mode {
                MeshMode::Cartesian => f(xi_in),
            };
            ndarray::ArrayD::from_shape_vec(IxDyn(&[v.len()]), v)?
        };
        let xi_west = match mode {
            MeshMode::Cartesian => &xi_cent - 0.5 * dxi,
        };
        let xi_east = match mode {
            MeshMode::Cartesian => &xi_cent + 0.5 * dxi,
        };

        let xi_cent_inv = 1.0 / &xi_cent;
        let xi_west_inv = 1.0 / &xi_west;
        let xi_east_inv = 1.0 / &xi_east;

        let h_xi_cent = match mode {
            MeshMode::Cartesian => &xi_cent / &xi_cent,
        };
        let h_xi_west = match mode {
            MeshMode::Cartesian => &xi_west / &xi_west,
        };
        let h_xi_east = match mode {
            MeshMode::Cartesian => &xi_east / &xi_east,
        };

        let h_eta_cent = match mode {
            MeshMode::Cartesian => &xi_cent / &xi_cent,
        };
        let h_eta_west = match mode {
            MeshMode::Cartesian => &xi_west / &xi_west,
        };
        let h_eta_east = match mode {
            MeshMode::Cartesian => &xi_east / &xi_east,
        };

        let h_phi_cent = match mode {
            MeshMode::Cartesian => &xi_cent / &xi_cent,
        };
        let h_phi_west = match mode {
            MeshMode::Cartesian => &xi_west / &xi_west,
        };
        let h_phi_east = match mode {
            MeshMode::Cartesian => &xi_east / &xi_east,
        };

        let sqrt_g = &h_xi_cent * &h_eta_cent * &h_phi_cent;

        let line_xi = &h_xi_cent * dxi;
        let line_xi_inv = 1.0 / &line_xi;

        let d_area_xi_deta_dphi_west = &h_eta_west * &h_phi_west;
        let d_area_xi_deta_dphi_east = &h_eta_east * &h_phi_east;

        let area_cell = match mode {
            MeshMode::Cartesian => &xi_east * &xi_east - &xi_west * &xi_west,
        };
        let area_west = &d_area_xi_deta_dphi_west * deta * dphi;
        let area_east = &d_area_xi_deta_dphi_east * deta * dphi;
        let volume = &sqrt_g * dxi * deta * dphi;
        let deta_dphi_d_volume = deta * dphi / (f64::EPSILON + &volume);

        let cell_width = &xi_east - &xi_west;
        let cell_width_inv = 1.0 / &cell_width;

        let cexe =
            0.5 * (&h_phi_east + &h_phi_west) * (&h_eta_east - &h_eta_west) * &deta_dphi_d_volume;
        let cpxp =
            0.5 * (&h_eta_east + &h_eta_west) * (&h_phi_east - &h_phi_west) * &deta_dphi_d_volume;

        let cxex = &xi_cent / &xi_cent;
        let cpep = &xi_cent / &xi_cent;
        let cxpx = &xi_cent / &xi_cent;
        let cepe = &xi_cent / &xi_cent;

        let minus_cexe = -1.0 * &cexe;
        let minus_cpxp = -1.0 * &cpxp;
        let minus_cxex = -1.0 * &cxex;
        let minus_cpep = -1.0 * &cpep;
        let minus_cxpx = -1.0 * &cxpx;
        let minus_cepe = -1.0 * &cepe;

        return Ok(Mesh {
            mode,
            is_logarithmic,
            is_linear,
            n_comp,
            n_gc,
            n_all,
            imin,
            imax,
            ixi_in,
            ixi_out,
            xi_in,
            xi_out,
            ratio_disk,
            n_disk,
            idisk_end,
            dxi,
            deta,
            dphi,
            xi_cent,
            xi_west,
            xi_east,
            xi_cent_inv,
            xi_west_inv,
            xi_east_inv,
            h_xi_cent,
            h_xi_west,
            h_xi_east,
            h_eta_cent,
            h_eta_west,
            h_eta_east,
            h_phi_cent,
            h_phi_west,
            h_phi_east,
            sqrt_g,
            line_xi,
            line_xi_inv,
            d_area_xi_deta_dphi_west,
            d_area_xi_deta_dphi_east,
            area_cell,
            area_west,
            area_east,
            volume,
            deta_dphi_d_volume,
            cell_width,
            cell_width_inv,
            cexe,
            cpxp,
            cxex,
            cpep,
            cxpx,
            cepe,
            minus_cexe,
            minus_cpxp,
            minus_cxex,
            minus_cpep,
            minus_cxpx,
            minus_cepe,
        });
    }
}
