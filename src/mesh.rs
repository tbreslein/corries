// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Mesh] struct that provides coordinates and geometric information

use crate::config::meshconfig::{MeshConfig, MeshMode};
use crate::data::{Data, DataName, StructAssociation};
use crate::errorhandling::Validation;
use crate::Collectable;
use color_eyre::eyre::{bail, ensure, Context};
use color_eyre::Result;
use ndarray::Array1;

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
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Mesh<const S: usize> {
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
    pub xi_cent: Array1<f64>,

    /// xi coordinates at the west border of each cell
    pub xi_west: Array1<f64>,

    /// xi coordinates at the east border of each cell
    pub xi_east: Array1<f64>,

    /// Inverse of xi_cent
    pub xi_cent_inv: Array1<f64>,

    /// Inverse of xi_west
    pub xi_west_inv: Array1<f64>,

    /// Inverse of xi_east
    pub xi_east_inv: Array1<f64>,

    /// Geometric scale along xi at the centre of each cell
    pub h_xi_cent: Array1<f64>,

    /// Geometric scale along xi at the west border of each cell
    pub h_xi_west: Array1<f64>,

    /// Geometric scale along xi at the east border of each cell
    pub h_xi_east: Array1<f64>,

    /// Geometric scale along eta at the centre of each cell
    pub h_eta_cent: Array1<f64>,

    /// Geometric scale along eta at the west border of each cell
    pub h_eta_west: Array1<f64>,

    /// Geometric scale along eta at the east border of each cell
    pub h_eta_east: Array1<f64>,

    /// Geometric scale along Phi at the centre of each cell
    pub h_phi_cent: Array1<f64>,

    /// Geometric scale along Phi at the west border of each cell
    pub h_phi_west: Array1<f64>,

    /// Geometric scale along Phi at the east border of each cell
    pub h_phi_east: Array1<f64>,

    /// Square root of the scalar metric
    pub sqrt_g: Array1<f64>,

    /// Line element along xi
    pub line_xi: Array1<f64>,

    /// Inverse of line_xi
    pub line_xi_inv: Array1<f64>,

    /// Shorthand: area_west / (deta / dphi)
    pub d_area_xi_deta_dphi_west: Array1<f64>,

    /// Shorthand: area_east / (deta / dphi)
    pub d_area_xi_deta_dphi_east: Array1<f64>,

    /// Surface area of each cell normal to xi x Phi
    pub area_cell: Array1<f64>,

    /// Surface area of each cell's west border normal to Phi x eta
    pub area_west: Array1<f64>,

    /// Surface area of each cell's east border normal to Phi x eta
    pub area_east: Array1<f64>,

    /// Volume of each cell
    pub volume: Array1<f64>,

    /// Shorthand: deta * dPhi / dVolume
    pub deta_dphi_d_volume: Array1<f64>,

    /// Distance between west and east border of each cell
    pub cell_width: Array1<f64>,

    /// Inverse of cell_width
    pub cell_width_inv: Array1<f64>,

    /// Commutator coefficient for eta and xi
    pub cexe: Array1<f64>,

    /// Commutator coefficient for Phi and xi
    pub cpxp: Array1<f64>,

    /// Commutator coefficient for xi and eta
    pub cxex: Array1<f64>,

    /// Commutator coefficient for Phi and eta
    pub cpep: Array1<f64>,

    /// Commutator coefficient for xi and Phi
    pub cxpx: Array1<f64>,

    /// Commutator coefficient for eta and Phi
    pub cepe: Array1<f64>,

    /// Shorthand: -1.0 * cexe
    pub minus_cexe: Array1<f64>,

    /// Shorthand: -1.0 * cpxp
    pub minus_cpxp: Array1<f64>,

    /// Shorthand: -1.0 * cxex
    pub minus_cxex: Array1<f64>,

    /// Shorthand: -1.0 * cpep
    pub minus_cpep: Array1<f64>,

    /// Shorthand: -1.0 * cxpx
    pub minus_cxpx: Array1<f64>,

    /// Shorthand: -1.0 * cepe
    pub minus_cepe: Array1<f64>,
}

unsafe impl<const S: usize> Send for Mesh<S> {}
unsafe impl<const S: usize> Sync for Mesh<S> {}

impl<const S: usize> Mesh<S> {
    /// Builds a new `color_eyre::Result<Mesh>` object.
    ///
    /// # Arguments
    ///
    /// * `meshconf` - `MeshConfig` containing configuration to build `Mesh` objects
    pub fn new(meshconf: &MeshConfig) -> Result<Mesh<S>> {
        let mode = meshconf.mode;
        let is_logarithmic = match mode {
            MeshMode::Cartesian => false,
        };
        let is_linear = !is_logarithmic;

        let n_all = S;
        let n_gc = 2;
        let n_comp = S - 2 * n_gc;

        let imin = 0;
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

        let xi_west = {
            let mut v = Array1::zeros(S);
            let mut f = |x: f64| {
                for i in imin..=imax {
                    v[i] = x + dxi * (i as i64 - ixi_in as i64) as f64;
                }
            };
            match mode {
                MeshMode::Cartesian => f(xi_in),
            };
            v
        };
        let xi_cent = match mode {
            MeshMode::Cartesian => xi_west.mapv(|xi| xi + 0.5 * dxi),
        };
        let xi_east = match mode {
            MeshMode::Cartesian => xi_west.mapv(|xi| xi + dxi),
        };

        let xi_cent_inv = xi_cent.mapv(f64::recip);
        let xi_west_inv = xi_west.mapv(f64::recip);
        let xi_east_inv = xi_east.mapv(f64::recip);

        let h_xi_cent = match mode {
            MeshMode::Cartesian => xi_cent.mapv(|_| 1.0),
        };
        let h_xi_west = match mode {
            MeshMode::Cartesian => xi_west.mapv(|_| 1.0),
        };
        let h_xi_east = match mode {
            MeshMode::Cartesian => xi_east.mapv(|_| 1.0),
        };

        let h_eta_cent = match mode {
            MeshMode::Cartesian => xi_cent.mapv(|_| 1.0),
        };
        let h_eta_west = match mode {
            MeshMode::Cartesian => xi_west.mapv(|_| 1.0),
        };
        let h_eta_east = match mode {
            MeshMode::Cartesian => xi_east.mapv(|_| 1.0),
        };

        let h_phi_cent = match mode {
            MeshMode::Cartesian => xi_cent.mapv(|_| 1.0),
        };
        let h_phi_west = match mode {
            MeshMode::Cartesian => xi_west.mapv(|_| 1.0),
        };
        let h_phi_east = match mode {
            MeshMode::Cartesian => xi_east.mapv(|_| 1.0),
        };

        let sqrt_g = &h_xi_cent * &h_eta_cent * &h_phi_cent;

        let line_xi = &h_xi_cent * dxi;
        let line_xi_inv = line_xi.mapv(f64::recip);

        let d_area_xi_deta_dphi_west = &h_eta_west * &h_phi_west;
        let d_area_xi_deta_dphi_east = &h_eta_east * &h_phi_east;

        let area_cell = match mode {
            MeshMode::Cartesian => (&xi_east - &xi_west) * (&xi_east - &xi_west),
        }
        .map(|a| a.abs());

        let area_west = &d_area_xi_deta_dphi_west * deta * dphi;
        let area_east = &d_area_xi_deta_dphi_east * deta * dphi;
        let volume = &sqrt_g * dxi * deta * dphi;
        let deta_dphi_d_volume = deta * dphi / (f64::EPSILON + &volume);

        let cell_width = &xi_east - &xi_west;
        let cell_width_inv = cell_width.mapv(f64::recip);

        let cexe = 0.5 * (&h_phi_east + &h_phi_west) * (&h_eta_east - &h_eta_west) * &deta_dphi_d_volume;
        let cpxp = 0.5 * (&h_eta_east + &h_eta_west) * (&h_phi_east - &h_phi_west) * &deta_dphi_d_volume;

        let cxex = xi_cent.mapv(|_| 1.0);
        let cpep = xi_cent.mapv(|_| 1.0);
        let cxpx = xi_cent.mapv(|_| 1.0);
        let cepe = xi_cent.mapv(|_| 1.0);

        let minus_cexe = cexe.mapv(|c| -1.0 * c);
        let minus_cpxp = cpxp.mapv(|c| -1.0 * c);
        let minus_cxex = cxex.mapv(|c| -1.0 * c);
        let minus_cpep = cpep.mapv(|c| -1.0 * c);
        let minus_cxpx = cxpx.mapv(|c| -1.0 * c);
        let minus_cepe = cepe.mapv(|c| -1.0 * c);

        let mesh = Mesh {
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
        };

        if cfg!(feature = "validation") {
            mesh.validate().context("Validating Mesh")?;
        }

        Ok(mesh)
    }
}

impl<const S: usize> Validation for Mesh<S> {
    fn validate(&self) -> Result<()> {
        // checks on doubles
        check_finite_multiple_doubles![self.dxi, self.deta, self.dphi];

        // checks on arrays
        check_nonempty_finite_multiple_arrayd![
            self.xi_cent,
            self.xi_west,
            self.xi_east,
            self.xi_cent_inv,
            self.xi_west_inv,
            self.xi_east_inv,
            self.h_xi_cent,
            self.h_xi_west,
            self.h_xi_east,
            self.h_eta_cent,
            self.h_eta_west,
            self.h_eta_east,
            self.h_phi_cent,
            self.h_phi_west,
            self.h_phi_east,
            self.sqrt_g,
            self.line_xi,
            self.line_xi_inv,
            self.d_area_xi_deta_dphi_west,
            self.d_area_xi_deta_dphi_east,
            self.area_cell,
            self.area_west,
            self.area_east,
            self.volume,
            self.deta_dphi_d_volume,
            self.cell_width,
            self.cell_width_inv,
            self.cexe,
            self.cpxp,
            self.cxex,
            self.cpep,
            self.cxpx,
            self.cepe,
            self.minus_cexe,
            self.minus_cpxp,
            self.minus_cxex,
            self.minus_cpep,
            self.minus_cxpx,
            self.minus_cepe
        ];

        Ok(())
    }
}

impl<const S: usize> Collectable for Mesh<S> {
    fn collect_data(&self, data: &mut Data, mesh_offset: usize) -> Result<()> {
        match (data.association, data.name) {
            (StructAssociation::Mesh, DataName::XiCent) => self.write_vector(&self.xi_cent.view(), data, mesh_offset),
            (StructAssociation::Mesh, DataName::XiWest) => self.write_vector(&self.xi_west.view(), data, mesh_offset),
            (StructAssociation::Mesh, DataName::XiEast) => self.write_vector(&self.xi_east.view(), data, mesh_offset),
            (StructAssociation::Mesh, x) => bail!("Tried associating {:?} with Mesh!", x),
            (StructAssociation::Physics, x) | (StructAssociation::TimeStep, x) => {
                bail!("name.association() for {:?} returned {:?}", x, data.association)
            },
        }?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    const S: usize = 1004;
    // TODO: more test submodules with different sized meshes?
    prop_compose! {
        fn arb_mesh(mode: MeshMode)
                    (xi_in in 0.1f64..100_000.0,
                    xi_out in 0.1f64..100_000.0,
                    ratio_disk in 0.1f64..=1.0) -> Result<Mesh<S>> {
            Mesh::new(&MeshConfig { mode, xi_in, xi_out, ratio_disk })
        }
    }

    mod cartesian_meshes {
        use super::*;
        use approx::assert_relative_eq;
        proptest! {
            #[test]
            fn cartesian_meshes_constructed_correctly(mesh in arb_mesh(MeshMode::Cartesian)) {
                prop_assume!(mesh.is_ok());
                let mesh = mesh.unwrap();
                prop_assume!(mesh.xi_in < mesh.xi_out);
                assert_relative_eq!(mesh.xi_in, mesh.xi_west[mesh.ixi_in], max_relative = 1.0e-12);
                assert_relative_eq!(mesh.xi_out, mesh.xi_east[mesh.ixi_out], max_relative = 1.0e-12);
            }
        }
    }
}
