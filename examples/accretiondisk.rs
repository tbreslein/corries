use color_eyre::Result;
use corries::config;
use corries::run_sim;

fn main() -> Result<()> {
    run_sim(config::CorriesConfig {
        name: "accretiondisk".to_string(),
        meshconf: config::meshconfig::MeshConfig {
            mode: config::meshconfig::MeshMode::Cartesian,
            n_comp: 100,
            n_gc: 2,
            xi_in: 0.1,
            xi_out: 100.0,
            ratio: 1.0,
        },
    })
}
