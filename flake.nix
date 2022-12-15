{
  description = "A corrosive 1D hydrodynamics solver";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    rust-overlay.url = "github:oxalica/rust-overlay";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, rust-overlay, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        overlays = [ (import rust-overlay) ];
        pkgs = import nixpkgs {
          inherit system overlays;
        };
      in
      with pkgs;
      {
        devShells.default = mkShell {
          buildInputs = [
            openssl
            pkg-config
            # (rust-bin.selectLatestNightlyWith (toolchain: toolchain.default.override {
            (rust-bin.stable.latest.default.override {
              extensions = [ "rust-src" "rust-analyzer" "rustfmt" ];
            })

            # external deps
            openblas

            # dev tools
            bacon
            nixpkgs-fmt
            rust-analyzer
            cargo-nextest # https://nexte.st/
            cargo-criterion # https://github.com/bheisler/cargo-criterion
          ];

          # shellHook = ''
          #   alias ls=exa
          #   alias find=fd
          # '';
        };
      }
    );
}
