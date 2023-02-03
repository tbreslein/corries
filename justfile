clean:
  cargo clean

format:
  rustfmt ./**/*.rs

clippy:
  cargo clippy

doc:
  cargo doc --document-private-items

build: format clippy doc
  cargo build -r

test: build
  cargo nextest run -E 'binary(corries)'
  cargo test --doc
  cargo nextest run -r -E 'kind(test)'

update:
  nix flake update
  cargo update

bench:
  cargo bench
