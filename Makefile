jobs := $(or $(jobs),$(j))
buildargs := $(or $(buildargs),$(b))
configargs := $(or $(configargs),$(c))
testargs := $(or $(testargs),$(t))

.PHONY: format
format:
	rustfmt ./**/*.rs

.PHONY: clippy
clippy:
	cargo clippy -- -A clippy::needless_return -A clippy::op_ref

.PHONY: 
check:
	cargo check

.PHONY: doc
doc:
	cargo doc --document-private-items

.PHONY: build
build:
	cargo build -r

.PHONY: unittest
unittest:
	cargo test --lib

.PHONY: reltest
test:
	cargo test -r

.PHONY: bench
bench:
	cargo test --bench -r

.PHONY: clean
clean:
	cargo clean

.PHONY: runlints
runlints: format check clippy unittest

.PHONY: fulltests
fulltests: runlints reltest

.PHONY: fullbuild
fullbuild: runlints build doc
