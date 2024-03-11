#!/bin/bash

# Install binaries of Rust language
mamba create -c conda-forge -n hyperex rust -y

# Activate the environment created
mamba activate hyperex

# Install the desired package
cargo install hyperex

# Deactivate the environemnt
mamba deactivate hyperex

# https://github.com/Ebedthan/hyperex