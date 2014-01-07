#!/usr/bin/env bash

# Make sure the package information is up-to-date
apt-get update

## Qt5 visualization
#apt-get install -y qt5-default
#apt-get install -y imagemagick

# Compilers
apt-get install -y g++-4.7
apt-get install -y gfortran-4.7
apt-get install -y clang-3.4

# Message Passing Interface
apt-get install -y libmpich2-dev

# BLAS and LAPACK
apt-get install -y libopenblas-dev
apt-get install -y liblapack-dev

# Configuration
apt-get install -y cmake

# Source control
apt-get install -y git

# Check out Elemental
git clone http://github.com/elemental/Elemental.git
chown -R vagrant Elemental
