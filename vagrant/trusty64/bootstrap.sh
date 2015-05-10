#!/usr/bin/env bash

# Make sure the package information is up-to-date
apt-get update

## Qt5 visualization
apt-get install -y qt5-default
#apt-get install -y imagemagick

# Python visualization
apt-get install -y python-matplotlib python-networkx

# Compilers
apt-get install -y gfortran-4.8 clang-3.5

# Message Passing Interface
apt-get install -y libmpich2-dev

# Configuration
apt-get install -y cmake

# Source control
apt-get install -y git

# Check out Elemental
git clone https://github.com/elemental/Elemental
chown -R vagrant Elemental
