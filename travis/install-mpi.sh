#!/bin/sh
set -e
case $1 in
  mpich) set -x;
    sudo apt-get install -q mpich libmpich-dev;;
  openmpi) set -x;
    sudo apt-get install openmpi-bin openmpi-dev;;
  *)
    echo "Unknown MPI implementation:" $1; exit 1;;
esac
