#!/bin/sh
set -e
case $1 in
  mpich2) set -x;
    sudo apt-get install -q mpich2 libmpich2-3 libmpich2-dev;;
  openmpi) set -x;
    sudo apt-get install openmpi-bin openmpi-dev;;
  *)
    echo "Unknown MPI implementation:" $1; exit 1;;
esac
