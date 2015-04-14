#!/bin/sh
set -e
case $1 in
  mpich2) set -x;
    sudo apt-get install -q mpich2 libmpich2-3 libmpich2-dev;;
  openmpi) set -x;
    sudo apt-get install alien;
    wget -q http://www.open-mpi.org/software/ompi/v1.8/downloads/openmpi-1.8.4-1.src.rpm;
    sudo alien -d ./openmpi-1.8.4-1.src.rpm;
    sudo dpkg -i ./openmpi-1.8.4-1.deb;
    rm -f ./openmpi-1.8.4-1.src.rpm;
    rm -f ./openmpi-1.8.4-1.deb;;
  *)
    echo "Unknown MPI implementation:" $1; exit 1;;
esac
