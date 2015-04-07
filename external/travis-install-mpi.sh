#!/bin/sh
set -e
case $1 in
  mpich2) set -x;
    sudo apt-get install -q mpich2 libmpich2-3 libmpich2-dev;;
  mpich3) set -x;
    sudo apt-get install -q gfortran libcr0 default-jdk;
    wget -q http://www.cebacad.net/files/mpich/ubuntu/mpich-3.1/mpich_3.1-1ubuntu_amd64.deb;
    sudo dpkg -i ./mpich_3.1-1ubuntu_amd64.deb;
    rm -f ./mpich_3.1-1ubuntu_amd64.deb;;
  openmpi) set -x;
    sudo apt-get install openmpi-bin openmpi-dev;;
  *)
    echo "Unknown MPI implementation:" $1; exit 1;;
esac
