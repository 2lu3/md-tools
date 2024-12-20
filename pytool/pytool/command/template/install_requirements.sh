#!/bin/bash
set -xe

# gcc & gfortran
sudo apt-get install -y gfortran gcc g++ build-essential

# openmpi
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.4.tar.gz -O /tmp/openmpi.tar.gz

mkdir /tmp/openmpi && tar -xvf /tmp/openmpi.tar.gz -C /tmp/openmpi --strip-components 1

mkdir -p ~/.local/
pushd /tmp/openmpi/
./configure --prefix=${HOME}/.local CC=gcc CXX=g++ F77=gfortran FC=gfortran
make -j$(fgrep 'cpu cores' /proc/cpuinfo | sort -u | sed 's/.*: //')
make install
popd

# lapack
sudo apt-get install -y liblapack-dev libblas-dev
