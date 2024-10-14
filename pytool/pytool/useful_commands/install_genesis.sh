#!/bin/bash
set -xe

# check if the sfotware directory exists
if [ ! -d software ]; then
    echo "Please create the directory software"
    exit 1
fi

genesis_name=genesis-2.1.0

mkdir -p $genesis_name

# configura parameters
# --enable-single
#   turn on single precisionc alculation. in this case, only spdyn is installed
# --enable-mixed
#   turn on mixed precision calculation. in this case, only spdyn is installed
# --enable-double (default)
#   turn on double precision calculation. in this case, all binaries are installed 
# --enable-gpu
#   turn on GPGPU calculation. In this case, only SPDYN is intalled
# --with--cuda=PATH
#   define path to the CUDA libraries manually



curl -fsSL "https://www.r-ccs.riken.jp/labs/cbrt/?sdm_process_download=1&download_id=25931" -o ${genesis_name}.tar.bz2
tar xvfj ${genesis_name}.tar.bz2 -C ./${genesis_name} --strip-components 1


cd ${genesis_name}
#autoconf
./configure "$@"
make -j$(fgrep 'cpu cores' /proc/cpuinfo | sort -u | sed 's/.*: //')
make install
