#!/bin/bash

sudo apt-get install software-properties-common
sudo add-apt-repository ppa:george-edison55/cmake-3.x
sudo apt-get update

sudo apt-get install cmake
sudo apt-get install libgsl-dev
sudo apt-get install r-base-dev

SCRIPTDIR="$( cd $( dirname ${BASH_SOURCE[0]} ) >/dev/null 2>&1 && pwd )"

${SCRIPTDIR}/packages.R
${SCRIPTDIR}/packages.R
