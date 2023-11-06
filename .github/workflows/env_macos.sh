#!/bin/bash


git clone https://github.com/spack/spack.git

#replace the default configuration file of spack by a simpler one without hash, compiler versions, tags and so on 
cp /Users/runner/work/gmds/gmds/.github/workflows/misc/config-0.20.3.yaml /Users/runner/work/gmds/gmds/spack/etc/spack/defaults/config.yaml
. ./spack/share/spack/setup-env.sh

git clone --branch gmds_temp --depth=1 https://github.com/LIHPC-Computational-Geometry/spack_recipes.git

spack repo add ./spack_recipes/meshing_repo
spack repo add ./spack_recipes/supersede_repo

spack external find cmake
spack install --only dependencies gmds+blocking~cgns
