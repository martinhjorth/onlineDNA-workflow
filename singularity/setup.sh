#!/usr/bin/env bash
set -ex
#allow passing on builddir from elsewhere, otherwise use default
builddir=${builddir:-$(realpath "./build")}
mkdir -p $builddir
pushd $builddir

#get miniconda
miniconda="miniconda.sh"
if [ ! -f ${miniconda} ]
then
  wget -q https://repo.anaconda.com/miniconda/Miniconda3-py37_4.9.2-Linux-x86_64.sh -O ${miniconda}
fi

popd
