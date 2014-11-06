#!/bin/bash
#
# Syntax:
#
# grip_install.sh path\file.md 
#
# Convert `path\file.md` markdown file to HTML file 'file.html' in build_tree.
#
# Use (and possibly install grip) for conversion.
# Install GRIP python package for conversion of GitHub Markdown files to HTML.


set -x

MD_FILE="$1"

BUILD_TREE="`pwd`/${0%/*}/../build_tree"
PREFIX="${BUILD_TREE}/local"

SAVE_DIR="`pwd`"
cd "${PREFIX}"

# if no system virtualenv, install our own
if ! virtualenv --help >/dev/null
then
  VERSION=virtualenv-1.9
  wget https://pypi.python.org/packages/source/v/virtualenv/${VERSION}.tar.gz
  tar xvfz ${VERSION}.tar.gz
  VIRTUALENV="python ${VERSION}/virtualenv.py"
else
  VIRTUALENV=virtualenv
fi

if ! "${PREFIX}/bin/grip" --help >/dev/null
then 
  # make virtual python environment
  python "${VIRTUALENV}" "${PREFIX}" 
  # second call since for first call we get some error under cygwin64
  # note also that there was a bug in Cygwin64, that make some install tools
  # cripple. Namely, pip does nothing. Mysterious workaround is to install libuuid-devel package.
  ${VIRTUALENV} "${PREFIX}" 
  source "${PREFIX}/bin/activate"
  "${PREFIX}/bin/pip"  install grip
fi

cd "${SAVE_DIR}"
"${PREFIX}/bin/grip" --export --gfm --context=https://github.com/flow123d/flow123d.git "${MD_FILE}"
TARGET="${MD_FILE%.md}.html"
mv "${TARGET}" "${BUILD_TREE}"