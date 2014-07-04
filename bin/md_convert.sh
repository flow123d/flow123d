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


MD_FILE="$1"

BUILD_TREE="`pwd`/${0%/*}/../build_tree"
PREFIX="${BUILD_TREE}/local"

OUR_PYTHONPATH="${PREFIX}/lib/python2.7/site-packages"
OUR_PATH="${PREFIX}/bin"

export PYTHONPATH="${PYTHONPATH}:${OUR_PYTHONPATH}"
export PATH="${PATH}:${OUR_PATH}"

if [ ! -d "${OUR_PATH}" ]
then
  mkdir -p "${OUR_PATH}"
fi
if [ ! -d "${OUR_PYTHONPATH}" ]
then
  mkdir -p "${OUR_PYTHONPATH}"
fi

SAVE_DIR="`pwd`"
cd "${PREFIX}"

# check and install easy_install
if ! easy_install --help >/dev/null
then
  # install
  wget http://peak.telecommunity.com/dist/ez_setup.py
  python ez_setup.py --prefix "${PREFIX}"
  rm -f ez_setup.py
  
  # check
  if ! easy_install --help >/dev/null
  then
    echo "can not install 'easy_install' ... GIVE UP"
    exit
  else
    echo "easy_install ... OK"
  fi  
else
  echo "easy_install ... OK"
fi  



# check and install GRIP
if ! grip --help >/dev/null
then
  # install
  wget http://bacula.nti.tul.cz/~jan.brezina/flow123d_libraries/grip-master.tar.gz
  tar -xzf grip-master.tar.gz
  cd grip-master
  python setup.py install --prefix "${PREFIX}"
  cd ..
  rm -rf grip-master
  rm grip-master.tar.gz
  
  # check
  if ! grip --help >/dev/null
  then
    echo "can not install 'grip' ... GIVE UP"
    exit
  else
    echo "grip... OK"
  fi  
else
  echo "grip ... OK"
fi  

cd "${SAVE_DIR}"

${PREFIX}/bin/grip --export --gfm --context=https://github.com/flow123d/flow123d.git "${MD_FILE}"
TARGET="${MD_FILE%.md}.html"
mv "${TARGET}" "${BUILD_TREE}"