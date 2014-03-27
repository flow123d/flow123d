#!/bin/bash

# syntax: 
# convert_data.sh file
#
# convert_data.sh -a 
#
# convert a given *.con file from  format with "bc_data", "bulk_data" pair of lists to the format
# using just one "data" list per equation.
#  
# second variant convert all con files in all directories under current one
#

set -x
function convert_file {

  FILE=$1

  # skip directories
  if [ -d $FILE ] ; then return; fi

  echo "Converting $FILE"
  TMP=tmp
  cp $FILE $TMP
  # substitute possible boundary between bc and bulk data lists
  #perl -i.orig -0pe  's/}[ \n]*\][ \n]*,[ \n]*\n( *)bc_data[ \n]*=[ \n]*\[[ \n]*{/},\n\1  {/g' $TMP
  #perl -i.orig -0pe  's/}[ \n]*\][ \n]*,[ \n]*\n( *)bulk_data[ \n]*=[ \n]*\[[ \n]*{/},\n\1  {/g' $TMP

  # substitute remaining data key
  #perl -i.orig -0pe  's/bc_data[ \n]*=/data=/g' $TMP
  #perl -i.orig -0pe  's/bulk_data[ \n]*=/data=/g' $TMP

  # dual porosity
  #perl -i.orig -0pe  's/alpha[ \n]*=/diffusion_rate_immobile=/g' $TMP
  #perl -i.orig -0pe  's/immob_porosity[ \n]*=/porosity_immobile=/g' $TMP
  
  #sorption
  perl -i.orig -0pe  's/solvent_dens[ \n]*=/solvent_density=/g' $TMP
  perl -i.orig -0pe  's/molar_masses[ \n]*=/molar_mass=/g' $TMP
  perl -i.orig -0pe  's/mult_coefs[ \n]*=/isotherm_mult=/g' $TMP
  perl -i.orig -0pe  's/second_params[ \n]*=/isotherm_other=/g' $TMP
  perl -i.orig -0pe  's/sorption_types[ \n]*=/adsorption_type=/g' $TMP
  
  mv tmp $FILE
}  

# set path to current script
THIS=`pwd`/${0} 
if [ "$1" == "-a" ] 
then
  find . -name "*.con" -execdir $THIS '{}' \;
else
  convert_file $1
fi