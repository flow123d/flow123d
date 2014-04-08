#!/bin/bash

# syntax: 
# convert_data.sh file
# - perform a substitution on given file
# 
#
# convert_data.sh -a MASK 
# - perform a substitution on all files with relative path matching given MASK
#

set -x
function convert_file {

  FILE=$1

  # skip directories
  if [ -d $FILE ] ; then return; fi

  echo "Converting $FILE"
  TMP=tmp
  cp $FILE $TMP
  
  ####################################################################
  # substitute possible boundary between bc and bulk data lists
  #perl -i.orig -0pe  's/}[ \n]*\][ \n]*,[ \n]*\n( *)bc_data[ \n]*=[ \n]*\[[ \n]*{/},\n\1  {/g' $TMP
  #perl -i.orig -0pe  's/}[ \n]*\][ \n]*,[ \n]*\n( *)bulk_data[ \n]*=[ \n]*\[[ \n]*{/},\n\1  {/g' $TMP

  # substitute remaining data key
  #perl -i.orig -0pe  's/bc_data[ \n]*=/data=/g' $TMP
  #perl -i.orig -0pe  's/bulk_data[ \n]*=/data=/g' $TMP

  #####################################################################
  # dual porosity
  #perl -i.orig -0pe  's/alpha[ \n]*=/diffusion_rate_immobile=/g' $TMP
  #perl -i.orig -0pe  's/immob_porosity[ \n]*=/porosity_immobile=/g' $TMP
  
  #sorption
  #perl -i.orig -0pe  's/solvent_dens[ \n]*=/solvent_density=/g' $TMP
  #perl -i.orig -0pe  's/molar_masses[ \n]*=/molar_mass=/g' $TMP
  #perl -i.orig -0pe  's/mult_coefs[ \n]*=/isotherm_mult=/g' $TMP
  #perl -i.orig -0pe  's/second_params[ \n]*=/isotherm_other=/g' $TMP
  #perl -i.orig -0pe  's/sorption_types[ \n]*=/adsorption_type=/g' $TMP
  
  #porosity
  #perl -i.orig -0pe  's/por_m[ \n]*=/porosity=/g' $TMP
  
  ########################################################################
  # substitute inf time -> 0 time for steady ref data
  perl -i.orig -0pe 's/inf /0 /g' $TMP
  
  mv tmp $FILE
  rm tmp.orig
}  

# set path to current script
THIS=`pwd`/${0} 
if [ "$1" == "-a" ] 
then
  find . -path "$2" -execdir $THIS '{}' \;
else
  convert_file $1
fi