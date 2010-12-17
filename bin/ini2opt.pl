#!/usr/bin/perl
# 
# Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
#
# Please make a following refer to Flow123d on your project site if you use the program for any purpose,
# especially for academic research:
# Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
#
# This program is free software; you can redistribute it and/or modify it under the terms
# of the GNU General Public License version 3 as published by the Free Software Foundation.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not,
# write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
#
# $Id$
# $Revision$
# $LastChangedBy$
# $LastChangedDate$
#

#
# A script for conversion from ".ini" files fo PETSc option database files
#
#
# ini2opt file1.ini [file2]
# 
# default file2=file1.opt
#

$InFName=shift;
if ($InFName eq '') {die "usage: ini2opt file1.ini [file2]\n";}
$OutFName=shift;
if ($OutFName eq '') {
  $OutFName=$InFName;
  $OutFName=~s/\.ini$/.opt/;
}

open(IN,"<$InFName") || die "Can't find input file $InFName.\n";
open(OUT,">$OutFName") || die "Can't open output file $OutFName.\n";

while ($line=<IN>) {
  chop $line;
  $line=~s/^\r//;
  if ($line=~ /\s*\[([a-zA-Z_ ]*)\]\s*(\;.*)?/) {
    # new section
    $section=$1;
    $section=~tr/ /_/;
    print OUT "#$line\n";
  }
  elsif ($line=~ /\s*(\w*)\s*=\s*([^;\r]*)(\;.*)?(\r?)/) {
    # value line
    if ($section eq '') {print "warning: parameter out of section\n";}
    else {
      $lineend=$3;
      $lineend=~s/\;/#/;
      print OUT "-".$section."_".$1." \"$2\"".$endline.$4."\n";}
  }
  else {
    print OUT "#$line\n";
  }
}

close IN;
close OUT;


