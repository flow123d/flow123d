#!/bin/bash

# make input files from *.geo

LARGE="large_mesh.geo"
SMALL="small_mesh.geo"
NGH=../../../bin/ngh/bin/ngh
BCD=../../../bin/bcd/bin/bcd

gmsh -2 -algo front2d -3 -algo front3d $LARGE
gmsh -2 -algo front2d -3 -algo front3d $SMALL

cat ngh.ini | sed "s/{MESH}/${LARGE%.geo}/" >ngh_large.ini
$NGH ngh_large.ini
cat ngh.ini | sed "s/{MESH}/${SMALL%.geo}/" >ngh_small.ini
$NGH ngh_small.ini

$BCD bcd.ini
