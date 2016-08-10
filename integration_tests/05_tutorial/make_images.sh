#!/bin/bash

# Tutorial 01

# make image of mesh, bc and flow field
gmsh make_flow.geo
pdfjam 01_mesh.pdf 01_bc.pdf 01_flow.pdf --nup 3x1 --papersize '{20cm,10cm}' --delta '2cm 0cm' --outfile 01_mesh_bc_flux.pdf
rm 01_mesh.pdf 01_flow.pdf


# Tutorial 02

# make plot of mass and mass flux
echo "
set terminal pdf color enhanced lw 2
set output '02_mass_plot.pdf'
set xlabel 'time [years]'
set ylabel 'mass in the rock [kg]'
set y2label 'mass flux through .surface/.tunnel [kg/s]'
set key top center out horizontal
set ytics nomirror
set y2tics
plot '<grep rock ref_out/02_column_transport/mass_balance.txt' u (\$1/86400/365):7 w l axes x1y1 t 'rock',\\
     '<grep .surface ref_out/02_column_transport/mass_balance.txt' u (\$1/86400/365):4 w l axes x1y2 t '.surface',\\
     '<grep .tunnel ref_out/02_column_transport/mass_balance.txt' u (\$1/86400/365):4 w l axes x1y2 t '.tunnel'" | gnuplot

# make image with concentrations
gmsh make_transport.geo
pdfjam 02_transport_1.pdf 02_transport_2.pdf 02_transport_3.pdf 02_transport_4.pdf --nup 4x1 --papersize '{20cm,10cm}' --outfile 02_transport.pdf
rm 02_transport_[1234].pdf