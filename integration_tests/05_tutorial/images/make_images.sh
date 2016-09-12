#!/bin/bash

# Tutorial 01

# make image of mesh, bc and flow field
gmsh make_01_flow.geo
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
plot '<grep rock ../ref_out/02_column_transport/mass_balance.txt' u (\$1/86400/365):7 w l axes x1y1 t 'rock',\\
     '<grep .surface ../ref_out/02_column_transport/mass_balance.txt' u (\$1/86400/365):4 w l axes x1y2 t '.surface',\\
     '<grep .tunnel ../ref_out/02_column_transport/mass_balance.txt' u (\$1/86400/365):4 w l axes x1y2 t '.tunnel'" | gnuplot

# make image with concentrations
gmsh make_02_transport.geo
pdfjam 02_transport_1.pdf 02_transport_2.pdf 02_transport_3.pdf 02_transport_4.pdf --nup 4x1 --papersize '{20cm,10cm}' --outfile 02_transport.pdf
rm 02_transport_[1234].pdf


# Tutorial 03

#make plot of mass flux
echo "
set terminal pdf color enhanced lw 2
set output '03_mass_flux.pdf'
set xdata time
set format x '%Y'
set xtics 31557600
unset mxtics
set ytics nomirror
set y2tics
set yrange [-11:-9]
set y2range [-18:-4]
set key top center out horizontal
set xlabel 'time [years]'
set ylabel '{/Symbol d}18O in tunnel [permil V-SMOW]'
set y2label '{/Symbol d}18O in precipitation [permil V-SMOW]'
set style circle radius graph 0.004
plot '03_conc_tunnel.txt' u (timecolumn(1,\"%d.%m.%Y\")):2 w l lw 0.1 lc \"00000000\" axes x1y1 not,\\
     '' u (timecolumn(1,\"%d.%m.%Y\")):2 w circles fs solid axes x1y1 t 'measured',\\
     '03_conc_surface.txt' u (timecolumn(1,\"%d.%m.%Y\")):2 w l lc rgb \"#d0ff0000\" axes x1y2 t 'precipitation',\\
     '<grep .tunnel ../ref_out/03_tunnel/mass_balance.txt' u ((timecolumn(1,\"%s\")+36*365.25)*86400):(-\$4) w l axes x1y1 t 'computed'
" | gnuplot

# make image of flux and pressure
gmsh make_03.geo



# Tutorial 05

#image with bc and mesh
gmsh make_05.geo
pdfjam 05_bc.pdf 05_mesh.pdf --nup 2x1 --papersize '{20cm,10cm}' --outfile 05_bcmesh.pdf
rm -rf 05_mesh.pdf

#make plot of production power
echo "
set terminal pdf color enhanced lw 2
set output '05_power.pdf'
set xlabel 'time [years]'
plot [1:] \"< sed -n '/\\\\.well[12]_surface\\\"/p' ../ref_out/05_heat/energy_balance.txt | awk 'BEGIN {t=0;v=0} {if (t==\$1) v+=\$4; else {print t, v; t=\$1;v=\$4}} END {print t,v}'\" u (\$1/365.25/86400):(-\$2/1e6) w lp t 'power [MW]',\\
          25 dt 2 not
" | gnuplot
