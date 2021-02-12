#!/usr/bin/gnuplot
set term png

##########################################################

#Mapa zrelaksowanego potencjalu dla omegi równej 1.0
set output "map_1.0.png"
set xlabel "x"
set ylabel "y"
set title "Mapa zrelaksowanego potencjału V(x,y) dla omegi równej 1.0"

set pm3d map
set palette defined (0 "green", 5 "white", 10 "red")
set size ratio -1

splot [0:15][0:10] "map_1.0_relglob.dat" i 0 u 1:2:3

#Mapa zrelaksowanego potencjalu dla omegi równej 0.6
set output "map_0.6.png"
set xlabel "x"
set ylabel "y"
set title "Mapa zrelaksowanego potencjału V(x,y) dla omegi równej 0.6"

set pm3d map
set palette defined (0 "green", 5 "white", 10 "red")
set size ratio -1

splot [0:15][0:10] "map_0.6_relglob.dat" i 0 u 1:2:3

##########################################################