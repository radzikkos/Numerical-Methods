#!/usr/bin/gnuplot
set term png

##########################################################

#Blad rozwiazania
set output "error_1.0.png"
set xlabel "x"
set ylabel "y"
set title "Błąd rozwiazania dla relaksacji globalnej i omegi równej 1.0"

set pm3d map
set palette rgbformulae 33,13,10
set size ratio -1

splot [0:15][0:10] "error_1.0_relglob.dat" i 0 u 1:2:3


set output "error_0.6.png"
set xlabel "x"
set ylabel "y"
set title "Błąd rozwiazania dla relaksacji globalnej i omegi równej 0.6"

set pm3d map
set palette rgbformulae 33,13,10
set size ratio -1

splot [0:15][0:10] "error_0.6_relglob.dat" i 0 u 1:2:3

##########################################################