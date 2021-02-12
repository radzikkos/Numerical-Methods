#!/usr/bin/gnuplot
set term png

set output "integral.png"
set xlabel "it"
set ylabel "S(it)"
set title "Całka funkcjonalna S(it)"

plot "integral.dat" u 2:3 w l t ""

set output "map_k_16.png"
set xlabel "x"
set ylabel "y"
reset
set pm3d map
set palette defined (-1 "blue", 0 "yellow", 1 "red")
set size ratio -1

splot [0:26][0:26] "maps.dat" i 0 u 1:2:3 notitle

reset
set output "map_k_8.png"
set xlabel "x"
set ylabel "y"
set pm3d map
set palette defined (0 "blue", 5 "yellow", 10 "red")
set size ratio -1

splot [0:26][0:26] "maps.dat" i 1 u 1:2:3 notitle


reset
set output "map_k_4.png"
set xlabel "x"
set ylabel "y"
set pm3d map
set palette defined (0 "blue", 5 "yellow", 10 "red")
set size ratio -1

splot [0:26][0:26] "maps.dat" i 2 u 1:2:3 notitle


reset
set output "map_k_2.png"
set xlabel "x"
set ylabel "y"
set pm3d map
set palette defined (0 "blue", 5 "yellow", 10 "red")
set size ratio -1

splot [0:26][0:26] "maps.dat" i 3 u 1:2:3 notitle


reset
set output "map_k_1.png"
set xlabel "x"
set ylabel "y"
set pm3d map
set palette defined (0 "blue", 5 "yellow", 10 "red")
set size ratio -1

splot [0:26][0:26] "maps.dat" i 4 u 1:2:3 notitle