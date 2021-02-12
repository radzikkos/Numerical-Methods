#!/usr/bin/gnuplot
set term png

set output "map1a.png"
set xlabel "x"
set ylabel "y"
set title "nx = ny = 50"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:5][0:5] "map1a.dat" i 0 u 2:1:3

reset

set output "map1b.png"
set xlabel "x"
set ylabel "y"
set title "nx = ny = 100"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:10][0:10] "map1b.dat" i 0 u 2:1:3

reset

set output "map1c.png"
set xlabel "x"
set ylabel "y"
set title "nx = ny = 200"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:20][0:20] "map1c.dat" i 0 u 2:1:3

reset

set output "map2a.png"
set xlabel "x"
set ylabel "y"
set title " eps1 = eps2 = 1"
set pm3d map
set palette defined (-1 "blue", 0 "white", 1 "red")
set size ratio -1

splot [0:10][0:10][-0.8:0.8] "map2a.dat" i 0 u 2:1:3


reset

set output "map2b.png"
set xlabel "x"
set ylabel "y"
set title " eps1 = 1 eps2 = 2"
set pm3d map
set palette defined (-1 "blue", 0 "white", 1 "red")
set size ratio -1

splot [0:10][0:10][-0.8:0.8] "map2b.dat" i 0 u 2:1:3


reset

set output "map2c.png"
set xlabel "x"
set ylabel "y"
set title "2c) eps1 = 1 eps2 = 10"
set pm3d map
set palette defined (-1 "blue", 0 "white", 1 "red")
set size ratio -1

splot [0:10][0:10][-0.8:0.8] "map2c.dat" i 0 u 2:1:3