#!/usr/bin/gnuplot
set term png

##########################################################

#Calka funkcjonalna dla relaksacji globalnej
set output "integral_global.png"
set logscale x
set xlabel "it"
set ylabel "S(it)"
set title "Całka funkcjonalna S(it) dla relaksacji globalnej"
plot "integral_1.0_relglob.dat" u 1:2 w l t "omega = 1.0" , "integral_0.6_relglob.dat" u 1:2 w l lw 2 t "omega = 0.6"

#Calka funkcjonalna dla relaksacji lokalnej
set output "integral_local.png"
set logscale x
set logscale x
set xlabel "it"
set ylabel "S(it)"
set title "Całka funkcjonalna S(it) dla relaksacji lokalnej"
plot "integral_1.0_relloc.dat" u 1:2 w l t "omega = 1.0" , "integral_1.4_relloc.dat" u 1:2 w l lw 2 t "omega = 1.4", "integral_1.8_relloc.dat" u 1:2 w l lw 2 t "omega = 1.8", "integral_1.9_relloc.dat" u 1:2 w l lw 2 t "omega = 1.9"


##########################################################