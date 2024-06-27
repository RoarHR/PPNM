set terminal svg enhanced background "white"
set output 'Convergence.rmin.svg'

set title "Convergence with respect to rmin"
set xlabel "rmin"
set ylabel "E"
set logscale x

f0(r) = -0.5

plot "Out.Convergence.rmin.txt" using 1:5 with lp title "E vs rmin", \
     f0(x) with lines title "Exact solution"

