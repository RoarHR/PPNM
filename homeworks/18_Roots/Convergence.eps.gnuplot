set terminal svg enhanced background "white"
set output 'Convergence.eps.svg'

set title "Convergence with respect to eps"
set xlabel "eps"
set ylabel "E"
set logscale x

f0(r) = -0.5

plot "Out.Convergence.eps.txt" using 4:5 with lp title "E vs eps", \
     f0(x) with lines title "Exact solution"

