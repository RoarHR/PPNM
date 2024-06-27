set terminal svg enhanced background "white"
set output 'Convergence.rmax.svg'

set title "Convergence with respect to rmax"
set xlabel "rmax"
set ylabel "E"
set yrange [-0.501:-0.482]

f0(r) = -0.5

plot "Out.Convergence.rmax.txt" using 2:5 with lp title "E vs rmax", \
     f0(x) with lines title "Exact solution"

