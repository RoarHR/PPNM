set terminal svg enhanced background "white"
set output 'Convergence.acc.svg'

set title "Convergence with respect to acc"
set xlabel "acc"
set ylabel "E"
set logscale x
set yrange [-0.50001:-0.49984]

f0(r) = -0.5

plot "Out.Convergence.acc.txt" using 3:5 with lp title "E vs acc", \
     f0(x) with lines title "Exact solution"

