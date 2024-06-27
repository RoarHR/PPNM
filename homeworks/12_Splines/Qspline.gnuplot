set terminal svg enhanced background "white"
set output 'Qspline.svg'

set xlabel 'x'
set ylabel 'y'

# Title
set title 'Linear interpolation of sine function'

f(x) = cos(x)

plot 'Out.Qspline.spline.txt' using 1:2 w histeps title 'Spline', \
     'Out.Qspline.xy.txt' using 1:2 w p title 'Points', \
     'Out.Qspline.spline.txt' using 1:3 w lp title 'Derivative', \
     f(x) title 'f(x) = cos(x)' lw 2

