set terminal svg enhanced background "white"
set output 'Linterp.svg'

set xlabel 'x'
set ylabel 'y'

# Title
set title 'Linear interpolation of Cosine function'

plot 'Out.Linterp.interp.txt' using 1:2 w histeps title 'Interpolation', \
     'Out.Linterp.xy.txt' using 1:2 w p title 'Points'
