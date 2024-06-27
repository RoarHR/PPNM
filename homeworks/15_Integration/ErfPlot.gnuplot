set terminal svg enhanced background "white"
set output 'ErfPlot.svg'

set xlabel 'x'
set ylabel 'erf(x)'

# Title
set title 'Error Function: Numerical vs Tabulated'

# Plot numerical and tabulated values
plot 'Out.Erf.txt' using 1:2 with lines title 'Numerical', \
     'err_tabulated.txt' using 1:2 with points title 'Tabulated'
