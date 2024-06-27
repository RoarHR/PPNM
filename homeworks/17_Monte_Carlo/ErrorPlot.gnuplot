set terminal svg enhanced background "white"
set output 'ErrorPlot.svg'

set xlabel 'N'
set ylabel 'Error'
set logscale x
set logscale y
set title 'Error vs. N for PlainMC and QuasiMC'

# Fit the plainMC error to a 1/sqrt(N) curve
f(N) = a / sqrt(N)
fit f(x) 'Out.ErrorPlot.txt' using 1:2 via a

# Plot the errors
plot 'Out.ErrorPlot.txt' using 1:2 title 'PlainMC Estimated Error' with linespoints, \
     'Out.ErrorPlot.txt' using 1:3 title 'PlainMC |Analytical - Numerical|' with linespoints, \
     'Out.ErrorPlot.txt' using 1:4 title 'QuasiMC Estimated Error' with linespoints, \
     'Out.ErrorPlot.txt' using 1:5 title 'QuasiMC |Analytical - Numerical|' with linespoints, \
     f(x) title 'PlainMC 1/sqrt(N) Fit' with lines
