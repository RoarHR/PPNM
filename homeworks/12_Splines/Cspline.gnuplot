set terminal svg enhanced background "white"
set output 'Cspline.svg'

set xlabel 'x'
set ylabel 'y'

# Title
set title 'Cubic spline of ugly dataset'

plot 'Out.Cspline.spline.txt' using 1:2 w lines title 'My Cubic Spline', \
     'Cspline.xy.txt' using 1:2 smooth csplines w p pt 7 ps 0.5 title 'Gnuplot Cspline', \
     'Cspline.xy.txt' using 1:2 w p title 'Data Points'

