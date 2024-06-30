set terminal svg size 800,600 enhanced font 'Verdana,10' background "white"
set output 'LAWD_37.svg'

set title "3D Plot of Original Data and Fitted Gaussian"
set xlabel "X-axis"
set ylabel "Y-axis"
set zlabel "Intensity"

set dgrid3d 30,30
set hidden3d
set style data lines

# Add transparency to the fitted Gaussian
set style fill transparent solid 0.5

# Rotate for better visibility
set view 60, 120

splot 'Out.LAWD_37.dat' with points pointtype 7 pointsize 0.5 linecolor rgb "blue" title 'Original Data', \
      'Out.LAWD_37_fit.dat' with lines linecolor rgb "red" title 'Fitted Gaussian'
