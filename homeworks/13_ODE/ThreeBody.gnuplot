set terminal svg enhanced background "white"
set output 'ThreeBody.svg'

set xlabel 'x'
set ylabel 'y'

set title 'Three-Body Problem'

# Plot the trajectories of the three bodies
plot 'Out.ThreeBody.txt' using 8:9 with lines title 'Body 1', \
     'Out.ThreeBody.txt' using 10:11 with lines title 'Body 2', \
     'Out.ThreeBody.txt' using 12:13 with lines title 'Body 3'

