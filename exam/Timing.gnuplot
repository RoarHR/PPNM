set terminal svg size 800,600 enhanced font 'Verdana,10' background "white"
set output 'Timing.svg'

# Set titles and labels
set title "Execution Time vs N"
set xlabel "N"
set ylabel "Time (seconds)"

# Set grid
set grid

# Plot the data
plot "Out.Timing.txt" using 1:2 with linespoints title "Execution Time" lw 2 pt 7

