set terminal svg enhanced background "white"
set output output_file

set xlabel "rmax"
set ylabel "epsilon_0"
set title "Convergence of epsilon_0 with varying rmax"

plot input_file using 1:3 title "epsilon_0" with linespoints

