set terminal svg enhanced background "white"
set output output_file

set xlabel "dr"
set ylabel "epsilon_0"
set title "Convergence of epsilon_0 with varying dr"

plot input_file using 2:3 title "epsilon_0" with linespoints

