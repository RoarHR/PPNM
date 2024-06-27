set datafile separator " "
set terminal svg enhanced background "white"
set output output_file

f(N) = a * N**3
a = 1.0
fit f(x) input_file using 1:2 via a

set xlabel "Matrix size N"
set ylabel "Time (seconds)"
set title "QR Factorization Timing"
set key at screen 0.4, 0.8
plot input_file using 1:2 title 'Data' with points pointtype 7, \
	f(x) title sprintf("Fit: a*N^3, a=%.3e", a) with lines

