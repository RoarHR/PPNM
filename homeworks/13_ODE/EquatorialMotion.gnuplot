set terminal svg enhanced background "white"
set output output_file

set xlabel 'x'
set ylabel 'y'
set title 'Equatorial Motion'

plot input_file using (1/$2)*cos($1):(1/$2)*sin($1) with lines title 'Orbit'

