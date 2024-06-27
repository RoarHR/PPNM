set terminal svg enhanced background "white"
set output 'FricOscillator.svg'

set xlabel 'Time'
set ylabel 'Angle'

# Title
set title 'Oscillator with friction (as in scipy example)'

plot 'Out.FricOscillator.txt' using 1:2 w p pt 6 title '{/Symbol q}'

