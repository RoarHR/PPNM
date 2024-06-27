set terminal svg enhanced background "white"
set output 'HarmOscillator.svg'

set xlabel 'time'
set ylabel 'x'

# Title
set title 'First ODE test'

plot 'Out.HarmOscillator.txt' using 1:2 w p pt 6 title 'Position'


