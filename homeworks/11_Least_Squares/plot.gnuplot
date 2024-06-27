set terminal svg enhanced background "white"
set output 'plot.svg'

set xlabel 'Time [days]'
set ylabel 'Activity'

# Title
set title 'Radioactive decay of ThX'

# Read coefficients from file
load 'Out.coefficients.txt'

# Define the exponential decay function
f(t) = a * exp(-lambda * t)
fm(t) = (a-da) * exp(- (lambda + dlambda) * t)
fp(t) = (a+da) * exp(- (lambda - dlambda) * t)

# Plot the data with error bars and the exponential decay function
plot f(x) title sprintf('Best fit: y = %.2f * exp(-%.2f * t)', a, lambda) lw 2, \
     fm(x) title '1 {/Symbol s} lower bound' lw 2, \
     fp(x) title '1 {/Symbol s} upper bound' lw 2, \
     'decay.data' using 1:2:3 with yerrorbars title 'Data with Error Bars'
