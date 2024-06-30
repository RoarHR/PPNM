set terminal svg size 800,600 enhanced font 'Verdana,10' background "white"
set output 'Fit.svg'

# Set titles and labels
set title "Plot"
set xlabel "x"
set ylabel "f(x)"

# Define the function
f(x) = cos(5 * x - 1) * exp(- x * x)

# Plot the function
plot "Out.Fit.txt" using 1:2 w l title 'Neural network', \
     f(x) title 'cos(5 * x - 1) * exp(- x * x)' lw 2
