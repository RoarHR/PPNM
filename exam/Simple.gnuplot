set terminal svg size 800,600 enhanced font 'Verdana,10' background "white"
set output 'Simple.svg'

# Set titles and labels
set title "Plot of f(x) = c[0] * sin(c[1] * x + c[3]) * cos(c[2] * x + c[3])"
set xlabel "x"
set ylabel "f(x)"

load 'Out.fitparameters.txt'

# Define the function
f(x) = fit0 * sin(fit1 * x + fit3) * cos(fit2 * x + fit3)

# Plot the function
plot "Out.Generated_points.data" using 1:2 w p title 'Data' pt 7, \
     f(x) title 'f(x) = 1 * sin(0.8 * x + 0.5) * cos(2 * x + 0.5)' lw 2

