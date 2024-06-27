set terminal svg enhanced background "white"
set output 'Wavefunction.svg'

# Set the labels and title
set title "Wavefunction Comparison"
set xlabel "r"
set ylabel "f(r)"

# Define the exact wavefunction
f0(r) = r * exp(-r)

# Plot the numerical wavefunction from the data file and the exact wavefunction
plot "Out.WaveFunction.txt" using 1:2 with lines title "Numerical Wavefunction", \
     f0(x) with lines title "Exact Wavefunction"

