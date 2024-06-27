set terminal svg enhanced background "white"
set output output_file

set xlabel "Radius (r/a0)"
set ylabel "Eigenfunctions"
set title "First Four Eigenfunctions of the Hydrogen Atom"

# Define constants
Z = 1   # Nuclear charge for hydrogen
a0 = 1  # Bohr radius

# Define the first four analytic wavefunctions
R21(r) = (1/sqrt(3.0)) * (Z/(2.0*a0))**(3.0/2.0) * (Z*r/a0) * exp(-Z*r/(2.0*a0))
R20(r) = 2.0 * (Z/(2.0*a0))**(3.0/2.0) * (1.0 - Z*r/(2.0*a0)) * exp(-Z*r/(2.0*a0))
R32(r) = (2.0*sqrt(2.0)/(27.0*sqrt(5.0))) * (Z/(3.0*a0))**(3.0/2.0) * (Z*r/a0)**2.0 * exp(-Z*r/(3.0*a0))
R31(r) = (4.0*sqrt(2.0)/3.0) * (Z/(3.0*a0))**(3.0/2.0) * (Z*r/a0) * (1.0 - Z*r/(6.0*a0)) * exp(-Z*r/(3.0*a0))

# Plot numerical and analytic wavefunctions
plot input_file using 1:2 title "Numerical Eigenfunction 1" with p, \
     input_file using 1:3 title "Numerical Eigenfunction 2" with p, \
     input_file using 1:4 title "Numerical Eigenfunction 3" with p, \
     input_file using 1:5 title "Numerical Eigenfunction 4" with p, \
     R21(x) title "Analytic R21" with lines, \
     R20(x) title "Analytic R20" with lines, \
     R32(x) title "Analytic R32" with lines, \
     R31(x) title "Analytic R31" with lines
