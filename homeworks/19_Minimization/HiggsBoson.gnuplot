set terminal svg enhanced background "white"
set output 'HiggsBoson.svg'

# Titles and labels
set title "Higgs Boson Signal"
set xlabel "Energy E [GeV]"
set ylabel "Signal σ(E) [units]"

# Plot data points with error bars and the fitted function
plot 'CERN.data' using 1:2:3 with errorbars title "Experimental Data", \
     'Out.HiggsBoson.txt' using 1:2 with lines title "Fitted Function"

