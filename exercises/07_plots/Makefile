all : errfunc.svg gammafunc.svg lngammafunc.svg

errfunc.svg: errfunc.txt err_tabulated.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "x" ;\
		set ylabel "erf(x)" ;\
		set autoscale ;\
		set offsets graph 0, 0.1, 0.1, 0 ;\
		set title "Error function" ;\
		plot "$<" using 1:2 w l title "Error function", "err_tabulated.txt" w p lt 6 lc 2 title "Tabulated values" ;\
		' | gnuplot


gammafunc.svg: gammafunc.txt factorials_tabulated.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "x" ;\
		set ylabel "gamma(x)" ;\
		set autoscale ;\
		set offsets graph 0, 0.1, 0.1, 0 ;\
		set title "Gamma function" ;\
		plot "$<" using 1:2 w l title "Gamma function", "factorials_tabulated.txt" w p lt 6 lc 2 title "Tabulated values" ;\
		' | gnuplot


lngammafunc.svg: lngammafunc.txt lnfactorials_tabulated.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "x" ;\
		set ylabel "ln(gamma(x))" ;\
		set autoscale ;\
		set offsets graph 0, 0.1, 0.1, 0 ;\
		set title "ln Gamma function" ;\
		plot "$<" using 1:2 w l title "ln Gamma function", "lnfactorials_tabulated.txt" w p lt 6 lc 2 title "Tabulated values" ;\
		' | gnuplot

errfunc.txt gammafunc.txt lngammafunc.txt : main.exe
	mono $< -out1:errfunc.txt -out2:gammafunc.txt -out3:lngammafunc.txt

main.exe : main.cs
	mcs -target:exe -out:$@ $<

clean :
	$(RM) *func.svg *func.txt main.exe
