CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS)

all: Out.txt plot.svg

plot.svg: plot.gnuplot decay.data Out.coefficients.txt
	gnuplot $<

Out.txt Out.coefficients.txt: main.exe decay.data
	echo "The table value for the halflife of 224Ra is 3.631 +/- 2 days" > $@
	mono $< < decay.data >> $@ 2> Out.coefficients.txt
	echo "The value from the Rutherford and Soddy experiment is thus outside the uncertainty of the table value. Perhaps the uncertainties were underestimated." >> $@
	cat $@

main.exe: main.cs linalg.dll
	$(MKEXE)

linalg.dll: ../09_Linear_Equations/matrix.cs ../09_Linear_Equations/vec.cs ../09_Linear_Equations/main.cs
	$(MKLIB) $^

clean:
	$(RM) *.exe *.dll *.log *.svg Out*
