CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS)

all: Out.txt Fit.svg

Fit.svg: Fit.gnuplot Out.Fit.txt
	gnuplot $<

Out.txt Out.Fit.txt: main.exe Makefile
	mono $< > $@
	cat $@

main.exe: main.cs linalg.dll
	$(MKEXE)

linalg.dll: ../09_Linear_Equations/matrix.cs ../09_Linear_Equations/vec.cs ../09_Linear_Equations/main.cs ../19_Minimization/main.cs
	$(MKLIB) $^

clean:
	$(RM) *.exe *.dll *.log *.svg Out*

.PHONY: all clean
