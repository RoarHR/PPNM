CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS)

all: Out.txt ErrorPlot.svg

Out.txt: main.exe Makefile
	mono $< > $@
	cat $@

Out.ErrorPlot.txt: main.exe Makefile
	rm -f $@
	for N in 100 125 158 200 251 316 398 501 631 794 1000 1259 1585 1995 2512 3162 3981 5011 6309 7943 10000 12589 15849 19953; do \
		mono $< $$N; \
	done

ErrorPlot.svg: ErrorPlot.gnuplot Out.ErrorPlot.txt
	gnuplot $< 2> gnuplot.log

main.exe: main.cs linalg.dll
	$(MKEXE)

linalg.dll: ../09_Linear_Equations/matrix.cs ../09_Linear_Equations/vec.cs
	$(MKLIB) $^

clean:
	$(RM) *.exe *.dll *.log *.svg Out* ErrorPlot.svg

.PHONY: all clean
