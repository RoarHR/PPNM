CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^)) 
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS)

all: Out.txt ErfPlot.svg

ErfPlot.svg: ErfPlot.gnuplot Out.Erf.txt err_tabulated.txt
	gnuplot $<

Out.txt Out.Erf.txt: main.exe integration_test.py
	mono $< > $@ 2> log.txt
	python integration_test.py >> $@
	cat $@

main.exe: main.cs
	$(MKEXE)

clean:
	$(RM) *.exe *.dll *.log log* *.svg Out*

.PHONY: all clean
