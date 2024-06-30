CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS)

all: Out.txt HiggsBoson.svg

HiggsBoson.svg: HiggsBoson.gnuplot Out.HiggsBoson.txt CERN.data
	gnuplot $<

Out.txt Out.HiggsBoson.txt: main.exe Makefile CERN.data
	mono $< < CERN.data > $@
	cat $@

main.exe: main.cs linalg.dll
	$(MKEXE)

linalg.dll: ../09_Linear_Equations/matrix.cs ../09_Linear_Equations/vec.cs ../09_Linear_Equations/main.cs
	$(MKLIB) $^

clean:
	$(RM) *.exe *.dll *.log *.svg Out*

.PHONY: all clean
