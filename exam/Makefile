CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS)
TIME = time --output=$@ --append --format "%e"

all: Out.txt LAWD_37.svg Simple.svg Timing.svg

LAWD_37.svg: LAWD_37.gnuplot Out.LAWD_37_fit.dat Out.LAWD_37.dat
	gnuplot $<

Simple.svg: Simple.gnuplot Out.Generated_points.data Out.fitparameters.txt
	gnuplot $<
	
Timing.svg: Timing.gnuplot Out.Timing.txt
	gnuplot $<

Out.txt Out.LAWD_37_fit.dat Out.LAWD_37.dat Out.Generated_points.data Out.fitparameters.txt: main.exe LAWD_37.data
	mono $< > $@
	cat Out.txt

Out.Timing.txt: main.exe
	rm -f Out.Timing.txt
	for N in 1000 10000 20000 30000 40000 50000 60000 70000; do \
		echo -n "$$N " >> Out.Timing.txt; \
		$(TIME) mono $< -N:$$N >> Out.Timing.txt 2>&1; \
	done

main.exe: main.cs linalg.dll
	$(MKEXE)

linalg.dll: ../homeworks/09_Linear_Equations/matrix.cs ../homeworks/09_Linear_Equations/vec.cs ../homeworks/09_Linear_Equations/main.cs
	$(MKLIB) $^

clean:
	$(RM) *.exe *.dll *.log *.svg Out*

.PHONY: all clean

