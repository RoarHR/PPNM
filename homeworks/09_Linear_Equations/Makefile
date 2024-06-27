CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)
TIME = time --output=$@ --append --format "$$matrixsize %e"

all: Out.txt Out.svg

Out.txt: main.exe
	mono $< >$@
	cat $@

Out.svg: Out.times.txt fit.gnuplot
	$(RM) $@
	gnuplot -e "input_file='$<'; output_file='$@'" fit.gnuplot

Out.times.txt: main.exe
	$(RM) log.txt $@
	for matrixsize in $(shell seq 10 60 900); do \
		echo "Decomposing matrix size $$matrixsize"; \
		$(TIME) mono main.exe -matrixsize:$$matrixsize >> timing.log; \
	done

main.exe: main.cs linalg.dll
	$(MKEXE)

linalg.dll: matrix.cs vec.cs
	$(MKLIB)

clean:
	$(RM) *.exe *.dll *.log *.svg Out*

