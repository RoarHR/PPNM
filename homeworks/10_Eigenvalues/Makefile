CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

all: Out.txt eigenfunctions.svg rmax_convergence.svg dr_convergence.svg

Out.txt: main.exe
	mono $< >$@
	cat $@

rmax = 10
dr = 0.1

eigenfunctions.svg: eigenfunctions.gnuplot eigenfunctions_radii.txt
	gnuplot -e "input_file='eigenfunctions_radii.txt'; output_file='$@'" $<

eigenfunctions_radii.txt: eigenfunctions.txt radii.txt
	paste radii.txt eigenfunctions.txt > $@

radii.txt: Makefile
	@awk 'BEGIN {rmax=$(rmax); dr=$(dr); for (i=1; i<rmax/dr; i++) printf "%.4f\n", i*dr; }' > $@

eigenfunctions.txt: main.exe Makefile
	mono $< -rmax:$(rmax) -dr:$(dr) -printfunctions:yep > $@

rmax_convergence.svg: rmax_convergence.gnuplot rmax_convergence.txt
	gnuplot -e "input_file='rmax_convergence.txt'; output_file='$@'" $<

rmax_convergence.txt: main.exe Makefile
	rm -f $@
	for rmax in $(shell seq 1 0.3 12); do \
		mono $< -rmax:$$rmax -dr:$(dr) >> $@; \
	done

dr_convergence.svg: dr_convergence.gnuplot dr_convergence.txt
	gnuplot -e "input_file='dr_convergence.txt'; output_file='$@'" $<

dr_convergence.txt: main.exe Makefile
	rm -f $@
	for dr in $(shell seq 0.02 0.02 0.3); do \
		mono $< -rmax:$(rmax) -dr:$$dr >> $@; \
	done

main.exe: main.cs linalg.dll
	$(MKEXE)

linalg.dll: ../09_Linear_Equations/matrix.cs ../09_Linear_Equations/vec.cs
	$(MKLIB) $^

clean:
	$(RM) *.exe *.dll *.log *.svg Out* eigenfunctions_radii.txt radii.txt eigenfunctions.txt rmax_convergence*.txt dr_convergence*.txt

.PHONY: all clean
