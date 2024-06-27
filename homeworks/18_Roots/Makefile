CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS)

default_rmin = 0.01
default_rmax = 8
default_acc = 0.01
default_eps = 0.01

rmin_values = 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2
rmax_values = 4 6 8 10 12 14 16 18 20
acc_values = 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1
eps_values = 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1

all: Out.txt WaveFunction.svg Convergence.rmin.svg Convergence.rmax.svg Convergence.acc.svg Convergence.eps.svg

WaveFunction.svg: WaveFunction.gnuplot Out.WaveFunction.txt
	gnuplot $< 2> gnuplot.log

Out.txt Out.WaveFunction.txt: main.exe
	mono $< > $@
	cat $@

Convergence.rmin.svg: Convergence.rmin.gnuplot Out.Convergence.rmin.txt
	gnuplot $< 2> gnuplot.log

Convergence.rmax.svg: Convergence.rmax.gnuplot Out.Convergence.rmax.txt
	gnuplot $< 2> gnuplot.log

Convergence.acc.svg: Convergence.acc.gnuplot Out.Convergence.acc.txt
	gnuplot $< 2> gnuplot.log

Convergence.eps.svg: Convergence.eps.gnuplot Out.Convergence.eps.txt
	gnuplot $< 2> gnuplot.log

Out.Convergence.rmin.txt: main.exe
	$(RM) $@
	for rmin in $(rmin_values); do \
		mono $< -filename:Out.Convergence.rmin.txt -rmin:$$rmin -rmax:$(default_rmax) -acc:$(default_acc) -eps:$(default_eps); \
	done

Out.Convergence.rmax.txt: main.exe
	$(RM) $@
	for rmax in $(rmax_values); do \
		mono $< -filename:$@ -rmin:$(default_rmin) -rmax:$$rmax -acc:$(default_acc) -eps:$(default_eps); \
	done

Out.Convergence.acc.txt: main.exe
	$(RM) $@
	for acc in $(acc_values); do \
		mono $< -filename:$@ -rmin:$(default_rmin) -rmax:$(default_rmax) -acc:$$acc -eps:$(default_eps); \
	done

Out.Convergence.eps.txt: main.exe
	$(RM) $@
	for eps in $(eps_values); do \
		mono $< -filename:$@ -rmin:$(default_rmin) -rmax:$(default_rmax) -acc:$(default_acc) -eps:$$eps; \
	done

main.exe: main.cs linalg.dll
	$(MKEXE)

linalg.dll: ../09_Linear_Equations/matrix.cs ../09_Linear_Equations/vec.cs ../09_Linear_Equations/main.cs ../13_ODE/main.cs ../12_Splines/main.cs
	$(MKLIB) $^

clean:
	$(RM) *.exe *.dll *.log log.* *.svg Out* Out.Convergence.*
