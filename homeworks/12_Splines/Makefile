CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS)

all: Out.txt Linterp.svg Qspline.svg Cspline.svg

Linterp.svg: Linterp.gnuplot Out.Linterp.xy.txt Out.Linterp.interp.txt
	gnuplot $<

Qspline.svg: Qspline.gnuplot Out.Qspline.xy.txt Out.Qspline.spline.txt
	gnuplot $<

Cspline.svg: Cspline.gnuplot Cspline.xy.txt Out.Cspline.spline.txt
	gnuplot $<

Out.txt Out.Linterp.xy.txt Out.Linterp.interp.txt Out.Qspline.xy.txt Out.Qspline.spline.txt Out.Cspline.spline.txt: main.exe Cspline.xy.txt
	mono $< > $@
	cat $@

main.exe: main.cs linalg.dll
	$(MKEXE)

linalg.dll: ../09_Linear_Equations/matrix.cs ../09_Linear_Equations/vec.cs 
	$(MKLIB) $^

clean:
	$(RM) *.exe *.dll *.log *.svg Out*
