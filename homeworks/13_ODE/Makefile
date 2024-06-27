CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS)

all: Out.txt HarmOscillator.svg FricOscillator.svg EquatorialMotion_Circular.svg EquatorialMotion_Elliptical.svg EquatorialMotion_Precession.svg ThreeBody.svg

HarmOscillator.svg: HarmOscillator.gnuplot Out.HarmOscillator.txt
	gnuplot $<

FricOscillator.svg: FricOscillator.gnuplot Out.FricOscillator.txt
	gnuplot $<

EquatorialMotion_Circular.svg: EquatorialMotion.gnuplot Out.EquatorialMotion_Circular.txt
	gnuplot -e "input_file='Out.EquatorialMotion_Circular.txt'; output_file='$@'" $<

EquatorialMotion_Elliptical.svg: EquatorialMotion.gnuplot Out.EquatorialMotion_Elliptical.txt
	gnuplot -e "input_file='Out.EquatorialMotion_Elliptical.txt'; output_file='$@'" $<

EquatorialMotion_Precession.svg: EquatorialMotion.gnuplot Out.EquatorialMotion_Precession.txt
	gnuplot -e "input_file='Out.EquatorialMotion_Precession.txt'; output_file='$@'" $<

ThreeBody.svg: ThreeBody.gnuplot Out.ThreeBody.txt
	gnuplot $<

Out.txt Out.HarmOscillator.txt Out.FricOscillator.txt Out.EquatorialMotion_Circular.txt Out.EquatorialMotion_Elliptical.txt Out.EquatorialMotion_Precession.txt Out.ThreeBody.txt: main.exe 
	mono $< > $@ 2> log.txt
	cat $@

main.exe: main.cs linalg.dll
	$(MKEXE)

linalg.dll: ../09_Linear_Equations/matrix.cs ../09_Linear_Equations/vec.cs ../12_Splines/main.cs
	$(MKLIB) $^

clean:
	$(RM) *.exe *.dll *.log *.svg Out*
