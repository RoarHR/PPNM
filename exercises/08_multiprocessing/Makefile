CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) # $(CODE) may be needed
TIME = time --output=$@ --append --format "$$nthreads %e %U"

all : Out.times.svg out.pitfalls.txt Makefile

Out.times.svg : out.times.txt Makefile
	echo '\
		set terminal svg background "white";\
		set out "$@";\
		set xlabel "Number of threads";\
		set ylabel "Real time [s]";\
		set title "Running time using multithreading";\
		plot "$<" with lp;\
		'| tee log-gpi | gnuplot

N = 3e8
out.times.txt : main.exe
	rm -f log.txt
	rm -f out.times.txt
	for nthreads in 1 2 3 4 5 6 7 8; do\
		$(TIME) mono main.exe -nthreads:$$nthreads -nterms:$N >> log.txt;\
	done

out.pitfalls.txt : main-for.exe main-threadLocal.exe Makefile
	echo "User time [s] Real time[s]\nUsing Parallel.For\n" > $@
	$(TIME) mono main-for.exe -nterms:$N >> $@
	echo "\n\nUsing ThreadLocal" >> $@
	$(TIME) mono main-threadLocal.exe -nterms:$N >> $@

main.exe : main.cs ; $(MKEXE)
main-for.exe : main-for.cs ; $(MKEXE)
main-threadLocal.exe : main-threadLocal.cs ; $(MKEXE)


clean:
	$(RM) *.exe *.dll [Oo]ut* [Ll]og*
