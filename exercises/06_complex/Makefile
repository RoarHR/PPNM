Out.txt : main.exe
	mono $< > $@
	cat $@

main.exe : main.cs cmath.dll
	mcs -reference:cmath.dll -target:exe -out:main.exe main.cs

cmath.dll : cmath.cs complex.cs
	mcs -target:library -out:$@ $^

clean:
	$(RM) *.exe *.dll [Oo]ut.txt
