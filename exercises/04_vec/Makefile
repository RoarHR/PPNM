Out.txt : main.exe
	mono $< > $@
	cat Out.txt

main.exe : main.cs vec.dll
	mcs -target:exe -out:$@ -reference:vec.dll $<

vec.dll: vec.cs
	mcs -target:library -out:$@ $<

clean:
	rm -f main.exe Out.txt vec.dll

