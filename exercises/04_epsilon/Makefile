Out.txt : main.exe
	mono $< > $@
	cat Out.txt

main.exe : main.cs
	mcs -target:exe -out:$@ $<

clean:
	rm -f main.exe Out.txt

