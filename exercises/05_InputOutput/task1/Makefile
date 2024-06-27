Out.txt : main.exe
	mono $< -numbers:0.1,0.4,0.6,1 > $@
	cat Out.txt

main.exe : main.cs
	mcs -target:exe -out:$@ $<

clean:
	rm -f main.exe Out.txt

