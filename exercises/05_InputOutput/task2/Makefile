Out.txt : main.exe input.txt
	echo 1 2 3 4 5 > input.txt
	mono $< -numbers:0.1,0.4,0.6,1 < input.txt > $@
	cat Out.txt

main.exe : main.cs
	mcs -target:exe -out:$@ $<

clean:
	rm -f main.exe Out.txt

