Out.txt : main.exe input.txt
	mono $< < input.txt > $@
	cat $@

main.exe : main.cs
	mcs $< -target:exe -out:$@
