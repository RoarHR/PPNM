Out.txt: main1.exe main2.exe main3.exe
	echo "Task 1:\n" > $@
	mono main1.exe -numbers:0.1,0.3,0.5,0.7,0.9 >> $@
	
	echo "\nTask 2:" >> $@
	echo "0.2 0.4	0.6\n0.8" > input.txt
	mono main2.exe < input.txt 2>> $@
	
	echo "\nTask 3:" >> $@
	mono main3.exe -input:input3.txt -output:$@
	cat $@

main1.exe: main1.cs
	mcs -target:exe -out:$@ $<

main2.exe: main2.cs
	mcs -target:exe -out:$@ $<

main3.exe: main3.cs
	mcs -target:exe -out:$@ $<

clean:
	$(RM) *.exe [Oo]ut.txt
