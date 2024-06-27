using System;
using System.IO;
using System.Collections.Generic;

class Program{
static void Main(string[] args){
	string input = null;
	string output = null;
	foreach(string arg in args){
		if (arg.Split(":")[0] == "-input"){input = arg.Split(":")[1];}
		else if(arg.Split(":")[0] == "-output"){output = arg.Split(":")[1];}
	}
	List<string> numbers = new List<string>();
	using (StreamReader reader = new StreamReader(input)){
		for (string line=reader.ReadLine();line!=null;line=reader.ReadLine()){
			numbers.Add(line);
		}
	}
	using (StreamWriter writer = new StreamWriter(output,append:true)){
		foreach (string number in numbers){
			double x=double.Parse(number);
			writer.WriteLine($"{x:e4}\t{Math.Sin(x):e4}\t{Math.Cos(x):e4}");
		}
	}
} // Main
} // Program
