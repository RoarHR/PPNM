using System;

class Program{
static void Main(){
	char[] split_delimiters = {' ','\t','\n'};
	var split_options = StringSplitOptions.RemoveEmptyEntries;
	for( string line = Console.ReadLine(); line != null; line = Console.ReadLine() ){
		//Console.WriteLine($"line = {line}, Here1");
		var numbers = line.Split(split_delimiters,split_options);
		//Console.WriteLine($"numbers = {numbers}, Here2");
		foreach(var number in numbers){
			//Console.WriteLine($"line = {line}\t number={number}");
			double x = double.Parse(number);
			Console.Error.WriteLine($"{x:e4}\t{Math.Sin(x):e4}\t{Math.Cos(x):e4}");
		}
		}
} // Main
} // Program
