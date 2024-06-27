using System;

class Program{
static void Main(string[] args){
	foreach(var arg in args){
		//Console.WriteLine(arg.Split(":")[0]);
		if(arg.Split(":")[0] == "-numbers"){
			//Console.WriteLine("True");
			string[] numbers;
			numbers = arg.Split(":")[1].Split(",");
			foreach(string number in numbers){
				double x = double.Parse(number);
				Console.WriteLine($"{x:e6}\t{Math.Sin(x):e6}\t{Math.Cos(x):e6}");
			}
		}
	}
	char[] split_delimiters = {' ','\t','\n'};
	var split_options = StringSplitOptions.RemoveEmptyEntries;
	for( string line = Console.ReadLine(); line != null; line = Console.ReadLine() ){
		Console.WriteLine($"line = {line}, Here1");
		//var numbers = line.Split(split_delimiters,split_options);
		//Console.WriteLine($"numbers = {numbers}, Here2");
		/*foreach(var number in numbers){
			Console.WriteLine($"line = {line}\t number={number}");
			double x = double.Parse(number);
			Console.Error.WriteLine($"{x} {Math.Sin(x)} {Math.Cos(x)}");
			}*/
		}
}
}
