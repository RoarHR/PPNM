using System;

class Program{
static void Main(string[] args){
	foreach(var arg in args){
		if(arg.Split(":")[0] == "-numbers"){
			string[] numbers;
			numbers = arg.Split(":")[1].Split(",");
			foreach(string number in numbers){
				double x = double.Parse(number);
				Console.WriteLine($"{x:e4}\t{Math.Sin(x):e4}\t{Math.Cos(x):e4}");
			}
		}
	}
}
}
