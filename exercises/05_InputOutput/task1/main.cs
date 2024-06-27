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
}
}
