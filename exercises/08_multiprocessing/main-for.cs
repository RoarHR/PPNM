using System;
using System.Threading.Tasks;

class Program {
static void Main(string[] args) {
	double sum = 0;

	int nterms=(int)1e7;
	foreach(string arg in args){
		Console.Write($"arg -> {arg}\n");
		var words = arg.Split(':');
		if(words[0]=="-nterms")nterms= (int)double.Parse(words[1]);
		//if(words[0]=="-nthreads")nthreads= (int)double.Parse(words[1]);
	}
	Console.Write($"Main: nterms={nterms}");
	
	System.Threading.Tasks.Parallel.For( 1, nterms+1, (int i) => sum+=1.0/i );

	Console.Write($"Total Harmonic Sum: {sum}\n");
} // Main
} // Program
