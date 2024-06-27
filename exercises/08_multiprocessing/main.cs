using System;
using System.Threading;

class Program{
public class harmdata{ public int a,b; public double sum;}

static void harm(object obj){
	harmdata d = (harmdata)obj;
	d.sum=0;
	for(int i=d.a; i<= d.b; i++){
		d.sum+= 1.0/i;
		//Console.Write($"harm: sum={d.sum}\n");
	}
}

public static int Main(string[] args){
	int nterms=(int)1e7,nthreads=1;
	foreach(string arg in args){
		Console.Write($"arg -> {arg}\n");
		var words = arg.Split(':');
		if(words[0]=="-nterms")nterms= (int)double.Parse(words[1]);
		if(words[0]=="-nthreads")nthreads= (int)double.Parse(words[1]);
	}
	Console.Write($"Main: nterms={nterms} nthreads={nthreads}\n");
	harmdata[] data = new harmdata[nthreads];
	int chunk=nterms/nthreads;
	for(int i=0;i<nthreads;i++){
		data[i] = new harmdata();
		data[i].a=i*chunk+1;
		data[i].b=data[i].a+chunk;
	}
	data[nthreads-1].b=nterms;
	var threads = new System.Threading.Thread[nthreads];
	Console.Write($"Main: Starting threads...\n");
	for(int i=0;i<nthreads;i++){
		threads[i] = new System.Threading.Thread(harm);
		threads[i].Start(data[i]);
	}

	foreach(var thread in threads) thread.Join();

	double total=0;
	foreach(var p in data) total+=p.sum;
	Console.WriteLine($"Total Harmonic Sum: {total}\n");
	return 0;
} // Main
} // Program
