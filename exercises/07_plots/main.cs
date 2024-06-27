using System;
using System.IO;
using static System.Math;

public class Program{
	static double erf(double x){
	/// single precision error function (Abramowitz and Stegun, from Wikipedia)
	if(x<0) return -erf(-x);
	double[] a={0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429};
	double t=1/(1+0.3275911*x);
	double sum=t*(a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*a[4]))));/* the right thing */
	return 1-sum*Exp(-x*x);
	}

	static double gamma(double x){
	///single precision gamma function (Gergo Nemes, from Wikipedia)
	if(x<0)return PI/Sin(PI*x)/gamma(1-x);
	if(x<9)return gamma(x+1)/x;
	double lngamma=x*Log(x+1/(12*x-1/x/10))-x+Log(2*PI/x)/2;
	return Exp(lngamma);
	}

	static double lngamma(double x){
	if(x<=0) throw new ArgumentException("lngamma: x<=0");
	if(x<9) return lngamma(x+1)-Log(x);
	return x*Log(x+1/(12*x-1/x/10))-x+Log(2*PI/x)/2;
	}

	static int factorial(int x){
		int result = 1;
		for (int i = 1; i <= x; i++){result = result * i;}
		return result;
	}
	
static void Main(string[] args){
	string out1 = null;
	string out2 = null;
	string out3 = null;
	foreach (string arg in args){
		if (arg.Split(":")[0] == "-out1"){out1 = arg.Split(":")[1];}
		else if (arg.Split(":")[0] == "-out2"){out2 = arg.Split(":")[1];}
		else if (arg.Split(":")[0] == "-out3"){out3 = arg.Split(":")[1];}
	}

	// Error function
	using (StreamWriter writer = new StreamWriter(out1)){
		for (double x = 0.001; x < 3.5; x = x + 0.005){
			writer.WriteLine($"{x:e2}\t{erf(x):e5}");
		}
	}

	// Gamma function
	using (StreamWriter writer = new StreamWriter(out2)){
		for (double x = 0.1; x < 6.2; x = x + 0.1){
			writer.WriteLine($"{x:e2}\t{gamma(x):e5}");
		}
	}
	using (StreamWriter writer = new StreamWriter("factorials_tabulated.txt")){
		for (int x = 1; x <= 6; x++){
			writer.WriteLine($"{x}\t{factorial(x - 1)}");
		}
	}

	// ln gamma function
	using (StreamWriter writer = new StreamWriter(out3)){
		for (double x = 0.1; x < 6.2; x = x + 0.1){
			writer.WriteLine($"{x:e2}\t{lngamma(x):e5}");
		}
	}
	using (StreamWriter writer = new StreamWriter("lnfactorials_tabulated.txt")){
		for (int x = 1; x <= 6; x++){
			writer.WriteLine($"{x}\t{Log(factorial(x - 1))}");
		}
	}
	
} // Main

} // Program
