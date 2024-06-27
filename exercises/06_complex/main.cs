using System;

public class Program{
static void Main(){
	int FailedTests = 0;
	double pi = 3.14159265358979323846;
	complex a = new complex(-1,0);
	complex sqa = cmath.sqrt(a);
	Console.WriteLine("\nTesting sqrt(-1) = +/- i...");
	complex i = new complex(0,1);
	if (complex.approx(sqa, i)){
		Console.WriteLine("Passed, sqrt(-1)=+i");}
	else if (complex.approx(sqa, -i)){
		Console.WriteLine("Passed, sqrt(-1)=-i");}
	else {Console.WriteLine("Failed"); FailedTests++;}

	Console.WriteLine("\nTesting sqrt(i) = 1/sqrt(2) + i/sqrt(2)...");
	complex test = cmath.sqrt(i);
	complex result = new complex(1/System.Math.Sqrt(2), 1/System.Math.Sqrt(2));
	if (complex.approx(test, result)){Console.WriteLine("Passed");}
	else {Console.WriteLine("Failed"); FailedTests++;}

	Console.WriteLine("\nTesting e^i = 0.54030230586813971740 0.84147098480789650665 i)");
	test = cmath.exp(i);
	result = new complex(0.54030230586813971740, 0.84147098480789650665);
	if (complex.approx(result, test)){Console.WriteLine("Passed");}
	else {Console.WriteLine("Failed"); FailedTests++;}

	Console.WriteLine("\nTesting e^pi i = -1");
	test = cmath.exp(i * pi);
	double dresult = -1;
	if (complex.approx(dresult, test)){Console.WriteLine("Passed");}
	else {Console.WriteLine("Failed"); FailedTests++;}

	Console.WriteLine("\nTesting i^i = e^(-pi/2)");
	test = cmath.pow(i, i);
	dresult = System.Math.Exp(-pi/2);
	if (complex.approx(dresult, test)){Console.WriteLine("Passed");}
	else {Console.WriteLine("Failed"); FailedTests++;}

	Console.WriteLine("\nTesting ln(i) = i pi/2");
	test = cmath.log(i);
	result = i * pi / 2;
	if (complex.approx(result, test)){Console.WriteLine("Passed");}
	else {Console.WriteLine("Failed"); FailedTests++;}

	Console.WriteLine("\nTesting sin(i pi) = 11.54873935725774837797 i");
	test = cmath.sin(i * pi);
	result = 11.54873935725774837797 * i;
	if (complex.approx(result, test)){Console.WriteLine("Passed");}
	else {Console.WriteLine("Failed"); FailedTests++;}

	if (FailedTests == 0){Console.WriteLine("\nPassed all tests");}
	else {Console.WriteLine($"\nFailed {FailedTests} tests");}
} // Main
} // Program
