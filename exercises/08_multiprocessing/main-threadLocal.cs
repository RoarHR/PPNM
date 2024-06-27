using System;
using System.Linq;
using System.Threading.Tasks;

class Program {
static void Main(string[] args) {
	int nterms = (int)1e7;

	foreach (string arg in args) {
		var words = arg.Split(':');
		if (words[0] == "-nterms") nterms = (int)double.Parse(words[1]);
	}
	Console.WriteLine($"Main: nterms={nterms}\n");

	var sum = new System.Threading.ThreadLocal<double>(() => 0, trackAllValues: true);
	System.Threading.Tasks.Parallel.For(1, nterms + 1, (int i) => sum.Value += 1.0 / i);
	double totalsum = sum.Values.Sum();

	Console.WriteLine($"Total Harmonic Sum: {totalsum}\n");
    }
}
