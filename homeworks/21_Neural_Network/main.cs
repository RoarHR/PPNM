using System;
using System.Collections.Generic;
using System.IO;
using static System.Math;

public class NelderMead {
	public static Vec Optimize(Func<Vec, double> costFunction, Vec initialGuess, double alpha = 1.0, double gamma = 2.0, double rho = 0.5, double sigma = 0.5, int maxfev = 1000, double tolerance = 1e-6) {
		int n = initialGuess.Length;
		List<Vec> simplex = new List<Vec>(n + 1);
		simplex.Add(initialGuess.Copy());

		// Initialize the simplex with perturbations
		for (int i = 0; i < n; i++) {
			Vec point = initialGuess.Copy();
			point[i] += 0.01; // Small perturbation
			simplex.Add(point);
		}

		int iterations = 0;
		while (iterations < maxfev) {
			// Debugging output
			simplex.Sort((a, b) => costFunction(a).CompareTo(costFunction(b)));
			Vec lowest = simplex[0];
			Vec highest = simplex[simplex.Count - 1];
			Vec secondHighest = simplex[simplex.Count - 2];

			// Calculate the centroid
			Vec centroid = new Vec(n);
			foreach (var point in simplex) {
				if (!point.Equals(highest))
					centroid += point;
			}
			centroid /= n;

			// Reflection
			Vec reflected = centroid + alpha * (centroid - highest);
			double reflectedCost = costFunction(reflected);

			if (reflectedCost < costFunction(lowest)) {
				// Expansion
				Vec expanded = centroid + gamma * (reflected - centroid);
				if (costFunction(expanded) < reflectedCost) {
					simplex[simplex.Count - 1] = expanded;
				} else {
					simplex[simplex.Count - 1] = reflected;
				}
			} else if (reflectedCost < costFunction(secondHighest)) {
				simplex[simplex.Count - 1] = reflected;
			} else {
				// Contraction
				Vec contracted = centroid + rho * (highest - centroid);
				if (costFunction(contracted) < costFunction(highest)) {
					simplex[simplex.Count - 1] = contracted;
				} else {
					// Reduction
					for (int i = 1; i < simplex.Count; i++) {
						simplex[i] = lowest + sigma * (simplex[i] - lowest);
					}
				}
			}

			// Check convergence
			double size = 0;
			foreach (var point in simplex)
				size += (point - centroid).Norm();

			if (size < tolerance)
				break;

			iterations++;
		}

		return simplex[0];
	}
}

public class Ann {
	private int n; // number of hidden neurons
	private Func<double, double> f; // activation function
	private Func<double,double> df; 	/* derivative of function */
	private Func<double,double> ddf; 	/* second derivative of function */
	private Func<double,double> If;		/* antiderivative of function */
	public Vec p; // network parameters

	public Ann(int n, Func<double, double> f = null, Func<double,double> df = null, Func<double,double> ddf = null, Func<double, double> If = null, Vec p0 = null) {
		this.n = n;
		this.f = f ?? (x => x * Math.Exp(-x * x)); // Default to Gaussian wavelet
		this.df = df ?? (x => (1 - 2 * x * x) * Math.Exp(-x * x));
		this.ddf = ddf ?? (x => 2 * Math.Exp(-x * x) * x * (2 * x * x - 3));
		this.If = If ?? (x => -  Math.Exp(-x * x) / 2);
		if (p0 == null) {
			p0 = new Vec(3 * n);
			p0.FillRandom(-0.1, 0.1); // Random initialization in the interval [0, 1]
		}
		this.p = p0;
	}

	public double Prompt(double x) {
		double sum = 0;
		for (int i = 0; i < n; i++) {
			double ai = p[3 * i];
			double bi = p[3 * i + 1];
			double wi = p[3 * i + 2];
			sum += f((x - ai) / bi) * wi;
		}
		return sum;
	}

	public double Cost(Vec xs, Vec ys) {
		double cost = 0;
		for (int i = 0; i < xs.Length; i++) {
			double diff = Prompt(xs[i]) - ys[i];
			//Console.WriteLine($"{Prompt(xs[i])}");
			cost += diff * diff;
		}
		return cost;
	}

	public void Train(Vec xs, Vec ys) {
		Func<Vec, double> costFunction = p => {
			this.p = p;
			return Cost(xs, ys);
		};
		this.p = NelderMead.Optimize(costFunction, p);
	}

	public double Derivative(double x) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		double ai = p[3 * i];
		double bi = p[3 * i + 1];
		double wi = p[3 * i + 2];
		double z = (x - ai) / bi;
		sum += df(z) * wi / bi;
		}
	return sum;
	}

	public double SecondDerivative(double x) {
		double sum = 0;
		for (int i = 0; i < n; i++) {
			double ai = p[3 * i];
			double bi = p[3 * i + 1];
			double wi = p[3 * i + 2];
			double z = (x - ai) / bi;
			sum += ddf(z) * wi / (bi * bi);
		}
		return sum;
	}

	public double Integral(double x) {
		double sum = 0;
		for (int i = 0; i < n; i++) {
			double ai = p[3 * i];
			double bi = p[3 * i + 1];
			double wi = p[3 * i + 2];
			double z = (x - ai) / bi;
			sum += If(z) * wi * bi;
		}
		return sum;
	}

}

public static class Program {
	public static void Main(string[] args) {
		Func<double, double> gaussianWavelet = x => x * Math.Exp(-x * x);
		// Func<double, double> gaussian = x => Math.Exp(-x * x);
		// Func<double, double> wavelet = x => Math.Cos(5 * x) * Math.Exp(-x * x);

		Func<double, double> activationFunction = gaussianWavelet;

		Func<double, double> targetFunction = x => Math.Cos(5 * x - 1) * Math.Exp(-x * x);

		int sampleSize = 100;
		Vec xs = new Vec(sampleSize);
		Vec ys = new Vec(sampleSize);
		double a = -1.5;
		double b = 1.5;

		double step = (b - a) / (sampleSize - 1);
		for (int i = 0; i < sampleSize; i++) {
			double x = a + i * step;
			xs[i] = x;
			ys[i] = targetFunction(x);
		}

		int hiddenNeurons = 5;
		Ann ann = new Ann(hiddenNeurons, activationFunction);
		ann.Train(xs, ys);


		double[] testPoints = new double[] { -1.2, 0.0, 0.6 };
		foreach (double x in testPoints) {
			double targetValue = targetFunction(x);
			double networkOutput = ann.Prompt(x);
			Console.WriteLine($"x = {x}, target = {targetValue}, network output = {networkOutput}");
		}

		using (StreamWriter sw = File.CreateText("Out.Fit.txt")) {
			for (int i = 0; i < 1000; i++) {
				double x = (double)i / 1000 * (b - a) + a;
				sw.WriteLine($"{x} {ann.Prompt(x)}");
			}
		}

		// Sample points to test derivatives and integrals
		foreach (double x in testPoints) {
			double networkOutput = ann.Prompt(x);
			double firstDerivative = ann.Derivative(x);
			double secondDerivative = ann.SecondDerivative(x);
			double integral = ann.Integral(x);

			Console.WriteLine($"x = {x}, F(x) = {networkOutput}, F'(x) = {firstDerivative}, F''(x) = {secondDerivative}, ∫F(x)dx = {integral}");
		}
	}
}
