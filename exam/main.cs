using System;
using System.Collections.Generic;
using System.IO;
using static System.Math;

public static class GaussNewton {
	public static (Matrix, Vec) Jacobian(List<Func<Vec, double>> r, Vec c) {
		// Calculates the Jacobian matrix and the vector r_i = r_i(c) at a point c for function that is a sum of squares
		int n = r.Count;
		int m = c.Length;
		Vec dc = new Vec(m);
		Vec copy = c.Copy();
		for (int k = 0; k < m; k++) {
			dc[k] = Max(Abs(c[k]), 1) * Pow(2, -26);
		}
		Matrix J = new Matrix(n, m);
		Vec rc = new Vec(n);
		for (int i = 0; i < n; i++) {
			Func<Vec, double> ri = r[i];
			rc[i] = ri(c);
			for (int k = 0; k < m; k++) {
				copy[k] += dc[k];
				J[i, k] = (ri(copy) - rc[i]) / dc[k];
				copy[k] -= dc[k];
			}
		}
		return (J, rc);
	}

	public static double Evaluate(List<Func<Vec, double>> r, Vec c) {
		// Evaluates the sum of squares of a list of functions at a point c
		double sumsquared = 0;
		for (int i = 0; i < r.Count; i++) {
			sumsquared += Pow(r[i](c), 2);
		}
		return sumsquared;
	}

	public static (Vec, Matrix, int) GNewton(List<Func<Vec, double>> r, Vec start, double acc = 1e-4, int maxfev = 1000, double Lmin = 1.0 / 32.0, double Lmax = 1.0) {
		// Applies the Gauss-Newton algorithm to minimize the sum of squares through the vector c
		int iterations = 0;
		Vec c = start.Copy();
		Matrix ccov = null;
		double phi = Evaluate(r, start);
		while (iterations < maxfev) {
			(Matrix J, Vec rc) = Jacobian(r, c);
			Matrix JT = J.Transpose();
			if ((2 * J.Transpose() * rc).Norm() < acc) {
				ccov = QRGS.Inverse(JT * J);
				break;
			}
			Vec Dc = - (QRGS.Inverse(JT * J) * (JT * rc));
			double L = Lmax;
			double phiLc = Evaluate(r, c + L * Dc);
			while (phiLc > phi && L > Lmin) {
				L /= 2;
				phiLc = Evaluate(r, c + L * Dc);
			}
			phi = phiLc;
			c += L * Dc;
			iterations++;
		}
		return (c, ccov, iterations);
	}
	
	public static List<Func<Vec, double>> RosenbrockFunction() {
		Func<Vec, double> r1 = (Vec c) => 1 - c[0];
		Func<Vec, double> r2 = (Vec c) => 10 * (c[1] - c[0] * c[0]);
		return new List<Func<Vec, double>> { r1, r2 };
	}

	public static List<Func<Vec, double>> HimmelblauFunction() {
		Func<Vec, double> r1 = (Vec c) => c[0] * c[0] + c[1] - 11;
		Func<Vec, double> r2 = (Vec c) => c[0] + c[1] * c[1] - 7;
		return new List<Func<Vec, double>> { r1, r2 };
	}

	public static void TestRosenBlau() {
		// Rosenbrock test
		List<Func<Vec, double>> rosenbrock = RosenbrockFunction();

		Vec start = new Vec(2) { [0] = -0.5, [1] = 2.0 };
		(Vec solution, _, int iters) = GNewton(rosenbrock, start);

		Console.WriteLine("\nRosenbrock's valley function test");
		Console.WriteLine($"Starting point: {start}");
		Console.WriteLine("True solution: [1, 1]");
		Console.WriteLine($"Solution: {solution}");
		Console.WriteLine($"Iterations: {iters}");

		// Himmelblau test
		List<Func<Vec, double>> himmelblau = HimmelblauFunction();

		Vec startH = new Vec(2) { [0] = 1, [1] = 5 };
		(solution, _, iters) = GNewton(himmelblau, startH);

		Console.WriteLine("\nHimmelblau's function test");
		Console.WriteLine($"Starting point: {startH}");
		Console.WriteLine("True (nearest) solution: [3, 2]");
		Console.WriteLine($"Solution: {solution}");
		Console.WriteLine($"Iterations: {iters}");
	}

	public static double Gaussian2D(double x, double y, Vec c) {
		double A = c[0]; // Amplitude
		double x0 = c[1]; // Center x
		double y0 = c[2]; // Center y
		double sigmax = c[3]; // Spread in x
		double sigmay = c[4]; // Spread in y
		double result = A * Exp(- (Pow(x - x0, 2) / (2 * Pow(sigmax, 2)) + Pow(y - y0, 2) / (2 * Pow(sigmay, 2))));
		return result;
	}

	static Func<Vec, double> ChiGauss(double x, double y, double dataPoint) {
		return (Vec c) => {
			if (c.Length < 5) {
				Console.WriteLine($"Vector c out of bounds: Length = {c.Length}");
				throw new IndexOutOfRangeException($"Vector c out of bounds: Length = {c.Length}");
			}
			return (Gaussian2D(x, y, c) - dataPoint) / Sqrt(dataPoint);
		};
	}

	static List<Func<Vec, double>> Gaussian2DchiFunctions(Matrix data) {
		List<Func<Vec, double>> chi = new List<Func<Vec, double>>();
		int n = data.Rows;
		int m = data.Cols;

		for (int i = 0; i < n; i++) {
			double x = (double)i;
			for (int j = 0; j < m; j++) {
				double y = (double)j;

				Func<Vec, double> func = ChiGauss(x, y, data[i, j]);
				chi.Add(func);
			}
		}
		return chi;
	}

	public static void FitStar(string data_file) {
		Matrix data = Matrix.Loadtxt(data_file);
		int n = data.Rows;
		int m = data.Cols;
		Vec c0 = new Vec(5);
		c0[0] = 10000;
		c0[1] = (double)n/2;
		c0[2] = (double)m/2;
		c0[3] = (double)n/4;
		c0[4] = (double)m/4;

		var func = Gaussian2DchiFunctions(data);
		(Vec fit, Matrix cov, int iters) = GNewton(func, c0);

		Console.WriteLine($"\nFitting observation of the star LAWD 37 (white dwarf) to 2D Gaussian function");
		Console.WriteLine($"Fitting section is {n} by {m} = {n*m} pixels");
		Console.WriteLine($"Initial guess for parameters: {c0}");
		Console.WriteLine($"Fitted star in {iters} steps");
		Console.WriteLine($"Amplitude: {fit[0]:F1} +/- {Sqrt(cov[0, 0]):F1}");
		Console.WriteLine($"Centroid x: {fit[1]:F4} +/- {Sqrt(cov[1, 1]):F4}");
		Console.WriteLine($"Centroid y: {fit[2]:F4} +/- {Sqrt(cov[2, 2]):F4}");
		Console.WriteLine($"Spread in x: {fit[3]:F4} +/- {Sqrt(cov[3, 3]):F4}");
		Console.WriteLine($"Spread in y: {fit[4]:F4} +/- {Sqrt(cov[4, 4]):F4}");

		using (StreamWriter writer = new StreamWriter("Out.LAWD_37.dat")) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					writer.WriteLine($"{i} {j} {data[i, j]}");
				}
			}
		}

		using (StreamWriter writer = new StreamWriter("Out.LAWD_37_fit.dat")) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					double fittedValue = Gaussian2D(i, j, fit);
					writer.WriteLine($"{i} {j} {fittedValue}");
				}
			}
		}
	}

	public static List<Func<Vec, double>> CreateChiFunction(Func<double, Vec, double> f, Vec x, Vec y, Vec dy = null) {
		List<Func<Vec, double>> chi = new List<Func<Vec, double>>();
		int n = x.Length;

		if (dy == null) {
			dy = new Vec(n);
			for (int i = 0; i < n; i++) {
				dy[i] = 1.0;
			}
		}

		for (int i = 0; i < n; i++) {
			int index = i;  // Capture the current index for the lambda expression
			Func<Vec, double> func = (Vec c) => (f(x[index], c) - y[index]) / dy[index];
			chi.Add(func);
		}
		return chi;
	}

	public static (Vec, Matrix) CurveFit(Func<double, Vec, double> f, Vec c0, Vec x, Vec y, Vec dy = null) {
		var funcs = CreateChiFunction(f, x, y, dy);
		(Vec c, Matrix ccov, _) = GNewton(funcs, c0);
		return (c, ccov);
	}

	public static double wave(double x, Vec c) {
		return c[0] * Sin(c[1] * x + c[3]) * Cos(c[2] * x + c[3]);
	}

	private static Random random = new Random();
	public static double NextGaussian(double std = 1, double mean = 0) {
		double u1 = random.NextDouble();
		double u2 = random.NextDouble();

		// Apply Box-Muller transform
		double randStdNormal = Sqrt(-2.0 * Log(u1)) * Sin(2.0 * PI * u2);
		return randStdNormal * std + mean;
	}

	public static void TestFitting(int n, bool write = false) {
		Vec c = new Vec(4);
		c[0] = 1;
		c[1] = 0.8;
		c[2] = 2;
		c[3] = 0.5;

		Vec c0 = new Vec(4);
		c0[0] = 1.1;
		c0[1] = 1.0;
		c0[2] = 2.5;
		c0[3] = 0.4;

		double a = -3.0;
		double b = 3.0;
		double std = 0.2;

		Vec x = new Vec(n);
		Vec y = new Vec(n);
		for (int i = 0; i < n; i++) {
			x[i] = (double)i / n * (b - a) + a;
			y[i] = wave(x[i], c);
			y[i] += NextGaussian(std: std);
		}

		var (solution, _) = CurveFit(wave, c0, x, y);

		if (write) {
			Console.WriteLine($"\nFitting f(x) = c[0] * Sin(c[1] * x + c[3]) * Cos(c[2] * x + c[3]) to {n} points in the range [{a}, {b}] with standard deviation {std}");
			Console.WriteLine("Original parameters:");
			Console.WriteLine($"c[0] = {c[0]}, c[1] = {c[1]}, c[2] = {c[2]}, c[3] = {c[3]}");
			Console.WriteLine("Fitted parameters:");
			Console.WriteLine($"c[0] = {solution[0]}, c[1] = {solution[1]}, c[2] = {solution[2]}, c[3] = {solution[3]}");

			using (StreamWriter sw = File.CreateText("Out.Generated_points.data")) {
				for (int i = 0; i < n; i++) {
					sw.WriteLine($"{x[i]} {y[i]}");
				}
			}
			using (StreamWriter sw = File.CreateText("Out.fitparameters.txt")) {
				sw.WriteLine($"fit0 = {solution[0]}");
				sw.WriteLine($"fit1 = {solution[1]}");
				sw.WriteLine($"fit2 = {solution[2]}");
				sw.WriteLine($"fit3 = {solution[3]}");
			}
		}
	}

	static void Main(string[] args) {
		int N = 0;
		foreach (var arg in args) {
			if (arg.StartsWith("-N:")) {
				N = Convert.ToInt32(arg.Substring("-N:".Length));
			}
		}
		if (N == 0) {
			TestRosenBlau();
			FitStar("LAWD_37.data");
			TestFitting(100, true);
		}
		else {
			TestFitting(N, false);
		}
	}
}
