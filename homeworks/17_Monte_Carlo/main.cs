using System;
using static System.Math;
using System.IO;
using System.Collections.Generic;

static class MonteCarlo {

static (double, double) PlainMC(Func<Vec, double> f, Vec a, Vec b, int N) {
	int dim = a.Length;
	double V = 1;
	for (int i = 0; i < dim; i++) {
		V *= b[i] - a[i];
	}

	double sum = 0, sum2 = 0;
	var x = new Vec(dim);
	var rnd = new Random();
	for (int i = 0; i < N; i++) {
		for (int k = 0; k < dim; k++) {
			x[k] = a[k] + rnd.NextDouble() * (b[k] - a[k]);
		}
		double fx = f(x);
		sum += fx;
		sum2 += fx * fx;
	}
	double mean = sum / N, sigma = Sqrt(sum2 / N - mean * mean);
	var result = (mean * V, sigma * V / Sqrt(N));
	return result;
}

static double Corput(int n, int b) {
	double q = 0, bk = 1.0 / b;
	while (n > 0) {
		q += (n % b) * bk;
		n /= b;
		bk /= b;
	}
	return q;
}

static Vec Halton(int n, int d, int baseShift = 0) {
	int[] baseArr = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61};
	int maxd = baseArr.Length - baseShift;
	if (d > maxd) throw new ArgumentException($"d should be <= {maxd}");
	Vec x = new Vec(d);
	for (int i = 0; i < d; i++) {
		x[i] = Corput(n + 1, baseArr[i] + baseShift);
	}
	return x;
}

static (double, double) QuasiMC(Func<Vec, double> f, Vec a, Vec b, int N) {
	int dim = a.Length;
	double V = 1;
	for (int i = 0; i < dim; i++) {
		V *= b[i] - a[i];
	}

	double sum1 = 0, sum2 = 0;
	Vec xVec = new Vec(dim);
	for (int i = 0; i < N / 2; i++) {
		Vec x1 = Halton(i, dim);
		Vec x2 = Halton(i + N / 2, dim);

		for (int k = 0; k < dim; k++) {
			xVec[k] = a[k] + x1[k] * (b[k] - a[k]);
		}
		double fx1 = f(xVec);
		sum1 += fx1;

		for (int k = 0; k < dim; k++) {
			xVec[k] = a[k] + x2[k] * (b[k] - a[k]);
		}
		double fx2 = f(xVec);
		sum2 += fx2;
	}

	double mean1 = sum1 / (N / 2);
	double mean2 = sum2 / (N / 2);
	double mean = (mean1 + mean2) / 2;

	double error = Abs(mean1 - mean2) * V / 2;

	var result = (mean * V, error);
	return result;
}

static (double, double) StratMC(Func<Vec, double> f, Vec a, Vec b, int N, int nmin) {
	// Inspired by the solution by Asbjørn Frost Teilmann
	if (N < nmin + 4) return PlainMC(f, a, b, N);

	int dim = a.Length;
	double V = 1;
	for (int i = 0; i < dim; i++) {
		V *= b[i] - a[i];
	}

	var x = new Vec(dim);
	var rnd = new Random();
	double sum = 0;
	double sum2 = 0;
	Matrix sums = new Matrix(2, dim);
	Matrix sums2 = new Matrix(2, dim);

	// Fills nmin points in this volume and calculates two sums, the second one for error estimation.
	// Divides the volume in two along all dimensions, checks which sub-volume each point is in.
	// Calculates the two sums for all dim*2 sub-volumes for later.
	for (int i = 0; i < nmin; i++) {
		for (int k = 0; k < dim; k++) {
			x[k] = a[k] + rnd.NextDouble() * (b[k] - a[k]);
		}
		double fx = f(x);
		sum += fx;
		sum2 += fx * fx;

		for (int k = 0; k < dim; k++) {
			if (x[k] < (a[k] + b[k]) / 2) {
				sums[0, k] += fx;
				sums2[0, k] += fx * fx;
			} else {
				sums[1, k] += fx;
				sums2[1, k] += fx * fx;
			}
		}
	}

	double mean = sum / nmin;
	double sigma = V * Sqrt(sum2 / nmin - mean * mean) / Sqrt(nmin);

	// Find the sub-volume with the highest internal variance.
	int maxDim = 0;
	double maxVar = 0;
	Matrix vars = new Matrix(2, dim);
	for (int k = 0; k < dim; k++) {
		vars[0, k] = (sums2[0, k] - sums[0, k] * sums[0, k] / nmin) / nmin;
		vars[1, k] = (sums2[1, k] - sums[1, k] * sums[1, k] / nmin) / nmin;
		if (vars[0, k] > maxVar) {
			maxVar = vars[0, k];
			maxDim = k;
		}
		if (vars[1, k] > maxVar) {
			maxVar = vars[1, k];
			maxDim = k;
		}
	}

	// Split the volume in two in the dimension of the subvolume with the highest variance.
	Vec a_new = a.Copy();
	Vec b_new = b.Copy();
	a_new[maxDim] = (a[maxDim] + b[maxDim]) / 2;
	b_new[maxDim] = (a[maxDim] + b[maxDim]) / 2;

	// Compare the variance of the two and distribute the remaining N-nmin points between them based on the relative variance.
	// There should, at minimum, be two points in each subvolume.
	int N_est = (int)Floor((N - nmin) / (1 + vars[0, maxDim] / vars[1, maxDim]));
	int N_up = Min(N - nmin - 2, Max(2, N_est));
	int N_down = N - nmin - N_up;

	// Call StratMC recursively. If the number of assigned points is fewer than nmin, PlainMC is returned.
	(double int_down, double sigma_down) = StratMC(f, a, b_new, N_down, nmin);
	(double int_up, double sigma_up) = StratMC(f, a_new, b, N_up, nmin);

	double grandint = ((int_down + int_up) * (N - nmin) + V * mean * nmin) / N;
	double grandsigma = Sqrt(Pow(sigma_down * (N - nmin) / N, 2) + Pow(sigma_up * (N - nmin) / N, 2) + Pow(sigma * nmin / N, 2));
	return (grandint, grandsigma);
}

static double CoolFunc(Vec x) {
	double result = (x[0] + x[1]) * Log(x[0] - x[0] * x[1] + x[1]);
	return result;
}

static Vec coolLimit0 = new Vec(2) { [0] = 0, [1] = 0 };
static Vec coolLimit1 = new Vec(2) { [0] = 1, [1] = 1 };
static double coolAnalRes = PI * PI / 3.0 - 7.0 / 2;

public static void DoTestMC(int N) {
	(double res, double err) = PlainMC(CoolFunc, coolLimit0, coolLimit1, N);

	Console.WriteLine($"\n∫_0^1∫_0^1 (x+y) ln(x - xy + y) dxdy = π²/3 + 7/2 ≈ {coolAnalRes}");
	Console.WriteLine($"Plain Monte Carlo (N = {N}):");
	Console.WriteLine($"Result = {res}");
	Console.WriteLine($"Error = {err}");
	Console.WriteLine($"|Analytical - Numerical| = {Abs(res - coolAnalRes)}\n");

	(res, err) = QuasiMC(CoolFunc, coolLimit0, coolLimit1, N);
	Console.WriteLine($"Quasi Monte Carlo (N = {N}):");
	Console.WriteLine($"Result = {res}");
	Console.WriteLine($"Error = {err}");
	Console.WriteLine($"|Analytical - Numerical| = {Abs(res - coolAnalRes)}\n");
}

public static void DoErrorTest(int N, string output_file) {
	(double plainRes, double plainErr) = PlainMC(CoolFunc, coolLimit0, coolLimit1, N);
	(double quasiRes, double quasiErr) = QuasiMC(CoolFunc, coolLimit0, coolLimit1, N);

	using (StreamWriter sw = File.AppendText(output_file)) {
		sw.WriteLine($"{N}\t{plainErr}\t{Abs(plainRes - coolAnalRes)}\t{quasiErr}\t{Abs(quasiRes - coolAnalRes)}");
	}
}

static double Sphere(Vec x) {
	double r2 = 0;
	for (int i = 0; i < x.Length; i++) {
		r2 += x[i] * x[i];
	}
	if (r2 <= 1) {
		return 1.0;
	} else {
		return 0.0;
	}
}

public static void DoVolumeTest(int N, int nmin) {
	int dim = 3;
	Vec a = new Vec(dim);
	Vec b = new Vec(dim);
	for (int i = 0; i < dim; i++) {
		a[i] = -1;
		b[i] = 1;
	}

	double analyticalVolume = 4.0 / 3.0 * PI;

	Console.WriteLine($"\nCalculating the volume of a unit {dim}-sphere = 4/3 π ≈ {analyticalVolume}");

	(double plainRes, double plainErr) = PlainMC(Sphere, a, b, N);
	Console.WriteLine($"Plain Monte Carlo (N = {N}):");
	Console.WriteLine($"Result = {plainRes}");
	Console.WriteLine($"Error = {plainErr}");
	Console.WriteLine($"|Analytical - Numerical| = {Abs(plainRes - analyticalVolume)}\n");

	(double quasiRes, double quasiErr) = QuasiMC(Sphere, a, b, N);
	Console.WriteLine($"Quasi Monte Carlo (N = {N}):");
	Console.WriteLine($"Result = {quasiRes}");
	Console.WriteLine($"Error = {quasiErr}");
	Console.WriteLine($"|Analytical - Numerical| = {Abs(quasiRes - analyticalVolume)}\n");

	(double stratRes, double stratErr) = StratMC(Sphere, a, b, N, nmin);
	Console.WriteLine($"Stratified Monte Carlo (N = {N}, nmin = {nmin}):");
	Console.WriteLine($"Result = {stratRes}");
	Console.WriteLine($"Error = {stratErr}");
	Console.WriteLine($"|Analytical - Numerical| = {Abs(stratRes - analyticalVolume)}\n");
}

static void Main(string[] args) {
	if (args.Length > 0) {
		int N = int.Parse(args[0]);
		DoErrorTest(N, "Out.ErrorPlot.txt");
	} else {
		DoTestMC(500);
		DoVolumeTest(10000, 200);
	}
}

}
