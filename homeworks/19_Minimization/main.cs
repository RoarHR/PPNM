using System;
using System.IO;
using static System.Math;

public static class Minimization {
	public static Vec Gradient(Func<Vec, double> f, Vec x) {
		Vec gf = new Vec(x.Length);
		double fx = f(x);
		for (int i = 0; i < x.Length; i++) {
			double dx = Max(Abs(x[i]), 1) * Pow(2, -26);
			x[i] += dx;
			gf[i] = (f(x) - fx) / dx;
			x[i] -= dx;
		}
		return gf;
	}

	public static Matrix Hessian(Func<Vec, double> f, Vec x) {
		Matrix H = new Matrix(x.Length, x.Length);
		Vec gfx = Gradient(f, x);
		for (int j = 0; j < x.Length; j++) {
			double dx = Max(Abs(x[j]), 1) * Pow(2, -13);
			x[j] += dx;
			Vec dgf = Gradient(f, x) - gfx;
			for (int i = 0; i < x.Length; i++) H[i, j] = dgf[i] / dx;
			x[j] -= dx;
		}
		return H;
	}

	public static (Matrix, Vec) HessianAndGradient(Func<Vec, double> f, Vec x) {
		int n = x.Length;
		Matrix H = new Matrix(n, n);
		Vec gf = new Vec(n);
		Vec dx = new Vec(n);

		// Initialize dx vector
		for (int i = 0; i < n; i++) {
			dx[i] = Max(Abs(x[i]), 1) * Pow(2, -26);
		}

		double fx = f(x);

		for (int j = 0; j < n; j++) {
			for (int k = j; k < n; k++) {
				x[j] += dx[j];
				x[k] += dx[k];
				double f1 = f(x);

				x[j] -= 2 * dx[j];
				x[k] -= 2 * dx[k];
				double f4 = f(x);

				if (j == k) {
					H[j, k] = (f1 + f4 - 2 * fx) / (4 * dx[j] * dx[k]);
					gf[j] = (f1 - f4) / (4 * dx[j]);

					// Reset x[j] and x[k]
					x[j] += dx[j];
					x[k] += dx[k];
				} else {
					x[k] += 2 * dx[k];
					double f2 = f(x);

					x[k] -= 2 * dx[k];
					x[j] += 2 * dx[j];
					double f3 = f(x);

					H[j, k] = (f1 - f2 - f3 + f4) / (4 * dx[j] * dx[k]);
					H[k, j] = H[j, k]; // Ensure symmetry

					// Reset x[j] and x[k]
					x[j] -= dx[j];
					x[k] += dx[k];
				}
			}
		}

		return (H, gf);
	}

	public static (Vec, int) Newton(Func<Vec, double> f, Vec x, double acc = 1e-3, int maxfev = 1000, double Lmin = 1.0 / 16, double Lmax = 1.0, bool CD = false) {
		if (CD) return NewtonCentralDifference(f, x, acc, maxfev, Lmin);
		int iterations = 0;
		double fx = f(x);
		while (iterations < maxfev) {
			Vec gf = Gradient(f, x);
			if (gf.Norm() < acc) break;
			Matrix H = Hessian(f, x);
			(Matrix HQ, Matrix HR) = QRGS.Decomp(H);
			//HQ.Print("HQ:");
			//HR.Print("HR:");
			Vec dx = QRGS.Solve(HQ, HR, -gf);
			double L = Lmax;
			double fLx = f(x + L * dx);
			while (fLx > fx && L > Lmin) {
				L /= 2;
				fLx = f(x + L * dx);
			}
			x.Print("x:");
			gf.Print("∇f");
			dx.Print("dx:");
			Console.WriteLine($"Stepping {L} · dx");
			fx = fLx;
			x += L * dx;
			iterations++;
		}
		return (x, iterations);
	}

	private static (Vec, int) NewtonCentralDifference(Func<Vec, double> f, Vec x, double acc = 1e-3, int maxfev = 1000, double Lmin = 1.0 / 16, double Lmax = 1.0) {
		int iterations = 0;
		double fx = f(x);
		while (iterations < maxfev) {
			(Matrix H, Vec gf) = HessianAndGradient(f, x);
			if (gf.Norm() < acc) break;
			(Matrix HQ, Matrix HR) = QRGS.Decomp(H);
			Vec dx = QRGS.Solve(HQ, HR, -gf);
			double L = Lmax;
			double fLx = f(x + L * dx);
			while (fLx > fx && L > Lmin) {
				L /= 2;
				fLx = f(x + L * dx);
			}
			x.Print("x:");
			gf.Print("∇f");
			dx.Print("dx:");
			Console.WriteLine($"Stepping {L} · dx");
			fx = fLx;
			x += L * dx;
			iterations++;
		}
		return (x, iterations);
	}

	public static void FindRosenbrocksMinimum() {
		Func<Vec, double> rosenbrock = v => Pow(1 - v[0], 2) + 100 * Pow(v[1] - v[0] * v[0], 2);
		Vec start = new Vec(2) { [0] = 1.5, [1] = 1.5 };

		Console.WriteLine($"Testing the forward and central difference methods (FD; CD) at (x, y) = ({start[0]:F2}, {start[1]:F2}) for Rosenbrock's valley function.");

		// Calculate gradient and Hessian using forward difference
		Vec gradFD = Gradient(rosenbrock, start);
		Matrix hessFD = Hessian(rosenbrock, start);
		
		Console.WriteLine("Forward Difference Method:");
		Console.WriteLine("Gradient (FD):");
		gradFD.Print();
		Console.WriteLine("Hessian (FD):");
		hessFD.Print();

		// Calculate gradient and Hessian using central difference
		(Matrix hessCD, Vec gradCD) = HessianAndGradient(rosenbrock, start);

		Console.WriteLine("Central Difference Method:");
		Console.WriteLine("Gradient (CD):");
		gradCD.Print();
		Console.WriteLine("Hessian (CD):");
		hessCD.Print();

		// Perform Newton's minimization using forward difference
		var (minimumFD, stepsFD) = Newton(rosenbrock, start, CD: false);
		Console.WriteLine($"Rosenbrock's function minimum (FD):\t{minimumFD}\tSteps:\t{stepsFD}");

		// Perform Newton's minimization using central difference
		var (minimumCD, stepsCD) = Newton(rosenbrock, start, CD: true);
		Console.WriteLine($"Rosenbrock's function minimum (CD):\t{minimumCD}\tSteps:\t{stepsCD}");
	}

	public static void FindHimmelblausMinimum() {
		Func<Vec, double> himmelblau = v => Pow(v[0] * v[0] + v[1] - 11, 2) + Pow(v[0] + v[1] * v[1] - 7, 2);
		Vec start = new Vec(2) { [0] = 3.0, [1] = 3.0 };

		var (minimumFD, stepsFD) = Newton(himmelblau, start, CD: false);
		Console.WriteLine($"Himmelblau's function minimum (FD):\t{minimumFD}\tSteps:\t{stepsFD}");

		var (minimumCD, stepsCD) = Newton(himmelblau, start, CD: true);
		Console.WriteLine($"Himmelblau's function minimum (CD):\t{minimumCD}\tSteps:\t{stepsCD}");
	}

	public static double BreitWigner(double E, double m, double width, double A) {
		double F = A / (Pow(E - m, 2) + Pow(width, 2) / 4);
		return F;
	}

	public static double BreitWignerDeviation(Vec x, Vec E, Vec signal, Vec dsignal) {
		double D = 0;
		for (int i = 0; i < E.Length; i++) {
			D += Pow((BreitWigner(E[i], x[0], x[1], x[2]) - signal[i]) / dsignal[i], 2);
		}
		return D;
	}

	public static (Vec, int) FitBreitWigner(Vec p0) {
		Matrix data = Matrix.Loadtxt();
		Vec E = data[0];
		Vec signal = data[1];
		Vec dsignal = data[2];
		var (x, steps) = Newton(nam => BreitWignerDeviation(nam, E, signal, dsignal), p0);
		return (x, steps);
	}

	public static void FindHiggsBoson(string out_file) {
		Vec p0 = new Vec(3);
		p0[0] = 125; // m[GeV/c²]
		p0[1] = 1; // Γ
		p0[2] = 1; // A
		var (result, steps) = FitBreitWigner(p0);
		Console.WriteLine($"Fitted the Breit-Wigner function in {steps} steps");
		Console.WriteLine($"m = {result[0]} GeV/c²");
		Console.WriteLine($"Γ = {result[1]}");
		Console.WriteLine($"A = {result[2]}");
		using (StreamWriter sw = File.CreateText(out_file)) {
			for (double E = 101; E < 159; E += 0.1) {
				sw.WriteLine($"{E} {BreitWigner(E, result[0], result[1], result[2])}");
			}
		}
	}

	static void Main() {
		FindRosenbrocksMinimum();
		FindHimmelblausMinimum();
		FindHiggsBoson("Out.HiggsBoson.txt");
	}
}
