using System;
using static System.Math;
using System.IO;
using System.Collections.Generic;

public class RootFinding {
	public static Matrix Jacobian(Func<Vec, Vec> f, Vec x, Vec fx = null, Vec dx = null) {
		if (dx == null) {
			dx = new Vec(x.Length);
			for (int i = 0; i < x.Length; i++) {
				dx[i] = Max(Abs(x[i]), 1) * Pow(2, -26);
			}
		}
		if (fx == null) fx = f(x);
		int n = x.Length;
		Matrix J = new Matrix(n, n);
		for (int j = 0; j < n; j++) {
			x[j] += dx[j];
			Vec df = f(x) - fx;
			for (int i = 0; i < n; i++) J[i, j] = df[i] / dx[j];
			x[j] -= dx[j];
		}
		return J;
	}

	public static Vec Newton(Func<Vec, Vec> f, Vec start, double acc = 1e-2, Vec dx = null) {
		Vec x = start.Copy();
		Vec fx = f(x);
		Vec z, fz;
		while (true) {
			if (fx.Norm() < acc) break;
			Matrix J = Jacobian(f, x, fx, dx);
			var (Q, R) = QRGS.Decomp(J);
			Vec Dx = QRGS.Solve(Q, R, -fx);
			double lambda = 1;
			while (true) {
				z = x + lambda * Dx;
				fz = f(z);
				if (fz.Norm() < (1 - lambda / 2) * fx.Norm()) break;
				if (lambda < Pow(2, -26)) break;
				lambda /= 2;
			}
			x = z;
			fx = fz;
		}
		return x;
	}

	public static double Newton(Func<double, double> f, double start, double acc = 1e-2, double? dx = null) {
		Func<Vec, Vec> fvec = x => {
			Vec fx = new Vec(1);
			fx[0] = f(x[0]);
			return fx;
		};
		Vec startvec = new Vec(1) { [0] = start };
		if (dx == null) return Newton(fvec, startvec, acc)[0];
		else {
			Vec dxvec = new Vec(1) { [0] = (double)dx };
			return Newton(fvec, startvec, acc, dxvec)[0];
		}
	}

	public static Vec NewtonLineSearch(Func<Vec, Vec> f, Vec start, double acc = 1e-2, Vec dx = null) {
		Vec x = start.Copy();
		Vec fx = f(x);
		Vec z, fz;
		while (true) {
			if (fx.Norm() < acc) break;
			Matrix J = Jacobian(f, x, fx, dx);
			var (Q, R) = QRGS.Decomp(J);
			Vec Dx = QRGS.Solve(Q, R, -fx);
			double lambda = 1;
			double phi0 = 0.5 * fx.Norm() * fx.Norm();
			double dphi0 = -phi0;
			
			while (true) {
				z = x + lambda * Dx;
				fz = f(z);
				double phi = 0.5 * fz.Norm() * fz.Norm();
				if (phi < (1 - lambda / 2) * phi0) break;
				if (lambda < Pow(2, -26)) break;
				
				double c = (phi - phi0 - dphi0 * lambda) / (lambda * lambda);
				lambda = -dphi0 / (2 * c);
			}
			
			x = z;
			fx = fz;
		}
		return x;
	}

	public static Vec GradientRosenbrock(Vec v) {
		double x = v[0], y = v[1];
		return new Vec(2) { [0] = 2 * (-1 + x + 200 * x * x * x - 200 * x * y), [1] = 200 * (y - x * x) };
	}

	public static Vec GradientHimmelblau(Vec v) {
		double x = v[0], y = v[1];
		return new Vec(2) { [0] = 2 * (x * (x * x - 3 * y + 11) + y - 11), [1] = 2 * (y * (y * y - 3 * x + 7) + x - 7) };
	}


	public static Vec SSchroedinger(double r, Vec f, double E) {
		Vec df = new Vec(2);
		df[0] = f[1];
		df[1] = -2 * (E + 1 / r) * f[0];
		return df;
	}

	public static (Vec, Vec) SolveSSchroedinger(double E, double rmin, double rmax, double acc = 0.01, double eps = 0.01, bool makespline = false, double dr = 0.01) {
		Vec f0 = new Vec(2);
		f0[0] = rmin - rmin * rmin;
		f0[1] = 1 - 2 * rmin;
		List<double> rlist;
		List<Vec> fslist;
		if (!makespline) (rlist, fslist) = ODE.Driver((r, f) => SSchroedinger(r, f, E), (rmin, rmax), f0, 0.01, acc, eps);
		else (rlist, fslist) = ODE.AlternativeDriver((r, f) => SSchroedinger(r, f, E), (rmin, rmax), f0, dr, 0.01, acc, eps);

		Vec rs = new Vec(rlist.Count), fs = new Vec(rlist.Count);
		for (int i = 0; i < rlist.Count; i++) {
			rs[i] = rlist[i];
			fs[i] = fslist[i][0];
		}
		return (rs, fs);
	}

	public static double ME(double E, double rmin, double rmax, double acc = 0.01, double eps = 0.01) {
		var (rs, fs) = SolveSSchroedinger(E, rmin, rmax, acc, eps);
		double ME = fs[fs.Length - 1];
		return ME;
	}

	public static double SolveForE(double rmin, double rmax, double acc = 0.01, double eps = 0.01) {
		double Estart = -1.0;
		double Esolution = Newton(E => ME(E, rmin, rmax, acc, eps), Estart);
		return Esolution;
	}

	public static void WriteWaveFunctionToFile(string fileName, Vec rs, Vec fs) {
		using (StreamWriter sw = File.CreateText(fileName)) {
			for (int i = 0; i < rs.Length; i++) {
				sw.WriteLine($"{rs[i]} {fs[i]}");
			}
		}
	}

	public static void InvestigateConvergence(string fileName, double rmin, double rmax, double acc, double eps) {
		double E = SolveForE(rmin, rmax, acc, eps);
		using (StreamWriter sw = File.AppendText(fileName)) {
			sw.WriteLine($"{rmin} {rmax} {acc} {eps} {E}");
		}
	}

	static void Main(string[] args) {
		string fileName = null;
		double? rmin = null, rmax = null, acc = null, eps = null;

		foreach (var arg in args) {
			if (arg.StartsWith("-filename:")) {
				fileName = arg.Substring("-filename:".Length);
			} else if (arg.StartsWith("-rmin:")) {
				rmin = Convert.ToDouble(arg.Substring("-rmin:".Length));
			} else if (arg.StartsWith("-rmax:")) {
				rmax = Convert.ToDouble(arg.Substring("-rmax:".Length));
			} else if (arg.StartsWith("-acc:")) {
				acc = Convert.ToDouble(arg.Substring("-acc:".Length));
			} else if (arg.StartsWith("-eps:")) {
				eps = Convert.ToDouble(arg.Substring("-eps:".Length));
			}
		}

		if (fileName != null && rmin.HasValue && rmax.HasValue && acc.HasValue && eps.HasValue) {
			InvestigateConvergence(fileName, rmin.Value, rmax.Value, acc.Value, eps.Value);
		}
	       	else {
			Vec startRosenbrock = new Vec(2) { [0] = 1.5, [1] = 1.5 };
			Vec extremumRosenbrock = Newton(GradientRosenbrock, startRosenbrock);
			Vec extremumRosenbrockLS = NewtonLineSearch(GradientRosenbrock, startRosenbrock);
			Console.WriteLine($"\nExtremum of Rosenbrock's function: {extremumRosenbrock}");
			Console.WriteLine($"Using quadratic interpolation line-search: {extremumRosenbrockLS}");

			Vec startHimmelblau = new Vec(2) { [0] = 1.5, [1] = 1.5 };
			Vec extremumHimmelblau = Newton(GradientHimmelblau, startHimmelblau);
			Vec extremumHimmelblauLS = NewtonLineSearch(GradientHimmelblau, startHimmelblau);
			Console.WriteLine($"\nExtremum of Himmelblau's function: {extremumHimmelblau}");
			Console.WriteLine($"Using quadratic interpolation line-search: {extremumHimmelblauLS}");

			double E = SolveForE(1e-3, 8);
			Console.WriteLine($"\nEnergy level for Hydrogen atom: {E}");

			var (rs, fs) = SolveSSchroedinger(E, 1e-3, 8, makespline: true);

			string fileNameOutput = "Out.WaveFunction.txt";
			WriteWaveFunctionToFile(fileNameOutput, rs, fs);
		}
	}
}
