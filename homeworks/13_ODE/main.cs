using System;
using System.IO;
using System.Collections.Generic;
using splines;

public class ODE {
	public static (Vec, Vec) RKstep23(
			Func<double, Vec, Vec> f,
			double x,
			Vec y,
			double h
			) {
		Vec k0 = f(x, y);
		Vec k1 = f(x + 1.0/2 * h, y + 1.0/2 * k0 * h);
		Vec k2 = f(x + 3.0/4 * h, y + 3.0/4 * k1 * h);
		Vec yh = y + 2.0/9 * k0 * h + 3.0/9 * k1 * h + 4.0/9 * k2 * h;
		Vec dy = yh - (y + k1 * h);

		return (yh, dy);
	}

	public static (List<double>, List<Vec>) Driver(
			Func<double, Vec, Vec> f,
			(double, double) interval,
			Vec ystart,
			double h = 0.125,
			double acc = 0.01,
			double eps = 0.01
			) {
		var (a, b) = interval;
		double x = a;
		Vec y = ystart.Copy();
		var xlist = new List<double> { x };
		var ylist = new List<Vec> { y };

		while (true) {
			if (x >= b) return (xlist, ylist);
			if (x + h > b) h = b - x;

			var (yh, dy) = RKstep23(f, x, y, h);
			double tol = (acc + eps * yh.Norm()) * Math.Sqrt(h / (b - a));
			double err = dy.Norm();

			if (err <= tol) {
				x += h;
				y = yh;
				xlist.Add(x);
				ylist.Add(y);
			}

			h *= Math.Min(Math.Pow(tol / err, 0.25) * 0.95, 2);
		}
	}

	public static (List<double>, List<Vec>) CInterpolator (List<double> xlist, List<Vec> ylist, double dx) {
		Vec x = new Vec(xlist.Count);
		for (int i = 0; i < xlist.Count; i++) x[i] = xlist[i];
		List<Vec> y = TransposeList(ylist);

		int N = (int)((x[x.Length-1] - x[0])/dx);
		Vec splinex = new Vec(N);
		List<Vec> spliney = new List<Vec>();
		for (int i = 0; i < N; i++) splinex[i] = x[0] + i * dx;
		for (int i = 0; i < y.Count; i++) {
			Vec newspliney = new Vec(N);
			Cspline spline = new Cspline(x, y[i]);
			for (int j = 0; j < N; j++) {
				newspliney[j] = spline.Evaluate(splinex[j]);
			}
			spliney.Add(newspliney);
		}
		List<double> splinexlist = new List<double>();
		for (int i = 0; i < N; i++) splinexlist.Add(splinex[i]);
		List<Vec> splineylist = TransposeList(spliney);
		return (splinexlist, splineylist);
	}

	public static (List<double>, List<Vec>) AlternativeDriver(
			Func<double, Vec, Vec> f,
			(double, double) interval,
			Vec ystart,
			double dx,
			double h = 0.125,
			double acc = 0.01,
			double eps = 0.01
			) {
		var (xlist, ylist) = Driver(f, interval, ystart, h, acc, eps);
		var (splinexlist, splineylist) = CInterpolator(xlist, ylist, dx);
		return (splinexlist, splineylist);
	}

	public static List<Vec> TransposeList(List<Vec> input) {
		int n = input.Count;
		int m = input[0].Length;
		List<Vec> result = new List<Vec>();

		for (int i = 0; i < m; i++) {
			Vec newVector = new Vec(n);
			for (int j = 0; j < n; j++) {
				newVector[j] = input[j][i];
			}
			result.Add(newVector);
		}

		return result;
	}

	public static void WriteToFile(List<double> xlist, List<Vec> ylist, string filename) {
		using (StreamWriter sw = File.CreateText(filename)) {
			for (int i = 0; i < xlist.Count; i++) {
				string line = $"{xlist[i]}";
				for (int j = 0; j < ylist[0].Length; j++) {
					line += $"\t{ylist[i][j]}";
				}
				sw.WriteLine(line);
			}
		}
	}

	public static Vec HarmOscillator(double x, Vec y) {
		// ODE implementation of u'' = -u or a = -x, where y[0] is position and y[1] is velocity
		Vec f = new Vec(2);
		f[0] = y[1];
		f[1] = -y[0];
		return f;
	}

	static void DoHarmOscillator(string filename) {
		Vec y0 = new Vec(2);
		y0[0] = 1.0;
		y0[1] = 0.0;
		var (xlist, ylist) = Driver(HarmOscillator, (0.0, 10.0), y0);
		WriteToFile(xlist, ylist, filename);
	}

	public static Vec FricOscillator(double x, Vec y) {
		// ODE implementation of theta' = omega and omega' = -b*omega - c*sin(theta)
		double b = 0.25, c = 5.0;
		Vec f = new Vec(2);
		f[0] = y[1];
		f[1] = -b * y[1] - c * Math.Sin(y[0]);
		return f;
	}

	static void DoFricOscillator(string filename) {
		Vec y0 = new Vec(2);
		y0[0] = Math.PI - 0.1; // Initial angle
		y0[1] = 0.0;		   // Initial angular velocity
		var (xlist, ylist) = Driver(FricOscillator, (0.0, 10.0), y0);
		WriteToFile(xlist, ylist, filename);
	}

	public static Vec EquatorialMotion(double phi, Vec u, double epsilon) {
		Vec f = new Vec(2);
		f[0] = u[1];
		f[1] = 1 - u[0] + epsilon * u[0] * u[0];
		return f;
	}

	static void DoEquatorialMotion(string filename, double epsilon, double initialU1) {
		Vec y0 = new Vec(2);
		y0[0] = 1.0; // Initial condition u(0)
		y0[1] = initialU1; // Initial condition u'(0)
		var (xlist, ylist) = AlternativeDriver((phi, u) => EquatorialMotion(phi, u, epsilon), (0.0, 100.0), y0, 0.01, 0.01, 0.01, 0.01);

		WriteToFile(xlist, ylist, filename);
	}

	public static Vec ThreeBody(double t, Vec z) {
		//      0    1    2    3    4    5    6   7   8   9   10  11
		// z = {x'1, y'1, x'2, y'2, x'3, y'3, x1, y1, x2, y2, x3, y3}
		Vec f = new Vec(12);
		double r12C = Math.Pow(Math.Pow(z[6] - z[8], 2) + Math.Pow(z[7] - z[9], 2), 3.0/2);
		double r13C = Math.Pow(Math.Pow(z[6] - z[10], 2) + Math.Pow(z[7] - z[11], 2), 3.0/2);
		double r23C = Math.Pow(Math.Pow(z[8] - z[10], 2) + Math.Pow(z[9] - z[11], 2), 3.0/2);

		f[0] = (z[8] - z[6])/r12C + (z[10] - z[6])/r13C;  // x1'
		f[1] = (z[9] - z[7])/r12C + (z[11] - z[7])/r13C;  // y1'
		f[2] = (z[6] - z[8])/r12C + (z[10] - z[8])/r23C;  // x2'
		f[3] = (z[7] - z[9])/r12C + (z[11] - z[9])/r23C;  // y2'
		f[4] = (z[6] - z[10])/r13C + (z[8] - z[10])/r23C; // x3'
		f[5] = (z[7] - z[11])/r13C + (z[9] - z[11])/r23C; // y3'

		for (int i = 6; i < 12; i++) {
		    f[i] = z[i - 6];
		}

		return f;
	}

	static void DoThreeBody(string filename) {
		Vec z0 = new Vec(12);

		z0[0] = 0.4662036850;   // x1'
		z0[1] = 0.4323657300;   // y1'
		z0[2] = -0.93240737;    // x2'
		z0[3] = -0.86473146;    // y2'
		z0[4] = 0.4662036850;   // x3'
		z0[5] = 0.4323657300;   // y3'

		z0[6] = -0.97000436;  // x1
		z0[7] = 0.24308753;   // y1
		z0[8] = 0;            // x2
		z0[9] = 0;            // y2
		z0[10] = 0.97000436;  // x3
		z0[11] = -0.24308753; // y3

		var (xlist, ylist) = AlternativeDriver(ThreeBody, (0.0, 10.0), z0, 0.1);

		WriteToFile(xlist, ylist, filename);
	}

	static void Main() {
		DoHarmOscillator("Out.HarmOscillator.txt");
		DoFricOscillator("Out.FricOscillator.txt");

		DoEquatorialMotion("Out.EquatorialMotion_Circular.txt", 0.0, 0.0);
		DoEquatorialMotion("Out.EquatorialMotion_Elliptical.txt", 0.0, -0.5);
		DoEquatorialMotion("Out.EquatorialMotion_Precession.txt", 0.01, -0.5);

		DoThreeBody("Out.ThreeBody.txt");
	}
}
