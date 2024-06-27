using System;
using System.IO;

namespace splines {
public static class Program {
	public static void Main() {
		DoLinterp(9, 0.1, "Out.Linterp.xy.txt", "Out.Linterp.interp.txt");
		DoQspline(9, 0.1, "Out.Qspline.xy.txt", "Out.Qspline.spline.txt");
		DoCspline(0.05, "Cspline.xy.txt", "Out.Cspline.spline.txt");
	}

	static void DoLinterp(int n, double dz, string xy_file, string interp_file) {
		Vec x = new Vec(n);
		Vec y = new Vec(n);
		using (StreamWriter sw = File.CreateText(xy_file)) {
			for (int i = 0; i < n; i++) {
				x[i] = (double)i;
				y[i] = Math.Cos(x[i]);
				sw.WriteLine($"{x[i]} {y[i]}");
			}
		}

		using (StreamWriter sw = File.CreateText(interp_file)) {
			for (double z = x[0]; z < x[n-1]; z += dz) {
				double yz = Linterp(x, y, z);
				sw.WriteLine($"{z} {yz}");
			}
		}
		double pi = 3.14159265;
		double integral1 = LinterpInteg(x, y, pi/2);
		double integral2 = LinterpInteg(x, y, pi);

		Console.WriteLine("\nLinear interpolation");
		Console.WriteLine("Integral of cos(x) from 0 to pi/2 (should be 1):");
		Console.WriteLine($"{integral1}");
		Console.WriteLine("Integral of cos(x) from 0 to pi (should be 0):");
		Console.WriteLine($"{integral2}");
	}

	static void DoQspline(int n, double dz, string xy_file, string spline_file) {
			Vec x = new Vec(n);
			Vec y = new Vec(n);
			using (StreamWriter sw = File.CreateText(xy_file)) {
				for (int i = 0; i < n; i++) {
					x[i] = (double)i;
					y[i] = Math.Sin(x[i]);
					sw.WriteLine($"{x[i]} {y[i]}");
				}
			}

			Qspline spline = new Qspline(x, y);

			using (StreamWriter sw = File.CreateText(spline_file)) {
				for (double z = x[0]; z < x[n - 1]; z += dz) {
					double yz = spline.Evaluate(z);
					double dydx = spline.Derivative(z);
					sw.WriteLine($"{z} {yz} {dydx}");
				}
			}
		double pi = 3.14159265;
		double integral1 = spline.Integrate(pi);
		double integral2 = spline.Integrate(2 * pi);

		Console.WriteLine("\nQuadratic spline");
		Console.WriteLine("Integral of sin(x) from 0 to pi (should be 2):");
		Console.WriteLine($"{integral1}");
		Console.WriteLine("Integral of sin(x) from 0 to 2*pi (should be 0):");
		Console.WriteLine($"{integral2}");
	}

	static void DoCspline(double dz, string xy_file, string spline_file) {
		Vec x = new Vec(0);
		Vec y = new Vec(0);

		using (StreamReader sr = File.OpenText(xy_file)) {
			string line;
			while ((line = sr.ReadLine()) != null) {
				string[] parts = line.Split(new char[] { ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries);
				if (parts.Length == 2) {
					double xVal = double.Parse(parts[0]);
					double yVal = double.Parse(parts[1]);
					x = x.Append(xVal);
					y = y.Append(yVal);
				}
			}
		}

		Cspline spline = new Cspline(x, y);

		using (StreamWriter sw = File.CreateText(spline_file)) {
			for (double z = x[0]; z < x[x.Length - 1]; z += dz) {
				double yz = spline.Evaluate(z);
				//double dydx = spline.Derivative(z);
				sw.WriteLine($"{z} {yz}");
			}
		}
	}

	public static int Binsearch(Vec x, double z) {
	if( z<x[0] || z>x[x.Length-1] ) throw new Exception("binsearch: bad z");
	int i=0, j=x.Length-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
	}

	public static double Linterp(Vec x, Vec y, double z){
        int i=Binsearch(x,z);
        double dx=x[i+1]-x[i]; if(!(dx>0)) throw new Exception("uups...");
        double dy=y[i+1]-y[i];
        return y[i]+dy/dx*(z-x[i]);
        }

	public static double LinterpInteg(Vec x, Vec y) {
		double integral = 0;
		for (int i = 0; i < x.Length - 1; i++) {
			integral += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i]);
		}
		return integral;
	}

	public static double LinterpInteg(Vec x, Vec y, double z) {
		int n = Binsearch(x, z) + 2;
		Vec xs = new Vec(n);
		Vec ys = new Vec(n);
		for (int i = 0; i < n - 1; i++) {
			xs[i] = x[i];
			ys[i] = y[i];
		}
		xs[n-1] = z;
		ys[n-1] = Linterp(x, y, z);
		
		double integral = LinterpInteg(xs, ys);
		return integral;
	}
}

public class Qspline{
	int n;
	Vec x, y, b, c;
	public Qspline(Vec xs, Vec ys) {
		n = xs.Length;
		x = xs.Copy();
		y = ys.Copy();
		b = new Vec(n - 1);
		c = new Vec(n - 1);

		Vec h = new Vec(n - 1);
		Vec p = new Vec(n - 1);
		
		for (int i = 0; i < n - 1; i++) {
			h[i] = x[i + 1] - x[i];
			p[i] = (y[i + 1] - y[i]) / h[i];
		}

		c[0] = 0;

		// Recursion up
		for (int i = 0; i < n - 2; i++) {
		    c[i + 1] = (p[i + 1] - p[i] - c[i] * h[i]) / h[i + 1];
		}

		c[n - 2] /= 2;

		// Recursion down
		for (int i = n - 3; i >= 0; i--) {
		    c[i] = (p[i + 1] - p[i] - c[i + 1] * h[i + 1]) / h[i];
		}

		for (int i = 0; i < n - 1; i++) {
		    b[i] = p[i] - c[i] * h[i];
		}
	}

	public double Evaluate(double z) {
		int i = Program.Binsearch(x, z);
		double h = z - x[i];
		return y[i] + h * (b[i] + h * c[i]);
	}

	public double Derivative(double z) {
		int i = Program.Binsearch(x, z);
		return 2 * c[i] * (z - x[i]) + b[i];
	}

	private double ydx(double xi, double xi1, double yi, double bi, double ci) {
		double result = -(1.0/6) * (xi - xi1) * (2 * xi*xi * ci - 3 * xi * bi - 4 * xi * xi1 * ci + 3 * xi1 * bi + 2 * xi1*xi1 * ci + 6 * yi);
		return result;
	}
	
	public double Integrate(double z) {
		int n = Program.Binsearch(x, z);
		double integral = 0;
		for (int i = 0; i < n; i++) {
			integral += ydx(x[i], x[i+1], y[i], b[i], c[i]);
		}
		integral += ydx(x[n], z, y[n], b[n], c[n]);
		return integral;
	}
}

public class Cspline {
	private int n;
	private Vec x, y, b, c, d;

	public Cspline(Vec xs, Vec ys) {
		n = xs.Length;
		x = xs.Copy();
		y = ys.Copy();
		b = new Vec(n);
		c = new Vec(n - 1);
		d = new Vec(n - 1);

		Vec h = new Vec(n - 1);
		Vec p = new Vec(n - 1);

		for (int i = 0; i < n - 1; i++) {
			h[i] = x[i + 1] - x[i];
			if (h[i] <= 0) throw new ArgumentException("x values must be in ascending order.");
			p[i] = (y[i + 1] - y[i]) / h[i];
		}

		Vec D = new Vec(n);
		Vec Q = new Vec(n - 1);
		Vec B = new Vec(n);

		D[0] = 2;
		for (int i = 0; i < n - 2; i++) {
			D[i + 1] = 2 * h[i] / h[i + 1] + 2;
		}
		D[n - 1] = 2;

		Q[0] = 1;
		for (int i = 0; i < n - 2; i++) {
			Q[i + 1] = h[i] / h[i + 1];
		}

		for (int i = 0; i < n - 2; i++) {
			B[i + 1] = 3 * (p[i] + p[i + 1] * h[i] / h[i + 1]);
		}
		B[0] = 3 * p[0];
		B[n - 1] = 3 * p[n - 2];

		for (int i = 1; i < n; i++) {
			D[i] -= Q[i - 1] / D[i - 1];
			B[i] -= B[i - 1] / D[i - 1];
		}

		b[n - 1] = B[n - 1] / D[n - 1];

		for (int i = n - 2; i >= 0; i--) {
			b[i] = (B[i] - Q[i] * b[i + 1]) / D[i];
		}

		for (int i = 0; i < n - 1; i++) {
			c[i] = (-2 * b[i] - b[i + 1] + 3 * p[i]) / h[i];
			d[i] = (b[i] + b[i + 1] - 2 * p[i]) / (h[i] * h[i]);
		}
	}

	public double Evaluate(double z) {
		if (z < x[0] || z > x[n - 1]) throw new ArgumentException("z is out of range.");
		int i = Program.Binsearch(x, z);
		double h = z - x[i];
		return y[i] + h * (b[i] + h * (c[i] + h * d[i]));
	}

	public double Derivative(double z) {
		int i = Program.Binsearch(x, z);
		return 2 * c[i] * (z - x[i]) + 3 * d[i] * (z - x[i])*(z - x[i]) + b[i];
	}

	private double ydx(double xi, double z, double yi, double bi, double ci, double di) {
		double result = (Math.Pow(xi, 2) * bi) / 2 
			- xi * yi 
			- (ci * Math.Pow(xi - z, 3)) / 3 
			+ (di * Math.Pow(xi - z, 4)) / 4 
			- xi * bi * z 
			+ yi * z 
			+ (bi * Math.Pow(z, 2)) / 2;
		return result;
	}
	
	public double Integrate(double z) {
		int n = Program.Binsearch(x, z);
		double integral = 0;
		for (int i = 0; i < n; i++) {
			integral += ydx(x[i], x[i+1], y[i], b[i], c[i], d[i]);
		}
		integral += ydx(x[n], z, y[n], b[n], c[n], d[n]);
		return integral;
	}
}
}
