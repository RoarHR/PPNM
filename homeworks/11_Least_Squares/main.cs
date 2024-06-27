using System;
using static System.Math;

public class Least_Squares {
	public static (Vec, Matrix) Lsfit (Func<double,double>[] fs, Vec x, Vec y, Vec dy) {
		Matrix A = new Matrix(x.Length, fs.Length);
		Vec b = new Vec(x.Length);
		for (int i = 0; i < x.Length; i++) {
			for (int k = 0; k < fs.Length; k++) {
				A[i, k] = fs[k](x[i])/dy[i];
			}
			b[i] = y[i] / dy[i];
		}
		
		(Matrix Q, Matrix R) = QRGS.Decomp(A);
		Vec cs = QRGS.Solve(Q, R, b);
		Matrix I = new Matrix(R.Cols, R.Cols);
		I.Identity();
		Matrix RI = QRGS.Inverse(I, R);
		Matrix Cov = RI * RI.Transpose();
		return (cs, Cov);
	}

	public static void Main() {
		Matrix data = Matrix.Loadtxt();
		Vec x = new Vec(data.Rows);
		Vec y = new Vec(data.Rows);
		Vec dy = new Vec(data.Rows);

		x = data[0];
		y = data[1];
		dy = data[2];

		Func<double, double>[] fs = new Func<double, double>[] {z => 1.0, z => -z};

		for (int i = 0; i < x.Length; i++) {
			dy[i] = dy[i]/y[i];
			y[i] = Math.Log(y[i]);
		}

		(Vec cs, Matrix Cov) = Lsfit(fs, x, y, dy);
		double a = Math.Exp(cs[0]);
		double lambda = cs[1];
		double da = a * Math.Sqrt(Cov[0, 0]);
		double dlambda = Math.Sqrt(Cov[1, 1]);

		Console.Error.WriteLine($"a = {a}");
		Console.Error.WriteLine($"lambda = {lambda}");
		Console.Error.WriteLine($"da = {da}");
		Console.Error.WriteLine($"dlambda = {dlambda}");

		double hl = Math.Log(2) / lambda;
		double dhl = Math.Log(2) * Math.Pow(lambda, -2) * dlambda;
		Console.WriteLine($"T½ = {hl:F2} +/- {dhl:F2} days");
	}

	}
