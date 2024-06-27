using System;
using System.IO;
using static System.Math;

public static class Integration {

	static (double, double) Integrate(Func<double, double> f, double a, double b, double acc = 0.001, double eps = 0.001, double f2 = Double.NaN, double f3 = Double.NaN) {
		if (a == double.NegativeInfinity || b == double.PositiveInfinity) {
			return InfIntegrate(f, a, b, acc, eps);
		}
		double h = b - a;
		if (Double.IsNaN(f2)) {
			f2 = f(a + 2.0 / 6 * h);
			f3 = f(a + 4.0 / 6 * h);
		}
		double f1 = f(a + 1.0 / 6 * h);
		double f4 = f(a + 5.0 / 6 * h);
		double Q = (2 * f1 + f2 + f3 + 2 * f4) / 6 * h;
		double q = (f1 + f2 + f3 + f4) / 4 * h;
		double err = Abs(Q - q);
		if (err <= acc + eps * Abs(Q)) return (Q, err);
		else {
			var (result1, err1) = Integrate(f, a, (a + b) / 2, acc / Sqrt(2), eps, f1, f2);
			var (result2, err2) = Integrate(f, (a + b) / 2, b, acc / Sqrt(2), eps, f3, f4);
			return (result1 + result2, err1 + err2);
		}
	}

	static Func<double, double> TransformToCos(Func<double, double> f, double a, double b) {
		return theta => {
			double x = (a + b) / 2 + (b - a) / 2 * Cos(theta);
			return f(x) * Sin(theta) * (b - a) / 2;
		};
	}

	static Func<double, double> TransformToCosWithPi(Func<double, double> f) {
		return theta => f(Cos(theta)) * Sin(theta);
	}

	static Func<double, double> TransformEq62(Func<double, double> f) {
		return t => (f((1 - t) / t) + f(-(1 - t) / t)) / (t * t);
	}

	static Func<double, double> TransformEq64(Func<double, double> f, double a) {
		return t => f(a + t / (1 - t)) / Pow(1 - t, 2);
	}

	static Func<double, double> TransformEq66(Func<double, double> f, double b) {
		return t => f(b - (1 - t) / t) / (t * t);
	}

	static (double, double) InfIntegrate(Func<double, double> f, double a, double b, double acc, double eps) {
		if (a == double.NegativeInfinity && b == double.PositiveInfinity) {
			Func<double, double> transformedF = TransformEq62(f);
			return Integrate(transformedF, 0, 1, acc, eps);
		} else if (a == double.NegativeInfinity) {
			Func<double, double> transformedF = TransformEq66(f, b);
			return Integrate(transformedF, 0, 1, acc, eps);
		} else if (b == double.PositiveInfinity) {
			Func<double, double> transformedF = TransformEq64(f, a);
			return Integrate(transformedF, 0, 1, acc, eps);
		} else {
			throw new ArgumentException("Both a and b must be finite for this integration method.");
		}
	}

	static (double, double) CCIntegrate(Func<double, double> f, double a, double b, double acc = 0.001, double eps = 0.001) {
		if (a == double.NegativeInfinity || b == double.PositiveInfinity) {
			return InfIntegrate(f, a, b, acc, eps);
		}
		if (a == -1 && b == 1) {
			Func<double, double> transformedF = TransformToCosWithPi(f);
			return Integrate(transformedF, 0, PI, acc, eps);
		} else {
			Func<double, double> transformedF = TransformToCos(f, a, b);
			return Integrate(transformedF, 0, PI, acc, eps);
		}
	}

	static double Erf(double z) {
		if (z < 0) {
			return -Erf(-z);
		} else if (z <= 1) {
			Func<double, double> f = x => Exp(-x * x);
			var (result, _) = CCIntegrate(f, 0, z);
			return 2 / Sqrt(PI) * result;
		} else {
			Func<double, double> f = t => Exp(-Pow(z + (1 - t) / t, 2)) / (t * t);
			var (result, _) = CCIntegrate(f, 0, 1);
			return 1 - 2 / Sqrt(PI) * result;
		}
	}

	static void DoErfTest(string filename) {
		using (StreamWriter sw = File.CreateText(filename)) {
			for (double x = -3; x <= 3; x += 0.1) {
				double erfValue = Erf(x);
				sw.WriteLine($"{x}\t{erfValue}");
			}
		}
	}

	static void DoIntegrationTest() {
		Func<double, double> f1 = x => Sqrt(x);
		Func<double, double> f2 = x => 4 * Sqrt(1 - x * x);

		var (result1, err1) = CCIntegrate(f1, 0, 1);
		var (result2, err2) = CCIntegrate(f2, 0, 1);

		Console.WriteLine($"∫_0^1 dx √(x) = {result1}, Expected: 2/3, Error: {err1}");
		Console.WriteLine($"∫_0^1 dx 4√(1-x²) = {result2}, Expected: π, Error: {err2}");
	}

	static void DoTransformationTest() {
		int calls = 0;
		Func<double, double> invSqrt = x => { calls++; return 1 / Sqrt(x); };
		Func<double, double> logInvSqrt = x => { calls++; return Log(x) / Sqrt(x); };

		calls = 0;
		var (resultCC1, errCC1) = CCIntegrate(invSqrt, 0, 1);
		int ccCalls1 = calls;

		calls = 0;
		var (resultCC2, errCC2) = CCIntegrate(logInvSqrt, 0, 1);
		int ccCalls2 = calls;

		calls = 0;
		var (resultRegular1, errRegular1) = Integrate(invSqrt, 0, 1);
		int regularCalls1 = calls;

		calls = 0;
		var (resultRegular2, errRegular2) = Integrate(logInvSqrt, 0, 1);
		int regularCalls2 = calls;

		Console.WriteLine($"∫_0^1 dx 1/√(x) = 2");
		Console.WriteLine($"∫_0^1 dx 1/√(x) Regular result: {resultRegular1}, calls: {regularCalls1}, Error: {errRegular1}");
		Console.WriteLine($"∫_0^1 dx 1/√(x) Clenshaw–Curtis result: {resultCC1}, calls: {ccCalls1}, Error: {errCC1}");
		Console.WriteLine($"∫_0^1 dx ln(x)/√(x) = -4");
		Console.WriteLine($"∫_0^1 dx ln(x)/√(x) Regular result: {resultRegular2}, calls: {regularCalls2}, Error: {errRegular2}");
		Console.WriteLine($"∫_0^1 dx ln(x)/√(x) Clenshaw–Curtis result: {resultCC2}, calls: {ccCalls2}, Error: {errCC2}");
	}

	public static void DoInfiniteIntegrationTest() {
		int calls = 0;
		Func<double, double> inverseSquare = x => { calls++; return 1 / (x * x); };

		calls = 0;
		var (resultRegular, errRegular) = Integrate(inverseSquare, 1, double.PositiveInfinity);
		int regularCalls = calls;

		Console.WriteLine($"∫_1^∞ dx 1/x² = 1");
		Console.WriteLine($"∫_1^∞ dx 1/x² Result: {resultRegular}, calls: {regularCalls}, Error: {errRegular}");
	}

	static void Main() {
		DoIntegrationTest();
		DoErfTest("Out.Erf.txt");
		DoTransformationTest();
		DoInfiniteIntegrationTest();
	}
}
