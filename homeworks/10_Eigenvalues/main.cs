using System;
using static System.Math;

public static class Program {
	public static void TimesJ(Matrix A, int p, int q, double theta) {
		double c = Cos(theta), s = Sin(theta);
		for (int i = 0; i < A.Rows; i++) {
			double Aip = A[i, p], Aiq = A[i, q];
			A[i, p] = c * Aip - s * Aiq;
			A[i, q] = s * Aip + c * Aiq;
		}
	}

	public static void Jtimes(Matrix A, int p, int q, double theta) {
		double c = Cos(theta), s = Sin(theta);
		for (int j = 0; j < A.Cols; j++) {
			double Apj = A[p, j], Aqj = A[q, j];
			A[p, j] = c * Apj + s * Aqj;
			A[q, j] = -s * Apj + c * Aqj;
		}
	}

	public static (Vec, Matrix) Cyclic(Matrix M) {
		Matrix A = M.Copy();
		Matrix V = new Matrix(M.Rows, M.Rows);
		V.Identity();
		Vec w = new Vec(M.Rows);

		bool changed;
		int n = M.Rows;
		do {
			changed = false;
			for (int p = 0; p < n - 1; p++) {
				for (int q = p + 1; q < n; q++) {
					double Apq = A[p, q], App = A[p, p], Aqq = A[q, q];
					double theta = 0.5 * Atan2(2 * Apq, Aqq - App);
					double c = Cos(theta), s = Sin(theta);
					double new_App = c * c * App - 2 * s * c * Apq + s * s * Aqq;
					double new_Aqq = s * s * App + 2 * s * c * Apq + c * c * Aqq;
					if (new_App != App || new_Aqq != Aqq) {
						changed = true;
						TimesJ(A, p, q, theta); // A ← A * J
						Jtimes(A, p, q, -theta); // A ← J^T * A
						TimesJ(V, p, q, theta); // V ← V * J
					}
				}
			}
		} while (changed);

		for (int i = 0; i < M.Rows; i++) w[i] = A[i, i];
		return (w, V);
	}

	public static (Vec, Matrix) HydrogenWavefunction(double rmax, double dr) {
		int npoints = (int)(rmax / dr) - 1;
		Vec r = new Vec(npoints);
		for (int i = 0; i < npoints; i++) r[i] = dr * (i + 1);

		Matrix H = new Matrix(npoints, npoints);
		for (int i = 0; i < npoints - 1; i++) {
			H[i, i] = -2 * (-0.5 / (dr * dr));
			H[i, i + 1] = 1 * (-0.5 / (dr * dr));
			H[i + 1, i] = 1 * (-0.5 / (dr * dr));
		}
		H[npoints - 1, npoints - 1] = -2 * (-0.5 / (dr * dr));
		for (int i = 0; i < npoints; i++) H[i, i] += -1 / r[i];

		(Vec eigenvalues, Matrix eigenvectors) = Cyclic(H);

		return (eigenvalues, eigenvectors);
	}

	public static void Main(string[] args) {
		double rmax = 10;
		double dr = 0.3;
		bool printfunctions = false;
		bool hydrogenMode = false;

		foreach (string arg in args) {
			var words = arg.Split(':');
			if (words[0] == "-rmax") {rmax = double.Parse(words[1]); hydrogenMode = true;}
			if (words[0] == "-dr") {dr = double.Parse(words[1]); hydrogenMode = true;}
			if (words[0] == "-printfunctions") {printfunctions = true; hydrogenMode = true;}
		}

		if (hydrogenMode) {
			if (printfunctions) {
				RunHydrogenNumericalCalculation(rmax, dr);
			}
			else RunHydrogenConvergence(rmax, dr);
		}
		else RunTestA();
	}

	public static void RunHydrogenNumericalCalculation(double rmax, double dr) {
		(Vec epsilon, Matrix f) = HydrogenWavefunction(rmax, dr); 
		// epsilon.Print("Eigenvalues");
		f.Print();
	}

	public static void RunHydrogenConvergence(double rmax, double dr) {
		(Vec epsilon, Matrix f) = HydrogenWavefunction(rmax, dr);
		Console.WriteLine($"{rmax} {dr} {epsilon[0]}");
	}	

	public static void RunTestA() {
		int N = 5;
		Matrix A = new Matrix(N, N);
		A.FillRandomSymmetric(-5, 10);
		(Vec w, Matrix V) = Cyclic(A);
		Matrix D = new Matrix(N, N);
		for (int i = 0; i < N; i++) D[i, i] = w[i];

		/// A.Print("A:");
		/// V.Print("V:");
		/// w.Print("w:");

		Matrix VT = V.Transpose();

		Matrix VTAV = VT * A * V;
		if (Matrix.approx(VTAV, D, 1e-7, 1)) Console.WriteLine("V^T * A * V == D -> Passed");
		else Console.WriteLine("V^T * A * V == D -> Failed");

		Matrix VDVt = V * D * VT;
		if (Matrix.approx(VDVt, A, 1e-7, 1)) Console.WriteLine("V * D * V^T == A -> Passed");
		else Console.WriteLine("V * D * V^T == A -> Failed");

		Matrix I = new Matrix(V.Cols, V.Cols);
		I.Identity();
		Matrix VTV = VT * V;
		if (Matrix.approx(VTV, I, 1e-9, 1)) Console.WriteLine("V^T * V == I -> Passed");
		else Console.WriteLine("V^T * V == I -> Failed");

		Matrix VVT = V * VT;
		if (Matrix.approx(VVT, I, 1e-9, 1)) Console.WriteLine("V * V^T == I -> Passed");
		else Console.WriteLine("V * V^T == I -> Failed");

		/// VTAV.Print("VTAV:");
		/// D.Print("D:");
		/// VDVt.Print("VDVt:");
		/// A.Print("A:");
		/// Console.WriteLine($"VTAV[0,3] = {VTAV[0,3]}");
		/// Console.WriteLine($"D[0,3] = {D[0,3]}");
	}
}
