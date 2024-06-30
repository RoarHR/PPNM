using System;
using static System.Math;

public class QRGS {
	public static (Matrix, Matrix) Decomp(Matrix A) {
		int n = A.Rows;
		int m = A.Cols;
		Matrix Q = A.Copy();
		Matrix R = new Matrix(m, m);
		for (int i=0; i<m; i++) {
			double norm = Q[i].Norm();
			R[i, i] = norm;
			Q[i] /= norm;
			for (int j = i+1; j<m; j++) {
				R[i, j]=Q[i].Dot(Q[j]);
				Q[j] = Q[j] - (Q[i] * R[i, j]);
				/// Console.WriteLine($"Q[{j}]':");
				/// Q[j].Print();
			}
		}
		return (Q, R);
	} // Decomp
	
	public static Vec Backsub(Matrix R, Vec x) {
		Vec c = x.Copy();
		for (int i = c.Length - 1; i >= 0; i--) {
			double sum = 0;
			for (int k = i + 1; k < c.Length; k++) sum += R[i, k] * c[k];
			c[i] = (c[i] - sum) / R[i, i];
		}
		return c;
	} // Backsub
	
	public static Vec Solve(Matrix Q, Matrix R, Vec b) {
		Vec x = Backsub(R, Q.Transpose() * b);
		return x;
	} //Solve
	
	public static double Det(Matrix A) {
		if (A.Rows != A.Cols) {
			throw new ArgumentException("Cannot find determinant for non-square matrix.");
		}
		double result = 1;
		for (int i=0; i<= A.Rows-1; i++) {
			result *= A[i, i];
		}
		return result;
	} // Det
	
	public static Matrix Inverse(Matrix Q, Matrix R) {
		int n = Q.Cols;
		Matrix B = new Matrix(n, n);
		Vec ei = new Vec(n);
		ei.Fill(0);
		for (int i=0; i<n; i++) {
			ei[i] = 1;
			B[i] = Solve(Q, R, ei);
			ei[i] = 0;
		}
		return B;
	} // Inverse
		
	public static Matrix Inverse(Matrix A) {
		(Matrix Q, Matrix R) = Decomp(A);
		Matrix B = Inverse(Q, R);
		return B;
	}
	
	public static int Main(string[] args) {
		bool timingMode = false;
		int matrixSize = 10;

		foreach (string arg in args) {
			var words = arg.Split(':');
			if (words[0] == "-matrixsize") {
				timingMode = true;
				matrixSize = int.Parse(words[1]);
			}
		}

		if (timingMode) {
			RunQRDecomposition(matrixSize);
		} else {
			RunTest1();
			RunTest2();
			RunTest3();
		}
		return 0;
	} // Main

	static void RunTest1() {
		Matrix A = new Matrix(4, 4);
		A.FillRandom(-10, 20);
		(Matrix Q, Matrix R) = Decomp(A);
		if (R.IsUpperTriangular()) Console.WriteLine("R is upper triangular -> Passed");
		else Console.WriteLine("R is upper triangular -> Failed");

		Matrix I = new Matrix(Q.Cols, Q.Cols);
		I.Identity();
		if (Matrix.approx(Q.Transpose() * Q, I)) Console.WriteLine("Q^TQ = I -> Passed");
		else Console.WriteLine("Q^TQ = I -> Failed");

		if (Matrix.approx(Q*R, A)) Console.WriteLine("QR = A -> Passed");
		else Console.WriteLine("QR = A -> Failed");

		/// Console.WriteLine("A:");
		/// A.Print();
		/// Console.WriteLine("A[0]:");
		/// A[0].Print();
		/// Console.WriteLine("Q:");
		/// Q.Print();
		/// Console.WriteLine("R:");
		/// R.Print();
		/// Console.WriteLine("QQ^T:");
		/// Matrix QQT = Q * Q.Transpose();
		/// QQT.Print();
		/// Console.WriteLine("QR:");
		/// Matrix QR = Q * R;
		/// QR.Print();
	} // Test 1
	
	static void RunTest2() {
		Matrix A = new Matrix(10, 10);
		Vec b = new Vec(10);
		A.FillRandom(0, 10);
		b.FillRandom(-1, 1);

		(Matrix Q, Matrix R) = Decomp(A);
		Vec x = Solve(Q, R, b);

		if (Vec.Approx(A*x, b)) Console.WriteLine("Ax = b -> Passed");
		else Console.WriteLine("Ax = b -> Failed");
	} // Test 2

	static void RunTest3() {
		int n = 10;
		Matrix A = new Matrix(n, n);
		A.FillRandom(-5, 5);
		(Matrix Q, Matrix R) = Decomp(A);
		Matrix B = Inverse(Q, R);
		Matrix I = new Matrix(n, n);
		I.Identity();

		if (Matrix.approx(A*B, I)) Console.WriteLine("AA⁻¹ = I -> Passed");
		else Console.WriteLine("AA⁻¹ = I -> Failed");
	} // Test 3
	
	static void RunQRDecomposition(int matrixSize) {
		Matrix A = new Matrix(matrixSize, matrixSize);
		A.FillRandom(-10, 10);
		(Matrix Q, Matrix R) = Decomp(A);
		if (Q[0, 0] + R[0, 0] == 0) Console.WriteLine(""); // To avoid warning
		Console.WriteLine($"Decomposed matrix N={matrixSize}");
	} // Time test
} // QRGS
