using System;

public class Vec {
	private int length;
	private double[] data;

	public Vec(int length) {
		this.length = length;
		this.data = new double[length];
	}

	public double this[int index] {
		get { return data[index]; }
		set { data[index] = value; }
	}

	public int Length {
		get { return length; }
	}

	public Vec Append(double z) {
		Vec x = new Vec(this.length + 1);
		for (int i = 0; i < this.length; i++) {
			x[i] = this[i];
		}
		x[this.length] = z;
		return x;
	}

	public Vec Copy() {
		Vec copy = new Vec(this.length);
		for (int i = 0; i < this.length; i++) {
			copy[i] = this[i];
		}
		return copy;
	}

	public void Print() {
		for (int i = 0; i < this.length; i++) {
			Console.Write($"{this[i],10:F4} ");
		}
		Console.WriteLine();
	}

	public void Print(string s) {
		Console.WriteLine(s);
		this.Print();
	}

	// Vector addition
	public static Vec operator +(Vec A, Vec B) {
		if (A.length != B.length) {
			throw new ArgumentException("Vector lengths do not match.");
		}

		Vec result = new Vec(A.length);
		for (int i = 0; i < A.length; i++) {
			result[i] = A[i] + B[i];
		}
		return result;
	}

	// Vector subtraction
	public static Vec operator -(Vec A) {
		Vec result = new Vec(A.length);
		for (int i = 0; i < A.length; i++) {
			result[i] = -A[i];
		}
		return result;
	}

	public static Vec operator -(Vec A, Vec B) {
		if (A.length != B.length) {
			throw new ArgumentException("Vector lengths do not match.");
		}

		Vec result = new Vec(A.length);
		for (int i = 0; i < A.length; i++) {
			result[i] = A[i] - B[i];
		}
		return result;
	}


	// Scalar multiplication
	public static Vec operator *(Vec A, double scalar) {
		Vec result = new Vec(A.length);
		for (int i = 0; i < A.length; i++) {
			result[i] = A[i] * scalar;
		}
		return result;
	}

	public static Vec operator *(double scalar, Vec A) {
		return A * scalar;
	}

	public static Vec operator /(Vec A, double scalar) {
		return A * (1 / scalar);
	}

	// Matrix-vector multiplication
	public static Vec operator *(Matrix matrix, Vec vec) {
		if (matrix.Cols != vec.Length) {
			throw new ArgumentException("Number of columns in the matrix must equal the length of the vector.");
		}

		Vec result = new Vec(matrix.Rows);
		for (int i = 0; i < matrix.Rows; i++) {
			double sum = 0;
			for (int j = 0; j < matrix.Cols; j++) {
				sum += matrix[i, j] * vec[j];
			}
			result[i] = sum;
		}
		return result;
	}

	public static Vec operator *(Vec vec, Matrix matrix) {
		if (vec.Length != matrix.Rows) {
			throw new ArgumentException("Length of the vector must equal the number of rows in the matrix.");
		}

		Vec result = new Vec(matrix.Cols);
		for (int i = 0; i < matrix.Cols; i++) {
			double sum = 0;
			for (int j = 0; j < matrix.Rows; j++) {
				sum += vec[j] * matrix[j, i];
			}
			result[i] = sum;
		}
		return result;
	}

	// Dot product
	public static double Dot(Vec A, Vec B) {
		if (A.length != B.length) {
			throw new ArgumentException("Vector lengths do not match.");
		}

		double sum = 0;
		for (int i = 0; i < A.length; i++) {
			sum += A[i] * B[i];
		}
		return sum;
	}

	public double Dot(Vec other) {
		return Dot(this, other);
	}

	// Norm
	public double Norm() {
		double normsquared = Dot(this, this);
		return Math.Sqrt(normsquared);
	}

	// Fill vector with a specific value
	public void Fill(double value) {
		for (int i = 0; i < length; i++) {
			data[i] = value;
		}
	}

	public void FillRandom(double boundaryA, double boundaryB) {
		Random rand = new Random();
		for (int i = 0; i < length; i++) {
			data[i] = rand.NextDouble() * (boundaryB - boundaryA) + boundaryA;
		}
	}

	// Check approximate equality
	public static bool Approx(Vec A, Vec B, double acc = 1e-9, double eps = 1e-9) {
		if (A.length != B.length) {
			return false;
		}

		for (int i = 0; i < A.length; i++) {
			if (Math.Abs(A[i] - B[i]) > acc) return false;
			if (Math.Abs(A[i] - B[i]) / (Math.Abs(A[i]) + Math.Abs(B[i])) > eps) return false;
		}
		return true;
	}

	public bool Approx(Vec other, double acc = 1e-9, double eps = 1e-9) {
		return Approx(this, other, acc, eps);
	}

	public override string ToString() {
		return $"[{string.Join(", ", data)}]";
	}

	// Method to extract a column vector from a matrix
	public static Vec FromMatrixColumn(Matrix matrix, int colIndex) {
		if (colIndex < 0 || colIndex >= matrix.Cols) {
			throw new ArgumentOutOfRangeException(nameof(colIndex), "Column index is out of range.");
		}
		Vec vec = new Vec(matrix.Rows);
		for (int i = 0; i < matrix.Rows; i++) {
			vec[i] = matrix[i, colIndex];
		}
		return vec;
	}

	// Method to set a column in a matrix from a vector
	public void ToMatrixColumn(Matrix matrix, int colIndex) {
		if (colIndex < 0 || colIndex >= matrix.Cols) {
			throw new ArgumentOutOfRangeException(nameof(colIndex), "Column index is out of range.");
		}
		if (this.length != matrix.Rows) {
			throw new ArgumentException("Vector length does not match the number of rows in the matrix.", nameof(matrix));
		}
		for (int i = 0; i < matrix.Rows; i++) {
			matrix[i, colIndex] = this[i];
		}
	}

	// Method to extract a row vector from a matrix
	public static Vec FromMatrixRow(Matrix matrix, int rowIndex) {
		if (rowIndex < 0 || rowIndex >= matrix.Rows) {
			throw new ArgumentOutOfRangeException(nameof(rowIndex), "Row index is out of range.");
		}
		Vec vec = new Vec(matrix.Cols);
		for (int i = 0; i < matrix.Cols; i++) {
			vec[i] = matrix[rowIndex, i];
		}
		return vec;
	}

	// Method to set a row in a matrix from a vector
	public void ToMatrixRow(Matrix matrix, int rowIndex) {
		if (rowIndex < 0 || rowIndex >= matrix.Rows) {
			throw new ArgumentOutOfRangeException(nameof(rowIndex), "Row index is out of range.");
		}
		if (this.length != matrix.Cols) {
			throw new ArgumentException("Vector length does not match the number of columns in the matrix.", nameof(matrix));
		}
		for (int i = 0; i < matrix.Cols; i++) {
			matrix[rowIndex, i] = this[i];
		}
	}
}
