using System;
using System.Globalization;
using System.Collections.Generic;

/// <summary>
/// Represents a matrix with arbitrary rows and columns.
/// </summary>
public class Matrix {
	private int rows;
	private int cols;
	private double[,] data;

	/// <summary>
	/// Initializes a new instance of the Matrix class with the specified dimensions.
	/// </summary>
	/// <param name="rows">The number of rows.</param>
	/// <param name="cols">The number of columns.</param>
	public Matrix(int rows, int cols) {
		this.rows = rows;
		this.cols = cols;
		this.data = new double[rows, cols];
	}

	/// <summary>
	/// Gets or sets the element at the specified row and column.
	/// </summary>
	/// <param name="row">The row index.</param>
	/// <param name="col">The column index.</param>
	/// <returns>The element at the specified row and column.</returns>
	public double this[int row, int col] {
		get { return data[row, col]; }
		set { data[row, col] = value; }
	}

	/// <summary>
	/// Gets the number of rows.
	/// </summary>
	public int Rows { get { return rows; } }

	/// <summary>
	/// Gets the number of columns.
	/// </summary>
	public int Cols { get { return cols; } }

	/// <summary>
	/// Gets or sets the column vector at the specified index.
	/// </summary>
	/// <param name="colIndex">The column index.</param>
	/// <returns>The vector representing the specified column.</returns>
	public Vec this[int colIndex] {
		get {
			return Vec.FromMatrixColumn(this, colIndex);
		}
		set {
			value.ToMatrixColumn(this, colIndex);
		}
	}
	
	/// <summary>
	/// Creates a deep copy of the matrix.
	/// </summary>
	/// <returns>A new matrix that is a copy of the current matrix.</returns>
	public Matrix Copy() {
		Matrix copy = new Matrix(this.rows, this.cols);
		for (int i = 0; i < this.rows; i++) {
			for (int j = 0; j < this.cols; j++) {
				copy[i, j] = this[i, j];
			}
		}
		return copy;
	}

	/// <summary>
	/// Prints the matrix to the console.
	/// </summary>
	public void Print() {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				Console.Write($"{data[i, j],10:F4} ");
			}
			Console.WriteLine();
		}
	}

	public void Print(string s) {
		Console.WriteLine(s);
		this.Print();
	}
	/// <summary>
	/// Adds two matrices.
	/// </summary>
	/// <param name="A">The first matrix.</param>
	/// <param name="B">The second matrix.</param>
	/// <returns>The result of adding the two matrices.</returns>
	/// <exception cref="ArgumentException">Thrown when the matrices have different dimensions.</exception>
	public static Matrix operator +(Matrix A, Matrix B) {
		if (A.rows != B.rows || A.cols != B.cols) {
			throw new ArgumentException("The shapes of the two matrices do not match, which is necessary for addition");
		}

		Matrix result = new Matrix(A.rows, A.cols);
		for (int i = 0; i < A.rows; i++) {
			for (int j = 0; j < A.cols; j++) {
				result[i, j] = A[i, j] + B[i, j];
			}
		}
		return result;
	}

	/// <summary>
	/// Multiplies two matrices.
	/// </summary>
	/// <param name="A">The first matrix.</param>
	/// <param name="B">The second matrix.</param>
	/// <returns>The result of multiplying the two matrices.</returns>
	public static Matrix operator *(Matrix A, Matrix B) {
		if (A.cols != B.rows) {
			throw new ArgumentException("Number of columns in A must equal number of rows in B.");
		}
		Matrix result = new Matrix(A.rows, B.cols);
		for (int i = 0; i < A.rows; i++) {
			for (int j = 0; j < B.cols; j++) {
				double sum = 0;
				for (int k = 0; k < A.cols; k++) {
					sum += A[i, k] * B[k, j];
				}
				result[i, j] = sum;
			}
		}
		return result;
	}

	/// <summary>
	/// Multiplies a matrix by a scalar.
	/// </summary>
	/// <param name="A">The matrix.</param>
	/// <param name="x">The scalar.</param>
	/// <returns>The result of multiplying the matrix by the scalar.</returns>
	public static Matrix operator *(Matrix A, double x) {
		Matrix result = new Matrix(A.rows, A.cols);
		for (int i = 0; i < A.rows; i++) {
			for (int j = 0; j < A.cols; j++) {
				result[i, j] = A[i, j] * x;
			}
		}
		return result;
	}

	/// <summary>
	/// Multiplies a scalar by a matrix.
	/// </summary>
	/// <param name="x">The scalar.</param>
	/// <param name="A">The matrix.</param>
	/// <returns>The result of multiplying the scalar by the matrix.</returns>
	public static Matrix operator *(double x, Matrix A) {
		return A * x;
	}

	/// <summary>
	/// Fills the matrix with a specific value.
	/// </summary>
	/// <param name="value">The value to fill the matrix with.</param>
	public void Fill(double value) {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				data[i, j] = value;
			}
		}
	}

	/// <summary>
	/// Fills the matrix with random values within the specified interval [lowerBoundary, upperBoundary].
	/// </summary>
	/// <param name="lowerBoundary">The lower boundary of the interval.</param>
	/// <param name="upperBoundary">The upper boundary of the interval.</param>
	public void FillRandom(double lowerBoundary, double upperBoundary) {
		Random rand = new Random();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				data[i, j] = rand.NextDouble() * (upperBoundary - lowerBoundary) + lowerBoundary;
			}
		}
	}

	/// <summary>
	/// Fills the matrix with random values within the specified interval [lowerBoundary, upperBoundary]
	/// ensuring the matrix is symmetric.
	/// </summary>
	/// <param name="lowerBoundary">The lower boundary of the interval.</param>
	/// <param name="upperBoundary">The upper boundary of the interval.</param>
	public void FillRandomSymmetric(double lowerBoundary, double upperBoundary) {
		if (rows != cols) {
			throw new InvalidOperationException("Matrix must be square to fill symmetrically.");
		}

		Random rand = new Random();
		for (int i = 0; i < rows; i++) {
			for (int j = i; j < cols; j++) {
				double value = rand.NextDouble() * (upperBoundary - lowerBoundary) + lowerBoundary;
				data[i, j] = value;
				if (i != j) {
					data[j, i] = value;
				}
			}
		}
	}

	/// <summary>
	/// Converts the matrix to an identity matrix.
	/// </summary>
	public void Identity() {
		if (rows != cols) {
			throw new ArgumentException("Cannot make non-square identity matrix");
		}

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				data[i, j] = (i == j) ? 1.0 : 0.0;
			}
		}
	}

	/// <summary>
	/// Transposes the matrix.
	/// </summary>
	/// <returns>The transposed matrix.</returns>
	public Matrix Transpose() {
		Matrix result = new Matrix(cols, rows);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				result[j, i] = data[i, j];
			}
		}
		return result;
	}

	public static Matrix Loadtxt() {
		List<double> values = new List<double>();
		int cols = 0, rows = 0;

		string line;
		while ((line = Console.ReadLine()) != null) {
			// Handle different delimiters
			char[] delimiters = { ' ', ',', '\t' };
			string[] parts = line.Split(delimiters, StringSplitOptions.RemoveEmptyEntries);

			if (cols == 0) {
				cols = parts.Length;
			} else if (parts.Length != cols) {
				throw new InvalidOperationException("Inconsistent number of columns in input file.");
			}

			foreach (var part in parts) {
				values.Add(double.Parse(part, CultureInfo.InvariantCulture));
			}

			rows++;
		}

		Matrix matrix = new Matrix(rows, cols);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				matrix[i, j] = values[i * cols + j];
			}
		}

		return matrix;
	}

	/// <summary>
	/// Checks if two doubles are approximately equal.
	/// </summary>
	/// <param name="x">The first double.</param>
	/// <param name="y">The second double.</param>
	/// <param name="acc">The absolute tolerance.</param>
	/// <param name="eps">The relative tolerance.</param>
	/// <returns>True if the doubles are approximately equal, otherwise false.</returns>
	public static bool approx(double x, double y, double acc = 1e-9, double eps = 1e-9) {
		if (Math.Abs(x - y) > acc) return false;
		if (x != 0 && y != 0 && Math.Abs(x - y) / (Math.Abs(x) + Math.Abs(y)) > eps) return false;
		return true;
	}

	/// <summary>
	/// Checks if two matrices are approximately equal.
	/// </summary>
	/// <param name="A">The first matrix.</param>
	/// <param name="B">The second matrix.</param>
	/// <param name="acc">The absolute tolerance.</param>
	/// <param name="eps">The relative tolerance.</param>
	/// <returns>True if the matrices are approximately equal, otherwise false.</returns>
	public static bool approx(Matrix A, Matrix B, double acc = 1e-9, double eps = 1e-9) {
		if (A.rows != B.rows || A.cols != B.cols) {
			throw new ArgumentException("The shapes of the matrices do not match, so they cannot be compared.");
		}
		for (int i = 0; i < A.rows; i++) {
			for (int j = 0; j < A.cols; j++) {
				if (!approx(A[i, j], B[i, j], acc, eps)) return false;
			}
		}
		return true;
	}

	/// <summary>
	/// Checks if this matrix is approximately equal to another matrix.
	/// </summary>
	/// <param name="other">The other matrix.</param>
	/// <param name="acc">The absolute tolerance.</param>
	/// <param name="eps">The relative tolerance.</param>
	/// <returns>True if the matrices are approximately equal, otherwise false.</returns>
	public bool approx(Matrix other, double acc = 1e-9, double eps = 1e-9) {
		return approx(this, other, acc, eps);
	}

	/// <summary>
	/// Checks if the matrix is approximately upper triangular within a given tolerance.
	/// </summary>
	/// <param name="tolerance">The tolerance for the approximation.</param>
	/// <returns>True if the matrix is approximately upper triangular, otherwise false.</returns>
	public bool IsUpperTriangular(double tolerance = 1e-9) {
		for (int i = 1; i < rows; i++) {
			for (int j = 0; j < i; j++) {
				if (Math.Abs(data[i, j]) > tolerance) {
					return false;
				}
			}
		}
		return true;
	}
	
	/// <summary>
	/// Checks if the matrix is approximately lower triangular within a given tolerance.
	/// </summary>
	/// <param name="tolerance">The tolerance for the approximation.</param>
	/// <returns>True if the matrix is approximately lower triangular, otherwise false.</returns>
	public bool IsLowerTriangular(double tolerance = 1e-9) {
		for (int i = 0; i < rows; i++) {
			for (int j = i + 1; j < cols; j++) {
				if (Math.Abs(data[i, j]) > tolerance) {
					return false;
				}
			}
		}
		return true;
	}
	
}

