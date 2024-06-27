import numpy as np
from scipy.integrate import quad

def count_calls(f):
	def wrapped_f(*args, **kwargs):
		wrapped_f.calls += 1
		return f(*args, **kwargs)
	wrapped_f.calls = 0
	return wrapped_f

# Define the functions
@count_calls
def inv_sqrt(x):
	return 1 / np.sqrt(x)

@count_calls
def log_inv_sqrt(x):
	return np.log(x) / np.sqrt(x)

# Perform the integrations using SciPy's quad
result1, _ = quad(inv_sqrt, 0, 1)
result2, _ = quad(log_inv_sqrt, 0, 1)

# Print the results and number of function calls
print("Python/numpy tests:")
print(f"∫_0^1 dx 1/√(x) result: {result1}, calls: {inv_sqrt.calls}")
print(f"∫_0^1 dx ln(x)/√(x) result: {result2}, calls: {log_inv_sqrt.calls}")

calls = 0
def inverse_square(x):
    global calls
    calls += 1
    return 1 / (x * x)

calls = 0
result, err = quad(inverse_square, 1, np.inf)
print(f"∫_1^∞ dx 1/x² = 1")
print(f"∫_1^∞ dx 1/x² Numpy result: {result}, calls: {calls}, Error: {err}")
