#Question 1

import numpy as np
import matplotlib.pyplot as plt

def true_function(x):
    return np.sin(x)

def taylor_approximation(x):
    return x - (x**3) / 6 + (x**5) / 120

def approx_1(x):
    return x / (1 + x**2 / 6 + 0.019 * x**4)

def approx_2(x):
    return x - 0.117 * x**3 / (1 + x**2 / 20)

def approx_3(x):
    return (x - 0.142 * x**3 )/ (1 + x**2 / 40)

# Define the range of x values
x_values = np.linspace(0, 5, 100)

# Compute the true values and approximation values
true_values = true_function(x_values)
taylor_values = taylor_approximation(x_values)
approx_1_values = approx_1(x_values)
approx_2_values = approx_2(x_values)
approx_3_values = approx_3(x_values)

# Calculate the errors
error_taylor = np.abs(true_values - taylor_values)
error_approx_1 = np.abs(true_values - approx_1_values)
error_approx_2 = np.abs(true_values - approx_2_values)
error_approx_3 = np.abs(true_values - approx_3_values)

# Plot the results
plt.figure(figsize=(10, 6))

plt.plot(x_values, true_values, label='True Function (sin(x))')
plt.plot(x_values, taylor_values, label='Taylor Approximation')
plt.plot(x_values, approx_1_values, label='Approximation 1')
plt.plot(x_values, approx_2_values, label='Approximation 2')
plt.plot(x_values, approx_3_values, label='Approximation 3')

plt.title('Comparison of Sin(x) Approximations')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()

plt.figure(figsize=(10, 6))

plt.semilogy(x_values, error_taylor, label='Taylor Approximation Error')
plt.semilogy(x_values, error_approx_1, label='Approximation 1 Error')
plt.semilogy(x_values, error_approx_2, label='Approximation 2 Error')
plt.semilogy(x_values, error_approx_3, label='Approximation 3 Error')

plt.title('Error Comparison of Sin(x) Approximations')
plt.xlabel('x')
plt.ylabel('Absolute Error')
plt.legend()

plt.show()




#Question 3

from scipy.integrate import quad

def f(x):
    return 1 / (1 + x**2)

def composite_trapezoidal(a, b, n):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    result = (f(a) + f(b)) / 2
    result += np.sum(f(x[1:-1]))
    result *= h
    return result

def composite_simpsons(a, b, n):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    result = f(a) + f(b)
    result += 4 * np.sum(f(x[1::2]))
    result += 2 * np.sum(f(x[2:-1:2]))
    result *= h / 3
    return result

# Function to find n for Trapezoidal Rule
def find_n_trapezoidal(tol):
    n = 1
    while True:
        result = composite_trapezoidal(-5, 5, n)
        error = abs(quad(f, -5, 5)[0] - result)
        if error < tol:
            return n
        n *= 2

# Function to find n for Simpson's Rule
def find_n_simpsons(tol):
    n = 2
    while True:
        result = composite_simpsons(-5, 5, n)
        error = abs(quad(f, -5, 5)[0] - result)
        if error < tol:
            return n
        n *= 2

# Finding n for both methods
n_trapezoidal = find_n_trapezoidal(1e-4)
n_simpsons = find_n_simpsons(1e-4)

# Computing values using both methods
result_trapezoidal = composite_trapezoidal(-5, 5, n_trapezoidal)
result_simpsons = composite_simpsons(-5, 5, n_simpsons)

# Comparing with scipy's quad routine
result_quad_default, _ = quad(f, -5, 5)
result_quad_custom, _ = quad(f, -5, 5, epsabs=1e-4)

print(f"Composite Trapezoidal Rule Result (n={n_trapezoidal}): {result_trapezoidal}")
print(f"Composite Simpson's Rule Result (n={n_simpsons}): {result_simpsons}")
print(f"Scipy Quad Default Result: {result_quad_default}")
print(f"Scipy Quad Custom Tolerance Result: {result_quad_custom}")

