import numpy as np
import matplotlib.pyplot as plt
#Part 1
from scipy.special import erf

# Constants
Ts = -15.0
Ti = 20.0
alpha = 0.138e-6
exposure_time_seconds = 60 * 24 * 3600
tolerance = 1e-13
max_iterations = 1000

# Define the function f(x) and its derivative f'(x)
def f(x):
    return Ts + (Ti - Ts) * erf(x / (2 * np.sqrt(alpha * exposure_time_seconds)))

def df(x):
    h = 1e-6  # Small step for numerical differentiation
    return (f(x + h) - f(x)) / h

# Create a range of x values for plotting
x_values = np.linspace(0, 100, 1000)  # Adjust the range as needed

# Calculate the corresponding f(x) values
f_values = f(x_values)

# Plot f(x)
plt.figure(figsize=(10, 6))
plt.plot(x_values, f_values, label='f(x)')
plt.xlabel('Depth (x)')
plt.ylabel('f(x)')
plt.title('Plot of f(x)')
plt.axhline(0, color='r', linestyle='--', linewidth=0.7)  # Horizontal line at y=0
plt.grid(True)
plt.legend()
plt.show()


# Bisection Method
def bisection_method(a, b, tol, max_iter):
    # Initialize variables
    iteration = 0
    x_approx = (a + b) / 2.0

    while abs(f(x_approx)) > tol and iteration < max_iter:
        if f(x_approx) * f(a) < 0:
            b = x_approx
        else:
            a = x_approx
        x_approx = (a + b) / 2.0
        iteration += 1

    return x_approx


# Newton's Method
def newton_method(x0, tol, max_iter):
    # Initialize variables
    iteration = 0

    while abs(f(x0)) > tol and iteration < max_iter:
        x0 = x0 - f(x0) / df(x0)
        iteration += 1

    return x0


# Calculate x_bar (the upper bound of the depth)
x_bar = 100.0  # Adjust as needed

# (c) Compute depth using Bisection Method
a_bisection = 0.0
b_bisection = x_bar
depth_bisection = bisection_method(a_bisection, b_bisection, tolerance, max_iterations)
print(f"Depth (Bisection Method): {depth_bisection} meters")

# (d) Compute depth using Newton's Method
x0_newton = 0.01  # Starting value for Newton's Method
depth_newton = newton_method(x0_newton, tolerance, max_iterations)
print(f"Depth (Newton's Method): {depth_newton} meters")

# Try starting Newton's Method with x0 = x_bar
depth_newton_x_bar = newton_method(x_bar, tolerance, max_iterations)
print(f"Depth (Newton's Method with x0 = x_bar): {depth_newton_x_bar} meters")

#Part 4

import sympy as sp

def newton_method(f, x0, tol, max_iter):
    # Initialize variables
    iteration = 0
    f_prime = sp.diff(f, x)
    while abs(f.subs(x0)) > tol and iteration < max_iter:
        x0 = x0 - f.subs(x0) / f_prime.subs(x0)
        iteration += 1

    return x0
x = np.linspace(0, 100, 1000)  
f = sp.exp(3*x) - 27*x**6 + 27*x**4*sp.exp(x) - 9*x**2*sp.exp(2*x)
root = newton_method(f, 4.0,tolerance,max_iterations)
#if root is not None:
    #print(f"Root found: {root}")
#else:
    #print("Method did not converge.")
#When uncommented, these lines result in a division by 0 error
def modified_newton_class(f, x0, tol=1e-6, max_iter=100):
    x = sp.symbols('x')
    f_prime = sp.diff(f, x)
    f_double_prime = sp.diff(f_prime, x)
    x_k = x0

    for _ in range(max_iter):
        f_val = f.subs(x, x_k)
        f_prime_val = f_prime.subs(x, x_k)

        if abs(f_val) < tol:
            return x_k

        if abs(f_prime_val) < tol:
            return None  # Avoid division by zero

        x_k = x_k - f_val / f_prime_val

    return None  # Convergence not achieved

# Usage example:
x = sp.symbols('x')
f = sp.exp(3*x) - 27*x**6 + 27*x**4*sp.exp(x) - 9*x**2*sp.exp(2*x)
root = modified_newton_class(f, x0=4.0)
if root is not None:
    print(f"Root found: {root}")
else:
    print("Method did not converge.")


def modified_newton_problem2(f, x0, tol=1e-6, max_iter=100):
    x = sp.symbols('x')
    f_prime = sp.diff(f, x)
    x_k = x0

    for _ in range(max_iter):
        f_val = f.subs(x, x_k)
        f_prime_val = f_prime.subs(x, x_k)

        if abs(f_val) < tol:
            return x_k

        x_k = x_k - f_val / f_prime_val

    return None  # Convergence not achieved

# Usage example:
x = sp.symbols('x')
f = sp.exp(3*x) - 27*x**6 + 27*x**4*sp.exp(x) - 9*x**2*sp.exp(2*x)
root = modified_newton_problem2(f, x0=4.0)
if root is not None:
    print(f"Root found: {root}")
else:
    print("Method did not converge.")


# Part 5
def f(x):
    return x**6 - x - 1

def f_prime(x):
    return 6 * x**5 - 1

# Newton's method
def newton_method(x0, tol, max_iter):
    x = x0
    error_history = []

    for i in range(max_iter):
        x_next = x - f(x) / f_prime(x)
        error = abs(x_next - x)
        error_history.append(error)

        if error < tol:
            break

        x = x_next

    return x, error_history

# Secant method
def secant_method(x0, x1, tol, max_iter):
    x_prev = x0
    x = x1
    error_history = []

    for i in range(max_iter):
        x_next = x - f(x) * (x - x_prev) / (f(x) - f(x_prev))
        error = abs(x_next - x)
        error_history.append(error)

        if error < tol:
            break

        x_prev = x
        x = x_next

    return x, error_history

# Initial values
x0_newton = 2
x0_secant = 2
x1_secant = 1
tolerance = 1e-10
max_iterations = 100

# Apply Newton's method and Secant method
root_newton, error_newton = newton_method(x0_newton, tolerance, max_iterations)
root_secant, error_secant = secant_method(x0_secant, x1_secant, tolerance, max_iterations)

# Display the results
print(f"Newton's Method: Approximated Root = {root_newton}")
print(f"Secant Method: Approximated Root = {root_secant}")

# (a) Create a table of errors
print("\nError Table:")
print("Iteration\tNewton's Error\tSecant Error")
for i in range(min(len(error_newton), len(error_secant))):
    print(f"{i + 1}\t\t{error_newton[i]}\t\t{error_secant[i]}")

# (b) Plot |xk+1 - α| vs |xk - α| on log-log axes
exact_root = 1.1347241386  # Exact root obtained from a numerical solver

plt.figure(figsize=(10, 6))
plt.loglog(error_newton[:-1], error_newton[1:], 'bo-', label="Newton's Method")
plt.loglog(error_secant[:-1], error_secant[1:], 'ro-', label="Secant Method")
plt.xlabel('|x_k+1 - α|')
plt.ylabel('|x_k - α|')
plt.title('Convergence Plot: |x_k+1 - α| vs |x_k - α|')
plt.legend()
plt.grid(True)

# Show the plot
plt.show()

