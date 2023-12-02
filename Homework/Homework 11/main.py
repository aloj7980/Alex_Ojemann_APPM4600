#Question 3 Part a

import numpy as np
from scipy.special import gamma
from scipy.integrate import trapz


def gamma_trapezoidal(x):
    # Define the upper limit for truncating the interval
    upper_limit = 50

    # Define the function to integrate
    def integrand(t):
        return t ** (x - 1) * np.exp(-t)

    # Create the array of points for trapezoidal rule
    t_values = np.linspace(0, upper_limit, 1000)
    y_values = integrand(t_values)

    # Compute the integral using composite trapezoidal rule
    result = trapz(y_values, t_values)

    return result


# Compute gamma function for x = 2, 4, 6, 8, 10 using trapezoidal rule
x_values = [2, 4, 6, 8, 10]
for x in x_values:
    gamma_approx = gamma_trapezoidal(x)
    gamma_scipy = gamma(x)
    print(f"Gamma({x}) (Trapezoidal): {gamma_approx:.8f}, SciPy Gamma({x}): {gamma_scipy:.8f}")


#Part b

from scipy.integrate import quad

# Define the function to integrate
def integrand(t, x):
    return t**(x-1) * np.exp(-t)

def gamma_adaptive_quad(x):
    # Perform adaptive quadrature using quad function
    result, _ = quad(integrand, 0, np.inf, args=(x,))
    return result

# Compute gamma function for x = 2, 4, 6, 8, 10 using quad
x_values = [2, 4, 6, 8, 10]
for x in x_values:
    gamma_approx = gamma_adaptive_quad(x)
    gamma_scipy = gamma(x)
    print(f"Gamma({x}) (Adaptive Quad): {gamma_approx:.8f}, SciPy Gamma({x}): {gamma_scipy:.8f}")


#Part c

from scipy.special import roots_laguerre


def gamma_gauss_laguerre(x):
    # Obtain weights and abscissae for Gauss-Laguerre quadrature
    weights, abscissae = roots_laguerre(x)

    # Define the function to integrate
    def integrand(t):
        return t ** (x - 1)

    # Perform Gauss-Laguerre quadrature
    result = sum(weights * integrand(abscissae))

    return result


# Compute gamma function for x = 2, 4, 6, 8, 10 using Gauss-Laguerre quadrature
x_values = [2, 4, 6, 8, 10]
for x in x_values:
    gamma_approx = gamma_gauss_laguerre(x)
    gamma_scipy = gamma(x)
    print(f"Gamma({x}) (Gauss-Laguerre): {gamma_approx:.8f}, SciPy Gamma({x}): {gamma_scipy:.8f}")

