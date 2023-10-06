import numpy as np
from numpy.linalg import norm, inv

def Iterate(x0, tol, Nmax):
    for its in range(Nmax):
        J = evalJ(x0)
        F = evalF(x0[0],x0[1])
        G = evalG(x0[0], x0[1])
        x1 = x0 - J.dot([F,G])
        print(x1)
        if norm(x1 - x0) < tol:
            xstar = x1
            ier = 0
            return [xstar, ier, its]
        x0 = x1
    xstar = x1
    ier = 1
    return [xstar, ier, its]


def evalJ(x):
    # Calculate and return the matrix at point x
    return np.array([[1/6,1/18],[0,1/6]])

def evalF(x,y):
    # Calculate and return the function vector at point x
    return 3*x**2 - y**2

def evalG(x,y):
    # Calculate and return the function vector at point x
    return 3*x*y**2 - x**3 - 1

# Initial guess, tolerance, and maximum iterations
initial_guess = np.array([1,1])
tolerance = 1e-6
max_iterations = 100

# Call Newton's method to find the root
root, error_code, iterations = Iterate(initial_guess, tolerance, max_iterations)

if error_code == 0:
    print(f"Approximate root: {root}")
    print(f"Iterations: {iterations}")
else:
    print("Method did not converge within the maximum number of iterations.")

def evalJ(x):
    # Calculate and return the matrix at point x
    return np.array([[6*x[0],-2*x[1]],[3*x[1]**2-3*x[1]**2,6*x[0]*x[1]]])
def Newton(x0,tol,Nmax):
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    for its in range(Nmax):
        J = evalJ(x0)
        Jinv = inv(J)
        F = evalF(x0[0], x0[1])
        G = evalG(x0[0], x0[1])
        x1 = x0 - Jinv.dot([F,G])
        print(x1)
        if (norm(x1-x0) < tol):
            xstar = x1
            ier =0
            return[xstar, ier, its]
        x0 = x1
    xstar = x1
    ier = 1
    return[xstar,ier,its]

# Initial guess, tolerance, and maximum iterations
initial_guess = np.array([1,1])
tolerance = 1e-6
max_iterations = 100

# Call Newton's method to find the root
root, error_code, iterations = Newton(initial_guess, tolerance, max_iterations)

if error_code == 0:
    print(f"Approximate root: {root}")
    print(f"Iterations: {iterations}")
else:
    print("Method did not converge within the maximum number of iterations.")