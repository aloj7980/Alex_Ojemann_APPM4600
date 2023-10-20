import numpy as np
import matplotlib.pyplot as plt

#Problem 1
def vandermonde_interpolation(x, y):
    V = np.vander(x, increasing=True)
    c = np.linalg.solve(V, y)
    return c

def eval_coef(xval, coefs, N):
    sum = 0
    for i in range(len(coefs)):
        sum = sum + coefs[i]*xval**i
    return sum

# Example usage:
# Define your data points
f = lambda x: (1/(1+(10*x)**2))

N = 18
a = -1
b = 1

xint = np.linspace(a, b, N + 1)

yint = f(xint)

# Perform interpolation
coefficients = vandermonde_interpolation(xint, yint)
print("Coefficients of the interpolating polynomial:", coefficients)

Neval = 1000
xeval = np.linspace(a, b, Neval + 1)
yeval = np.zeros(Neval + 1)

for kk in range(Neval + 1):
    yeval[kk] = eval_coef(xeval[kk], coefficients,  N)

fex = f(xeval)

plt.figure()
plt.plot(xeval, fex, 'o')
plt.plot(xeval, yeval, 'o')
plt.legend()
plt.show()

#Problem 2

def eval_lagrange(xeval, xint, yint, N):
    lj = np.ones(N + 1)

    for count in range(N + 1):
        for jj in range(N + 1):
            if (jj != count):
                lj[count] = lj[count] * (xeval - xint[jj]) / (xint[count] - xint[jj])

    yeval = 0.

    for jj in range(N + 1):
        yeval = yeval + yint[jj] * lj[jj]

    return (yeval)

f = lambda x: (1/(1+(10*x)**2))

N = 20
''' interval'''
a = -1
b = 1

''' create equispaced interpolation nodes'''
xint = np.linspace(a, b, N + 1)

''' create interpolation data'''
yint = f(xint)


Neval = 1000
xeval = np.linspace(a, b, Neval + 1)
yeval = np.zeros(Neval + 1)

''' evaluate polynomial '''
for kk in range(Neval + 1):
    yeval[kk] = eval_lagrange(xeval[kk], xint, yint, N)

''' create vector with exact values'''
fex = f(xeval)

plt.figure()
plt.plot(xeval, fex, 'o')
plt.plot(xeval, yeval, 'o')
plt.legend()
plt.show()

#Problem 3

def get_chebyshev_points(N):
    j_values = np.arange(1, N+2)
    x_values = np.cos((2 * j_values - 1) * np.pi / (2 * N))
    return x_values

''' create chebyshev interpolation nodes'''
xint = get_chebyshev_points(N)

''' create interpolation data'''
yint = f(xint)


Neval = 1000
xeval = np.linspace(a, b, Neval + 1)
yeval = np.zeros(Neval + 1)


''' evaluate polynomial using Chebyshev points'''
coefficients = vandermonde_interpolation(xint, yint)
for kk in range(Neval + 1):
    yeval[kk] = eval_coef(xeval[kk], coefficients,  N)

''' create vector with exact values'''
fex = f(xeval)

plt.figure()
plt.plot(xeval, fex, 'o')
plt.plot(xeval, yeval, 'o')
plt.legend()
plt.show()


