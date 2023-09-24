import numpy as np
def bisection(f,a,b,tol):
    # Inputs:
    # f,a,b - function and endpoints of initial interval
    # tol - bisection stops when interval length < tol
    # Returns:
    # astar - approximation of root
    # ier - error message
    # - ier = 1 => Failed
    # - ier = 0 == success
    # first verify there is a root we can find in the interval
    fa = f(a)
    fb = f(b)
    if (f(a)*f(b)>0):
        ier = 1
        astar = a
        return [astar, ier]
    # verify end points are not a root
    if (f(a) == 0):
        astar = a
        ier = 0
        return [astar, ier]
    if (f(b) == 0):
        astar = b
        ier = 0
        return [astar, ier]
    count = 0
    d = 0.5 * (a + b)
    while (abs(d - a) > tol):
        fd = f(d)
        if (fd == 0):
            astar = d
            ier = 0
            print(count)
            return [astar, ier]
        if (fa * fd < 0):
            b = d
        else:
            a = d
            fa = fd
        d = 0.5 * (a + b)
        count = count + 1
        #print('abs(d-a) = ', abs(d-a))
    astar = d
    ier = 0
    print(count)
    return [astar, ier]


print(bisection(lambda x: 2*x - 1 - np.sin(x),0,3.14,10**-8))
print(bisection(lambda x: (x-5)**-9,4.82,5.2,10**-4))
print(bisection(lambda x: x**9 - 45*x**8 +945*x**7 -12600*x**6+117649*x**5-823543*x**4+4477458*x**3-17887890*x**2+50328438*x-97656250,4.82,5.2,10**-4))
print(bisection(lambda x: x**3 + x - 4,1,4,10**-3))

import matplotlib.pyplot as plt

def f(x):
    return x - 4 * np.sin(2 * x) - 3

x = np.linspace(-2 * np.pi, 4 * np.pi, 400)  # Adjust the range as needed

y = f(x)

plt.figure(figsize=(8, 6))
plt.plot(x, y, label='f(x) = x - 4sin(2x) - 3', color='b')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Plot of f(x) = x - 4sin(2x) - 3')
plt.grid(True)
plt.legend()
plt.axhline(0, color='r', linestyle='--', linewidth=0.7)  # Horizontal line at y=0
plt.axvline(0, color='r', linestyle='--', linewidth=0.7)  # Vertical line at x=0
plt.show()


def fixedpt(f,x0,tol,Nmax):
    ''' x0 = initial guess'''
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''
    count = 0
    while (count <Nmax):
        count = count +1
        x1 = f(x0)
        if (abs(x1-x0) <tol):
            xstar = x1
            ier = 0
            return [xstar,ier]
        x0 = x1
    xstar = x1
    ier = 1
    return [xstar, ier]

print(fixedpt(lambda x: x - 4 * np.sin(2 * x) - 3,1,0.0000001, 500))

