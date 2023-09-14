# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

f = lambda x: x**2*(x-1)


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
    return [astar, ier]
#Success
print(bisection(f,0.5,2, 0.000001))
#This is unsuccessful because f of the start point and f of the end point are both negative
print(bisection(f,-1,0.5, 0.000001))
#Success
print(bisection(f,-1,2, 0.000001))


#Part 3
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

print(fixedpt(lambda x:x*(1+((7-x**5)/x**2))**3,1,0.0000001, 500))
#It breaks because the derivative is too large