import numpy as np
import math
import time
from numpy.linalg import inv
from numpy.linalg import norm

def evalF(x):
    F = np.zeros(3)
    F[0] = 3*x[0]-math.cos(x[1]*x[2])-1/2
    F[1] = x[0]-81*(x[1]+0.1)**2+math.sin(x[2])+1.06
    F[2] = np.exp(-x[0]*x[1])+20*x[2]+(10*math.pi-3)/3
    return F

def evalJ(x):
    J = np.array([[3.0, x[2]*math.sin(x[1]*x[2]), x[1]*math.sin(x[1]*x[2])],
        [2.*x[0], -162.*(x[1]+0.1), math.cos(x[2])],
        [-x[1]*np.exp(-x[0]*x[1]), -x[0]*np.exp(-x[0]*x[1]), 20]])
    return J

def LazyNewton(x0,tol,Nmax):
    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    J = evalJ(x0)
    Jinv = inv(J)
    for its in range(Nmax):
        F = evalF(x0)
        x1 = x0 - Jinv.dot(F)
        if (norm(x1-x0) < tol):
            xstar = x1
            ier =0
            return[xstar, ier,its]
        x0 = x1
    xstar = x1
    ier = 1
    return[xstar,ier,its]

def SlackerNewton(x0,tol,Nmax):
    ''' Slacker Newton = use only the inverse of the Jacobian when '''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    J = evalJ(x0)
    Jinv = inv(J)
    F = [0,0,0]
    jcount = 0
    for its in range(Nmax):
        prevF = F
        F = evalF(x0)
        if(norm(prevF - F) > 0.1):
            J = evalJ(x0)
            Jinv = inv(J)
            jcount = jcount + 1
        x1 = x0 - Jinv.dot(F)
        if (norm(x1-x0) < tol):
            xstar = x1
            ier =0
            return[xstar, ier,its,jcount]
        x0 = x1
    xstar = x1
    ier = 1
    return[xstar,ier,its,jcount]

def driver():
    x0 = np.array([0.1, 0.1, -0.1])
    Nmax = 100
    tol = 1e-10
    t = time.time()
    for j in range(20):
        [xstar,ier,its] = LazyNewton(x0,tol,Nmax)
    elapsed = time.time()-t
    print(xstar)
    print('Lazy Newton: the error message reads:',ier)
    print('Lazy Newton: took this many seconds:',elapsed/20)
    print('Lazy Newton: number of iterations is:',its)
    print('Lazy Newton: number of inversions is: 0')
    t = time.time()
    for j in range(20):
        [xstar,ier,its,jcount] = SlackerNewton(x0, tol,Nmax)
    elapsed = time.time()-t
    print(xstar)
    print('Slacker Newton: the error message reads:',ier)
    print('Slacker Newton: took this many seconds:',elapsed/20)
    print('Slacker Newton: number of iterations is:',its)
    print('Slacker Newton: number of inversions is:', jcount)

driver()

def Newton(x0,tol,Nmax):
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    for its in range(Nmax):
        J = evalJ(x0)
        Jinv = inv(J)
        F = evalF(x0)
        x1 = x0 - Jinv.dot(F)
        if (norm(x1-x0) < tol):
            xstar = x1
            ier =0
            return[xstar, ier, its]
        x0 = x1
    xstar = x1
    ier = 1
    return[xstar,ier,its]

def FNewton(x0,tol,Nmax):
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    for its in range(Nmax):
        J = evalJ(x0)
        Jinv = fdiff(J)
        F = evalF(x0)
        x1 = x0 - Jinv.dot(F)
        if (norm(x1-x0) < tol):
            xstar = x1
            ier =0
            return[xstar, ier, its]
        x0 = x1
    xstar = x1
    ier = 1
    return[xstar,ier,its]

def fdiff(J,h):
    #Return jacobian approximation using forward differencing
    for(i in range(len(J))):
        for(j in range(len(j[i]))):
            J[i,j] = (F(x1+h,x2+h)-F(x1,x2))/h