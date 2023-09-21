# import libraries
import numpy as np
def driver():
    # test functions
    f1 = lambda x: (10/(x+4))**(0.5)
    # fixed point is alpha1 = 1.4987....
    # fixed point is alpha2 = 3.09...
    Nmax = 100
    tol = 1e-6
    # test f1 '''
    x0 = 0.0
    myList1 = fixedpt(f1,x0,tol,Nmax)
    print("Fixed Point Method: ",myList1)
    #test f2 '''
    x0 = 0.0
    myList2 = vectorApprox(myList1,tol,Nmax)
    print("Aitken's Method: ",myList2)
# define routines
def fixedpt(f,x0,tol,Nmax):
    ''' x0 = initial guess'''
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''
    count = 0
    iterationList = np.array([x0])
    while (count < Nmax):
        count = count +1
        x1 = f(x0)
        iterationList = np.append(iterationList, x1)
        if (abs(x1-x0) <tol):
            return iterationList
        x0 = x1
    iterationList = ['Error']
    return iterationList

def vectorApprox(sequence, tol, nmax):
    newVector = np.array([sequence[0]])
    count = 0
    while(count < nmax):
        x1 = sequence[count] - ((sequence[count+1] - sequence[count])**2)/(sequence[count+2] - 2*sequence[count+1] + sequence[count])
        count = count + 1
        newVector = np.append(newVector,x1)
        if (abs(x1 - newVector[len(newVector) - 2]) < tol):
            return newVector
    newVector = ['Error']
    return newVector

driver()



