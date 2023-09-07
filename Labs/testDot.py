import numpy as np
import numpy.linalg as la
import math
def driver():
    m = 2
    n = 2
    x = np.linspace(0,np.pi,n)
    # this is a function handle. You can use it to define
    # functions instead of using a subroutine like you
    # have to in a true low level language.
    y = np.matrix([[1, 2],
        [3, 4]])
    w = np.matrix([[5, 6],
        [7, 8]])
    # evaluate the dot product of y and w
    dp = matrixMultiply(y,w,m,n)
    # print the output
    print('the matrix product is : ', dp)
    return

def dotProduct(x,y,n):
    # Computes the dot product of the n x 1 vectors x and y
    dp = 0.
    for j in range(n):
        dp = dp + x[0,j]*y[0,j]
    return dp
def matrixMultiply(x,y,m,n):
    matrix = np.array([np.zeros(n) for i in range(m)])
    for i in range(0,m):
        for j in range(0,n):
            matrix[i][j] = dotProduct(x[i,:], y[:,j].reshape(1,-1),n)
    return matrix
driver()