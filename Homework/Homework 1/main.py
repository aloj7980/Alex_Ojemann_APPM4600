# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

def p(x:float):
    return x**9 - 18*x**8 + 144*x**7 - 672*x**6 + 2016*x**5 - 4032*x**4 + 5376*x**3 - 4608*x**2 + 2304*x - 512

def q(x:float):
    return (x-2)**9

import numpy as np
x = np.linspace(1.920 , 2.080 , 161)
print(x)
import matplotlib.pyplot as plt
plt.plot(x,p(x))
plt.plot(x,q(x))
plt.legend(['(x-2)^9 expanded', '(x-2)^9 not expanded'])
plt.show()
