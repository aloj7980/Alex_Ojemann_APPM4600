import numpy as np

#3.1
x = np.linspace(0, 9, 10)
y = np.arange(10)
print(x)
print(y)

#3.2 and 3.3
print('The first three entries of x are', x[:3])

#3.4 and 3.5
w = 10**(-np.linspace(1,10,10))
s = 3 * w
x = np.linspace(1, len(w), len(w))
import matplotlib.pyplot as plt
plt.semilogy(x, w)
plt.semilogy(x, s)
plt.show()
