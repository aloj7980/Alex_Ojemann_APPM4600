
#Part a

import numpy as np

t = np.arange(0, np.pi + np.pi/30, np.pi/30)

y = np.cos(t)

S = np.sum(t * y)

print(f"The sum is: {S}")

#Part b

import matplotlib.pyplot as plt
import random

R = 1.2
delta_r = 0.1
f = 15
p = 0

theta = np.linspace(0, 2 * np.pi, 1000)

x = R * (1 + delta_r * np.sin(f * theta + p)) * np.cos(theta)
y = R * (1 + delta_r * np.sin(f * theta + p)) * np.sin(theta)

plt.figure(1)
plt.plot(x, y)
plt.title("Parametric Curve for R=1.2, δr=0.1, f=15, p=0")
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True)

plt.figure(2)
for i in range(1, 11):
    R_i = i
    delta_r_i = 0.05
    f_i = 2 + i
    p_i = random.uniform(0, 2)

    x_i = R_i * (1 + delta_r_i * np.sin(f_i * theta + p_i)) * np.cos(theta)
    y_i = R_i * (1 + delta_r_i * np.sin(f_i * theta + p_i)) * np.sin(theta)

    plt.plot(x_i, y_i, label=f"Curve {i}")

plt.title("Parametric Curves for 10 Sets of Parameters")
plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.grid(True)

plt.show()

