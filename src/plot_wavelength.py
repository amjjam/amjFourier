import numpy as np
import matplotlib.pyplot as plt

# Polynomial coefficients (highest power first)
coeffs = np.array([6.95694016e-06,-5.20602460e-03, 1.46461339e+00,-1.96798243e+02,
       1.28465801e+04])   # 2x^2 - 3x + 1

# x range
x = np.linspace(140,256, 400)

# Evaluate polynomial
y = np.polyval(coeffs, x)

# Plot
plt.figure()
plt.plot(x, y)
plt.xlabel("x")
plt.ylabel("p(x)")
plt.title("Polynomial plot")
plt.grid(True)
plt.show()
