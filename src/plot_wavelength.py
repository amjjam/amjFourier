# Plot polynomial function provided by Paolo, and also plots values from
# amjFourier::wavelength, to validate that the function produces the same

import numpy as np
import matplotlib.pyplot as plt
f=open("wavelength.dat")
d=np.fromfile(f,dtype=np.dtype([('i',np.int32),('w',np.float64)]))

# Load data amjFourier::wavelength output from wavelength.dat


# Polynomial coefficients (highest power first)

#
# 9 October e-mail from Paolo:
#coeffs =  np.array([6.95694016e-06,-5.20602460e-03, 1.46461339e+00,
#                    -1.96798243e+02,1.28465801e+04])

#
# 19 December e-mail from Paolo:
coeffs = np.array([-3.63299018e-06,3.46705393e-03,-1.17335426e+00,
                  1.56363989e+02,-4.71516112e+03])

# x range
x = np.linspace(140,256, 400)

# Evaluate polynomial
y = np.polyval(coeffs, x)

# Plot
plt.figure()
plt.plot(x, y)
plt.plot(d['i'],d['w'])
plt.xlabel("x")
plt.ylabel("p(x)")
plt.title("Polynomial plot")
plt.grid(True)
plt.show()
