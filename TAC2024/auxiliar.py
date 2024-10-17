import numpy as np

xs = 10.00
ys = 20.00
xi = -10.00
yi = -20.00

alpha = -67.521 * np.pi / 180.0

etas = ys * np.cos(alpha) - xs * np.sin(alpha) 
etai = yi * np.cos(alpha) - xi * np.sin(alpha)
print('etas =',etas)
print('etai =', etai)
h = np.abs(etas - etai)
print('h =', h)
epss = -2.60e-3
epsi = -0.52826e-3
x = epss/(epss-epsi)*h
print('x =', x)

