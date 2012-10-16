from numpy import *
# Create data
x = linspace(0,10,40)
y = copy(x)
# z = zeros((len(x), len(y)))
# for xi in arange(len(x)):
#     for yi in arange(len(y)):
#         z[xi, yi] = 1.3 * sin(x[xi]) + sin(y[yi])
# Plot data with mayavi
z = loadtxt('solution_000038.txt')
from mayavi import mlab
mlab.surf(x, y, z)
mlab.axes()
# With the following command, a GUI tool to edit plot parameters is opened:
mlab.show()