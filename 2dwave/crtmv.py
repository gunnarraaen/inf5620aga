import os, sys
#from matplotlib import *
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
Nx = 1000
i = 0; j = 0;
u = array([i/1000.0 for i in xrange(Nx*Nx)])
u.shape = Nx,Nx

files = []
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = linspace(0,1,Nx);
y = linspace(0,1,Nx);
X,Y = meshgrid(x,y);
frames = 100
for i in range(frames):
	ax.cla()
	lines = open('test.d%07d'%i,'r').readlines();
	for k in xrange(0,Nx):
		numbers = lines[k].split(" ")
		for l in xrange(0,Nx):
			u[k][l] = double(numbers[l]);
	fname = '_tmp%07d.png'%i
#	surf(X,Y,u,axis=[0,1,0,1,-3,3], savefig=fname);
	ax.plot_surface(X,Y,u, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=True);
	ax.set_zlim3d(-3,3);
	ax.set_axis_on();
	print 'Saving frame', fname
	fig.savefig(fname)
	files.append(fname)

os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=3 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o antimation.mpg");
