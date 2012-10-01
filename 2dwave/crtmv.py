
import matplotlib
matplotlib.use('Agg');
import os, sys

from matplotlib.ticker import LinearLocator, FormatStrFormatter

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from numpy import *

movname = sys.argv[2]
frames = int(sys.argv[1])
dt = 0.003;
Lx = 4;
Nx = 100
#i = 0; j = 0;
#u = array([i/1000.0 for i in xrange(Nx*Nx)])
#u.shape = Nx,Nx

files = []
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = linspace(0,Lx,Nx);
y = linspace(0,Lx,Nx);
X,Y = meshgrid(x,y);

for i in range(frames):
	ax.cla()
#	lines = open('test.d%07d'%i,'r').readlines();
#	for k in xrange(0,Nx):
#		numbers = lines[k].split(" ")
#		for l in xrange(0,Nx):
#			u[k][l] = double(numbers[l]);

	u = loadtxt('test.d%07d'%i);
	fname = '_tmp%07d.png' % i
#	surf(X,Y,u,axis=[0,1,0,1,-3,3], savefig=fname);
	ax.plot_surface(X,Y,u, rstride=2, cstride=2, cmap=cm.jet, linewidth=0, antialiased=True);
#	ax.plot_surface(X,Y,u,cmap=cm.jet,antialiased=True);
	ax.set_zlim3d(-2,2);
	ax.set_axis_on();
	ax.set_title('t = %s' %(i*dt));
	print 'Saving frame', fname
	fig.savefig(fname)
	files.append(fname)

os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=20 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s.mpg" % movname);
