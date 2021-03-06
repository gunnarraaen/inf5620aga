import os, sys
import numpy as np
#from enthought.mayavi import mlab as ml
from mayavi import  mlab as ml
import pylab as pl

"""
movname = sys.argv[2]
frames = int(sys.argv[1])
dt = 0.001;
Lx = 4;
Nx = 400

x = np.linspace(0,Lx,Nx);
y = np.linspace(0,Lx,Nx);
X,Y = np.meshgrid(x,y);
size = 500,500
"""

def mcrtmv(frames, mvname, dt,Lx,Ly,Nx,Ny):
	x = np.linspace(0,Lx,Nx);
	y = np.linspace(0,Ly,Ny);
	X,Y = np.meshgrid(x,y);
	size = 500,500
	
	fig = ml.figure(size= size, bgcolor=(1.,1.,1.));

	#fig.scene.anti_aliasing_frames=07

	#extent = [0,Nx-1,0,Ny-1,-30,30]
	
	ml.clf(figure=fig)
	u = np.loadtxt('test.d%07d'%1);
	fname = '_tmp%07d.png' % 1
	s = ml.surf(x,y,u,figure=fig,vmin=-0,vmax=1)
	ml.axes(extent=[0,Lx,0,Ly,-2,2])
	ml.colorbar()
	ml.xlabel('x position')
	ml.ylabel('y position')
	ml.zlabel('wave amplitude')
	"""
	pl.ion()
	arr = ml.screenshot()
	img = pl.imshow(arr)
	pl.axis('off')
	"""
	for i in range(2,frames):
		#arr = ml.screenshot()
		#img.set_array(arr)
		u = np.loadtxt('test.d%07d'%i);
		s.mlab_source.scalars = u
		fname = '_tmp%07d.png' % i
		#pl.savefig(filename=fname)#,figure=fig)
		#print 'Saving frame', fname
		#pl.draw()

	#os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=20 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s.mpg" % mvname);
if __name__ == "__main__":
	mcrtmv(int(sys.argv[1]),sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),int(sys.argv[6]),int(sys.argv[7]));
