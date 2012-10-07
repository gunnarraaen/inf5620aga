import os, sys
import numpy as np
#from enthought.mayavi import mlab as ml
from mayavi import  mlab as ml

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
	y = np.linspace(0,Lx,Nx);
	X,Y = np.meshgrid(x,y);
	size = 500,500
	
	fig = ml.figure(size= size);

	fig.scene.anti_aliasing_frames=07

	extent = [0,Nx-1,0,Ny-1,-30,30]
	for i in range(1,frames):
		ml.clf(figure=fig)
		u = np.loadtxt('solution_%06d.txt'%i);
		u[0][0] = 7
		u[0][1] = -7
		fname = '_tmp%07d.png' % i
		ml.surf(x,y,u,extent=extent,figure=fig,vmax=1,vmin=-0)
	#	ml.axes(figure=fig,xlabel='x',ylabel='y',zlabel='z',extent=extent)
	#	print  ml.view()
	#	ml.view(45,54,547,(100,100,0), figure=fig)
	#	fig.scene.parallel_projection = True
	#	fig.scene.camera.parallel_scale = 5;
		ml.title('t = %s' % (i*dt),figure=fig)
		ml.savefig(filename=fname,figure=fig)

		print 'Saving frame', fname

	fig.scene.disable_render = False
	os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=20 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s.mpg" % mvname);
