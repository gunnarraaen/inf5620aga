import os, sys
import numpy as np
from mayavi import mlab as ml
#from mayavi import *
#from mayavi.api import OffScreenEngine
#e = OffScreenEngine()
#e.start()
movname = sys.argv[2]
frames = int(sys.argv[1])
dt = 0.003;
Lx = 4;
Nx = 200

x = np.linspace(0,Lx,Nx);
y = np.linspace(0,Lx,Nx);
X,Y = np.meshgrid(x,y);
#ml.options.offscreen = True
size = 500,500
fig = ml.figure(size= size);
#fig.scene.anti_aliasing_frames=0
#fig = gcf()
extent = [0,Nx-1,0,Nx-1,-7,7]
for i in range(frames):
	ml.clf(figure=fig)
	u = np.loadtxt('test.d%07d'%i);
	fname = '_tmp%07d.png' % i
	ml.surf(x,y,u,extent=extent,figure=fig,vmax=0.7,vmin=-0.7)
#	ml.axes(figure=fig,xlabel='x',ylabel='y',zlabel='z',extent=extent)
#	print  ml.view()
	ml.view(45,54,547,(100,100,0), figure=fig)
#	fig.scene.parallel_projection = True
#	fig.scene.camera.parallel_scale = 5;
	ml.title('t = %s' % (i*dt),figure=fig)
	ml.savefig(filename=fname,figure=fig)

#	files.append(fname)
	print 'Saving frame', fname

fig.scene.disable_render = False
os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=20 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s.mpg" % movname);