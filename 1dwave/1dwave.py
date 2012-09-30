import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # For plotting graphs.
import numpy as np
import subprocess                 # For issuing commands to the OS.
import os
import sys                        # For determining the Python version.
import matplotlib.patches as mpatches

t0 = 0.0
T = 200.0
dt = 0.1
dx = 0.1
Nt = int((T-t0)/dt)
x0 = -20.
x1 = 20.
Nx = int((x1-x0)/dx)
t = np.linspace(t0,T,Nt)


xplot = np.linspace(x0,x1,Nx)
x = np.linspace(x0,x1,Nx)
xp = np.zeros(Nx)
xpp = np.zeros(Nx)

#initial positions
sigma = 2.0
xinit = np.exp(-x*x/(sigma*sigma))

def factor(x,t):
    factor = np.logical_and((x-x0)/(x1-x0) >= 2./3.,((x-x0)/(x1-x0) <= 4./5.))*0.4 + \
        np.logical_or((x-x0)/(x1-x0) <= 2./3., (x-x0)/(x1-x0)>= 4./5.)*0.3
    return factor*factor*dt*dt/(dx*dx)

#first step
for i in range(1,Nx-1):
	x[i] =  xinit[i] + (factor(xplot[i], t[0])/2.)*(xinit[i+1] - 2*xinit[i] + xinit[i-1])
x[0] = 0.
x[-1] = 0.
xpp = np.copy(xinit)
xp = np.copy(x)

"""
plt.plot(xplot,dxdtinit)
plt.hold('on')
plt.show()
"""

# Some rectangles for visualization
# first from 0 to 2/3: 
xy = np.array([x0,-1.3])
height = 2.6*(x1-x0)
width = (2./3.)*(x1-x0)
q1 = mpatches.Rectangle(xy, width, height, facecolor='yellow', edgecolor='yellow')

# then from 2/3 to 4/5:
xy = np.array([x0+2./3.*(x1-x0),-1.3])
height = 2.6*(x1-x0)
width = (4./5.-2./3.)*(x1-x0)
q2 = mpatches.Rectangle(xy, width, height, facecolor='cyan', edgecolor='cyan')

#finally from 4/5 to 1
xy = np.array([x0+4./5.*(x1-x0),-1.3])
height = 2.6*(x1-x0)
width = (1.-4./5.)*(x1-x0)
q3 = mpatches.Rectangle(xy, width, height, facecolor='yellow', edgecolor='yellow')


for n in range(Nt):
#    for i in range(1,Nx-1):
#        x[i] = 2*xp[i] - xpp[i] + factor*(xp[i-1] - 2*xp[i] + xp[i+1])
    x[1:-1] = (2*xp[1:-1] - xpp[1:-1] + factor(xplot[1:-1],t[n])*(xp[:-2] - 2*xp[1:-1] + xp[2:]))
    x[0] = x[1] #(2*(1-factor(xplot[0], t[n]) ))*xp[0] - xpp[0] + factor(xplot[0], t[n])*xp[1];
    x[-1] = x[-2] #(2*(1-factor(xplot[-1],t[n]) ))*xp[-1] - xpp[-1] + factor(xplot[0],t[n])*xp[-2];
    xpp = np.copy(xp)
    xp = np.copy(x)
    if np.mod(n,10) == 0:
	    #
	    # The next four lines are just like MATLAB.
	    #
            plt.gca().add_patch(q1)
            plt.gca().add_patch(q2)
            plt.gca().add_patch(q3)
	    plt.plot(xplot,x,'b',hold=True)
	    plt.axis((x0,x1,-1.3, 1.3))
	    plt.xlabel('x-position')
	    plt.ylabel('y-position')

	    #
	    # Notice the use of LaTeX-like markup.
	    #
	    plt.title(r'$t=%06d$' % n, fontsize=20)

	    #
	    # The file name indicates how the image will be saved and the
	    # order it will appear in the movie.  If you actually wanted each
	    # graph to be displayed on the screen, you would include commands
	    # such as show() and draw() here.  See the matplotlib
	    # documentation for details.  In this case, we are saving the
	    # images directly to a file without displaying them.
	    #
	    filename = str('%06d' % n) + '.png'
	    plt.savefig(filename, dpi=100)

	    #
	    # Let the user know what's happening.
	    #
	    print 'Wrote file', filename

	    #
	    # Clear the figure to make way for the next image.
	    #
	    plt.clf()

#
# Now that we have graphed images of the dataset, we will stitch them
# together using Mencoder to create a movie.  Each image will become
# a single frame in the movie.
#
# We want to use Python to make what would normally be a command line
# call to Mencoder.  Specifically, the command line call we want to
# emulate is (without the initial '#'):
# mencoder mf://*.png -mf type=png:w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o output.avi
# See the MPlayer and Mencoder documentation for details.
#

command = ('mencoder',
           'mf://*.png',
           '-mf',
           'type=png:w=800:h=600:fps=25',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           'output.avi')

#os.spawnvp(os.P_WAIT, 'mencoder', command)

print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
subprocess.check_call(command)

print "\n\n The movie was written to 'output.avi'"

print "\n\n You may want to delete *.png now.\n\n"
