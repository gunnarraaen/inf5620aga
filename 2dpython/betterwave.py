"""
2D wave equation solved by finite differences::

  dt, cpu_time, T = solver(I, V, f, c, Lx, Ly, Nx, Ny, dt, T,
                        user_action=None, version='scalar',
                        dt_safety_factor=1)

Solve the 2D wave equation u_tt + u_t = u_xx + u_yy + f(x,t) on (0,L) with
du/dn = 0 at the boundary

Nx and Ny are the total number of mesh cells in the x and y
directions. The mesh points are numbered as (0,0), (1,0), (2,0),
..., (Nx,0), (0,1), (1,1), ..., (Nx, Ny).

dt is the time step. If dt<=0, an optimal time step is used.
T is the stop time for the simulation.

I, V, f are functions: I(x,y), V(x,y), f(x,y,t). V and f
can be specified as None or 0, resulting in V=0 and f=0.

user_action: function of (u, x, y, t, n) called at each time
level (x and y are one-dimensional coordinate vectors).
This function allows the calling code to plot the solution,
compute errors, etc.
"""
import time
from scitools.std import *
import numpy as np
import os

def solver(I, V, f, q, beta, Lx, Ly, Nx, Ny, dt, T,
           user_action=None, version='scalar',
           dt_safety_factor=1):

    if version == 'scalar':
        advance = advance_scalar
    elif version == 'vectorized':
        advance = advance_vectorized
    elif version == 'cpp':
    	advance = advance_cpp

    	    
    import time                # measure CPU time
    t0 = time.clock()

    x = linspace(0, Lx, Nx+1)  # mesh points in x dir
    y = linspace(0, Ly, Ny+1)  # mesh points in y dir
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    u_2  = zeros((Nx+1,Ny+1))
    u_1  = zeros((Nx+1,Ny+1))
    u  = zeros((Nx+1,Ny+1))
    xv = x[:,newaxis]          # for vectorized function evaluations
    yv = y[newaxis,:]

    stability_limit = (1./sqrt(np.max(abs(q(xv,yv))))*(1/sqrt(1/dx**2 + 1/dy**2)))

    if dt <= 0:                 # max time step?
        dt = dt_safety_factor*stability_limit
    elif dt > stability_limit:
        print 'error: dt=%g exceeds the stability limit %g' % \
              (dt, stability_limit)

    N = int(round(T/float(dt)))
    t = linspace(0, N*dt, N+1)    # mesh points in time
    T = N*dt
    #Cx2 = (c*dt/dx)**2;  Cy2 = (c*dt/dy)**2    # help variables
    dt2 = dt**2

    # Allow f and V to be None or 0
    if f is None or f == 0:
        f = (lambda x, y, t: 0) if version == 'scalar' else \
            lambda x, y, t: zeros((x.shape[0], y.shape[1]))
        # or simpler: x*y*0
    if V is None or V == 0:
        V = (lambda x, y: 0) if version == 'scalar' else \
            lambda x, y: zeros((x.shape[0], y.shape[1]))

    # Set initial condition

    if version == 'scalar':
        for i in range(0, Nx+1):
            for j in range(0, Ny+1):
                u_1[i,j] = I(x[i], y[j])
    else: # use vectorized version
        u_1[:,:] = I(xv, yv)

    if user_action is not None:
        user_action(u_1, x, xv, y, yv, t, 0)

    # Special formula for first time step
    n = 0
    if version == 'scalar':
        for i in range(0, Nx+1):
            for j in range(0, Ny+1):
                ileft = i-1
                iright = i+1
                jleft = j-1
                jright = j+1
                if i-1 < 0:
                    ileft = iright
                if i+1 > Nx:
                    iright = ileft
                if j-1 < 0:
                    jleft = jright
                if j+1 > Ny:
                    jright = jleft
                Cx2l = q(x[i]-dx/2., y[j])*(dt/dx)**2;  Cx2r = q(x[i]+dx/2.,y[j])*(dt/dx)**2    # help variables
                Cy2l = q(x[i], y[j]-dy/2.)*(dt/dy)**2;  Cy2r = q(x[i],y[j]+dy/2.)*(dt/dy)**2    # help variables
                u[i,j] = u_1[i,j] + (1-beta*dt)*dt*V(x[i], y[j]) + \
                0.5*(-Cx2l*(u_1[i,j] - u_1[ileft,j]) + Cx2r*(u_1[iright,j] - u_1[i,j])) + \
                0.5*(-Cy2l*(u_1[i,j] - u_1[i,jleft]) + Cy2r*(u_1[i,jright] - u_1[i,j])) + \
                0.5*dt2*f(x[i], y[j], t[n])
        
        
        
    else:  # use vectorized version
        f_a = f(xv, yv, t[n])               # precompute, size as u
        qv = q(xv+dx/2., yv+dy/2.)               # precompute, size as u
        V_a = V(xv, yv)
        Cx2 = qv[:,:]*(dt/dx)**2
        Cy2 = qv[:,:]*(dt/dy)**2

        Nx = u.shape[0];  Ny = u.shape[1]
        xfactor = dt*dt/(dx*dx)
        yfactor = dt*dt/(dy*dy)
        #u_ghost = np.zero(Nx+1,Ny+1)
        u_ghost_1 = np.zeros((Nx+2,Ny+2))
        u_ghost_1[1:-1,1:-1] = u_1[:,:] 
        u_ghost_1[0,1:-1] = u_1[1,:]; u_ghost_1[-1,1:-1] = u_1[-2,:];
        u_ghost_1[1:-1,0] = u_1[:,1]; u_ghost_1[1:-1,-1] = u_1[:,-2];
        u_ghost_1[0,0] = u_1[1,1]; u_ghost_1[0,-1] = u_1[1,-2];
        u_ghost_1[-1,0] = u_1[-2,1]; u_ghost_1[-1,-1] = u_1[-2,-2]
        
        u[:,:] = (u_ghost_1[1:-1,1:-1] + (1-beta*dt)*dt*V_a[:,:] + \
            0.5*(xfactor*(qv[:,:]*(u_ghost_1[2:,1:-1] - u_ghost_1[1:-1,1:-1]) - \
            qv[:,:]*(u_ghost_1[1:-1,1:-1] - u_ghost_1[0:-2,1:-1]))) + \
            0.5*(yfactor*(qv[:,:]*(u_ghost_1[1:-1,2:] - u_ghost_1[1:-1,1:-1]) - \
            qv[:,:]*(u_ghost_1[1:-1,1:-1] - u_ghost_1[1:-1,0:-2]))) + \
            0.5*dt**2*f_a[:,:])
        
    if user_action is not None:
        user_action(u, x, xv, y, yv, t, 1)

    u_2[:,:] = u_1; u_1[:,:] = u

    if version == 'scalar':
        # use f(x,y,t) function
        for n in range(1,N+1):
        	u = advance(u, u_1, u_2, f, beta, x, y, t, n, q, dx, dy, dt)
        	filename = 'solution_%06d.txt' % n
        	np.savetxt(filename, u, fmt='%.18e', delimiter=' ', newline='\n')
        	u_2[:,:], u_1[:,:] = u_1, u
    elif version == 'vectorized':
        for n in range(1,N+1):
            f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
            qv[:,:] = q(xv, yv)         # precompute, size as u
            u = advance(u, u_1, u_2, f_a, beta, x, y, t, n, qv, dx, dy, dt)
            filename = 'solution_%06d.txt' % n
            np.savetxt(filename, u, fmt='%.18e', delimiter=' ', newline='\n')
            u_2[:,:], u_1[:,:] = u_1, u
        #np.savetxt(filename, u, fmt='%.18e', delimiter=' ', newline='\n')

    elif version == 'cpp':
    	file1 = 'uone.txt'
    	np.savetxt(file1, u_1, fmt='%.18e', delimiter=' ', newline='\n')
    	file2 = 'utwo.txt'
    	np.savetxt(file1, u_2, fmt='%.18e', delimiter=' ', newline='\n')
    	os.system('wave2d'+' '+file1+' '+file2)
        file3 = 'parameters.txt'
    t1 = time.clock()
    # dt might be computed in this function so return the value
    return dt, t1 - t0, T

def advance_scalar(u, u_1, u_2, f, beta, x, y, t, n, q, dx, dy, dt):
    Nx = u.shape[0]-1;  Ny = u.shape[1]-1
    for i in range(0, Nx+1):
        for j in range(0, Ny+1):
            Cx2l = q(x[i]-dx/2., y[j])*(dt/dx)**2;  Cx2r = q(x[i]+dx/2.,y[j])*(dt/dx)**2    # help variables
            Cy2l = q(x[i], y[j]-dy/2.)*(dt/dy)**2;  Cy2r = q(x[i],y[j]+dy/2.)*(dt/dy)**2    # help variables
            ileft = i-1
            iright = i+1
            jleft = j-1
            jright = j+1
            if i-1 < 0:
                ileft = iright
            if i+1 > Nx:
                iright = ileft
            if j-1 < 0:
                jleft = jright
            if j+1 > Ny:
                jright = jleft
            u[i,j] = 1./(1+dt*beta/2.)*(2*u_1[i,j] - (1-dt*beta/2.)*u_2[i,j] + \
                (-Cx2l*(u_1[i,j] - u_1[ileft,j]) + Cx2r*(u_1[iright,j] - u_1[i,j])) + \
                (-Cy2l*(u_1[i,j] - u_1[i,jleft]) + Cy2r*(u_1[i,jright] - u_1[i,j])) + \
                dt**2*f(x[i], y[j], t[n]))
    return u

def advance_vectorized(u, u_1, u_2, f_a, beta, x, y, t, n, qv, dx, dy, dt):
    Nx = u.shape[0];  Ny = u.shape[1]
    xfactor = dt*dt/(dx*dx)
    yfactor = dt*dt/(dy*dy)
    #u_ghost = np.zero(Nx+1,Ny+1)
    u_ghost_1 = np.zeros((Nx+2,Ny+2))
    u_ghost_1[1:-1,1:-1] = u_1[:,:] 
    u_ghost_1[0,1:-1] = u_1[1,:]; u_ghost_1[-1,1:-1] = u_1[-2,:];
    u_ghost_1[1:-1,0] = u_1[:,1]; u_ghost_1[1:-1,-1] = u_1[:,-2];
    u_ghost_1[0,0] = u_1[1,1]; u_ghost_1[0,-1] = u_1[1,-2];
    u_ghost_1[-1,0] = u_1[-2,1]; u_ghost_1[-1,-1] = u_1[-2,-2]

    u[:,:] = (1./(1+dt*beta/2.))*(2*u_ghost_1[1:-1,1:-1] - (1-dt*beta/2.)*u_2[:,:] + \
       xfactor*(qv[:,:]*(u_ghost_1[2:,1:-1] - u_ghost_1[1:-1,1:-1]) - \
       qv[:,:]*(u_ghost_1[1:-1,1:-1] - u_ghost_1[0:-2,1:-1])) + \
       yfactor*(qv[:,:]*(u_ghost_1[1:-1,2:] - u_ghost_1[1:-1,1:-1]) - \
       qv[:,:]*(u_ghost_1[1:-1,1:-1] - u_ghost_1[1:-1,0:-2])) + \
       dt**2*f_a[:,:])

    return u


def run_Gaussian(plot_method=1, version='vectorized'):
    """
    Initial Gaussian bell in the middle of the domain.
    plot_method=1 makes a movie of the output with mayavi, =0 means no plot.
    """
    # Clean up plot files
    for name in glob('tmp_*.png'):
        os.remove(name)

    Lx = 10
    Ly = 10
    def c(x,y): 
        s = (np.size(x),np.size(y))
        c = np.zeros(s)
        c = c + 1
        return c

    def I(x, y):
        """Gaussian peak at (Lx/2, Ly/2)."""
        return 3*exp(-0.5*(x-Lx/2.0)**2 - 0.5*(y-Ly/2.0)**2)

    Nx = 40; Ny = 40; T = 20
    dt,time = solver(I, None, None, c, 0.0, Lx, Ly, Nx, Ny, 0, T, version=version)
    print dt
    if plot_method == 1:
    	mvname = 'test'
    	from mcrtmv import mcrtmv
    	N = int(round(T/float(dt)))
    	mcrtmv(N, dt,Lx,Ly,Nx,Ny, savemovie=False, mvname=mvname)



if __name__ == '__main__':
    import sys
    run_Gaussian()#version='scalar')
    #from scitools.misc import function_UI
    #cmd = function_UI([test_quadratic, run_efficiency_tests,
    #                   run_Gaussian, ], sys.argv)
    #eval(cmd)
