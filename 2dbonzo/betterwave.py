"""
2D wave equation solved by finite differences::

  dt, cpu_time = solver(I, V, f, c, Lx, Ly, Nx, Ny, dt, T,
                        user_action=None, version='scalar',
                        dt_safety_factor=1)

Solve the 2D wave equation u_tt = u_xx + u_yy + f(x,t) on (0,L) with
u=0 on the boundary and initial condition du/dt=0.

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

def solver(I, V, f, c, Lx, Ly, Nx, Ny, dt, T,
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

    xv = x[:,newaxis]          # for vectorized function evaluations
    yv = y[newaxis,:]

    stability_limit = (1/float(c))*(1/sqrt(1/dx**2 + 1/dy**2))
    if dt <= 0:                # max time step?
        dt = dt_safety_factor*stability_limit
    elif dt > stability_limit:
        print 'error: dt=%g exceeds the stability limit %g' % \
              (dt, stability_limit)
    N = int(round(T/float(dt)))
    t = linspace(0, N*dt, N+1)    # mesh points in time
    Cx2 = (c*dt/dx)**2;  Cy2 = (c*dt/dy)**2    # help variables
    dt2 = dt**2

    # Allow f and V to be None or 0
    if f is None or f == 0:
        f = (lambda x, y, t: 0) if version == 'scalar' else \
            lambda x, y, t: zeros((x.shape[0], y.shape[1]))
        # or simpler: x*y*0
    if V is None or V == 0:
        V = (lambda x, y: 0) if version == 'scalar' else \
            lambda x, y: zeros((x.shape[0], y.shape[1]))


    order = 'Fortran' if version == 'f77' else 'C'
    u   = zeros((Nx+1,Ny+1), order=order)   # solution array
    u_1 = zeros((Nx+1,Ny+1), order=order)   # solution at t-dt
    u_2 = zeros((Nx+1,Ny+1), order=order)   # solution at t-2*dt
    f_a = zeros((Nx+1,Ny+1), order=order)   # for compiled loops

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
        for i in range(1, Nx):
            for j in range(1, Ny):
                u[i,j] = u_1[i,j] + dt*V(x[i], y[j]) + \
                0.5*Cx2*(u_1[i-1,j] - 2*u_1[i,j] + u_1[i+1,j]) + \
                0.5*Cy2*(u_1[i,j-1] - 2*u_1[i,j] + u_1[i,j+1]) + \
                0.5*dt2*f(x[i], y[j], t[n])
        j = 0
        for i in range(0, Nx): u[i,j] = 0
        j = Ny
        for i in range(0, Nx): u[i,j] = 0
        i = 0
        for j in range(0, Ny): u[i,j] = 0
        i = Nx
        for j in range(0, Ny): u[i,j] = 0
    else:  # use vectorized version
        f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
        V_a = V(xv, yv)
        u[1:-1,1:-1] = u_1[1:-1,1:-1] + dt*V_a[1:-1,1:-1] + \
        0.5*Cx2*(u_1[:-2,1:-1] - 2*u_1[1:-1,1:-1] + u_1[2:,1:-1]) +\
        0.5*Cy2*(u_1[1:-1,:-2] - 2*u_1[1:-1,1:-1] + u_1[1:-1,2:]) +\
        0.5*dt2*f_a[1:-1,1:-1]
        # Boundary condition u=0
        u[: ,0] = 0
        u[:,Ny] = 0
        u[0 ,:] = 0
        u[Nx,:] = 0

    if user_action is not None:
        user_action(u, x, xv, y, yv, t, 1)

    u_2[:,:] = u_1; u_1[:,:] = u

    if version == 'scalar':
        # use f(x,y,t) function
        for n in range(1,N):
        	u = advance(u, u_1, u_2, f, x, y, t, n, Cx2, Cy2, dt2)
        	filename = 'solution_%06d.txt' % n
        	np.savetxt(filename, u, fmt='%.18e', delimiter=' ', newline='\n')
        	u_2[:,:], u_1[:,:] = u_1, u
    elif version == 'vectorized':
        for n in range(1,N):
        	f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
        	u = advance(u, u_1, u_2, f_a, x, y, t, Cx2, Cy2, dt2)
        	filename = 'solution_%06d.txt' % n
        	np.savetxt(filename, u, fmt='%.18e', delimiter=' ', newline='\n')
        	u_2[:,:], u_1[:,:] = u_1, u

    elif version == 'cpp':
    	file1 = 'uone.txt'
    	np.savetxt(file1, u_1, fmt='%.18e', delimiter=' ', newline='\n')
    	file2 = 'utwo.txt'
    	np.savetxt(file1, u_2, fmt='%.18e', delimiter=' ', newline='\n')
    	os.system('wave2d'+' '+file1+' '+file2)
    t1 = time.clock()
    # dt might be computed in this function so return the value
    return dt, t1 - t0

def advance_scalar(u, u_1, u_2, f, x, y, t, n, Cx2, Cy2, dt2):
    Nx = u.shape[0]-1;  Ny = u.shape[1]-1
    for i in range(1, Nx):
        for j in range(1, Ny):
            u[i,j] = 2*u_1[i,j] - u_2[i,j] + \
                     Cx2*(u_1[i-1,j] - 2*u_1[i,j] + u_1[i+1,j]) + \
                     Cy2*(u_1[i,j-1] - 2*u_1[i,j] + u_1[i,j+1]) + \
                     dt2*f(x[i], y[j], t[n])
    # Boundary condition u=0
    j = 0
    for i in range(0, Nx+1): u[i,j] = 0
    j = Ny
    for i in range(0, Nx+1): u[i,j] = 0
    i = 0
    for j in range(0, Ny+1): u[i,j] = 0
    i = Nx
    for j in range(0, Ny+1): u[i,j] = 0
    return u

def advance_vectorized(u, u_1, u_2, f_a, x, y, t, Cx2, Cy2, dt2):
    u[1:-1,1:-1] = 2*u_1[1:-1,1:-1] - u_2[1:-1,1:-1] + \
         Cx2*(u_1[:-2,1:-1] - 2*u_1[1:-1,1:-1] + u_1[2:,1:-1]) + \
         Cy2*(u_1[1:-1,:-2] - 2*u_1[1:-1,1:-1] + u_1[1:-1,2:]) + \
         dt2*f_a[1:-1,1:-1]
    # Boundary condition u=0
    Nx = u.shape[0]-1;  Ny = u.shape[1]-1
    u[: ,0] = 0
    u[:,Ny] = 0
    u[0 ,:] = 0
    u[Nx,:] = 0
    return u

    import nose.tools as nt

def test_quadratic(Nx=4, Ny=5):
    def exact_solution(x, y, t):
        return x*(Lx - x)*y*(Ly - y)*(1 + 0.5*t)

    def I(x, y):
        return exact_solution(x, y, 0)

    def V(x, y):
        return 0.5*exact_solution(x, y, 0)

    def f(x, y, t):
        return 2*c**2*(1 + 0.5*t)*(y*(Ly - y) + x*(Lx - x))

    Lx = 3;  Ly = 3
    c = 1.5
    dt = -1 # use longest possible steps
    T = 18

    # Note: problem with nosetests - some module import problem
    def assert_no_error(u, x, xv, y, yv, t, n):
        u_e = exact_solution(xv, yv, t[n])
        diff = abs(u - u_e).max()
        #print n, version, diff
        nt.assert_almost_equal(diff, 0, places=12)

    for version in 'scalar', 'vectorized', 'cython', 'f77', 'c_cy', 'c_f2py':
        print 'testing', version
        dt, cpu = solver(I, V, f, c, Lx, Ly, Nx, Ny, dt, T,
                         user_action=assert_no_error,
                         version=version)


def run_efficiency_tests(nrefinements=4):
    def I(x, y):
        return sin(pi*x/Lx)*sin(pi*y/Ly)

    Lx = 10;  Ly = 10
    c = 1.5
    T = 100
    versions = ['scalar', 'vectorized', 'cpp']
    print ' '*15, ''.join(['%-13s' % v for v in versions])
    for Nxy in 15, 30, 60, 120:
        cpu = {}
        for version in versions:
            dt, cpu_ = solver(I, None, None, c, Lx, Ly, Nxy, Nxy,
                              -1, T, user_action=None,
                              version=version)
            cpu[version] = cpu_
        cpu_min = min(list(cpu.values()))
        if cpu_min < 1E-6:
            print 'Ignored %dx%d grid (too small execution time)' \
                  % (Nxy, Nxy)
        else:
            cpu = {version: cpu[version]/cpu_min for version in cpu}
            print '%-15s' % '%dx%d' % (Nxy, Nxy),
            print ''.join(['%13.1f' % cpu[version] for version in versions])

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
    c = 1.0

    def I(x, y):
        """Gaussian peak at (Lx/2, Ly/2)."""
        return exp(-0.5*(x-Lx/2.0)**2 - 0.5*(y-Ly/2.0)**2)
        

    Nx = 40; Ny = 40; T = 20
    dt,time = solver(I, None, None, c, Lx, Ly, Nx, Ny, 0, T, version=version)

    if plot_method == 1:
    	mvname = 'test'
    	from mcrtmv import mcrtmv
    	print dt
    	print type(dt)
    	N = int(round(T/float(dt)))
    	mcrtmv(N, mvname, dt,Lx,Ly,Nx,Ny)



if __name__ == '__main__':
    import sys
    run_Gaussian()
    #from scitools.misc import function_UI
    #cmd = function_UI([test_quadratic, run_efficiency_tests,
    #                   run_Gaussian, ], sys.argv)
    #eval(cmd)
