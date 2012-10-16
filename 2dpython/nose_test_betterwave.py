import nose.tools as nt
import sys, os
sys.path.insert(0, os.pardir)
import betterwave as bw
import numpy as np
import matplotlib.pyplot as plt
from mcrtmv import mcrtmv

def test_constant_solution():
	"""compare to the constant solution u = 1.2"""
	Lx = 10
	Ly = 10
	def c(x,y): 
		s = (np.size(x),np.size(y))
		c = np.zeros(s)
		c = c + 1
		return c

	def I(x, y):
		"""constant solution u = 1.2, should be exact"""
		return 1.2


	Nx = 40; Ny = 40; T = 20
	u_exact = 1.2*np.ones((Nx+1,Ny+1))
	for version in ['scalar','vectorized']:
		dt,time, T2 = bw.solver(I, None, None, c, 0.0, Lx, Ly, Nx, Ny, 0.1, T, version=version)
		n = int(T/dt - 1)
		fname = 'solution_%06d.txt' % n
		u = np.loadtxt(fname)
		diff = np.abs(u_exact - u).max()
		nt.assert_almost_equal(diff, 0, delta=1E-14)
		#mcrtmv(n, dt,Lx,Ly,Nx,Ny,savemovie=True, mvname=version)
		
def test_1Dsolution():
	Lx = 10
	Ly = 10
	def c(x,y): 
		s = (np.size(x),np.size(y))
		c = np.zeros(s)
		c = c + 1
		return c

	def I(x, y):
		"""Tests the plug wave solution, should be reproduced the exact solution!"""
		return np.logical_and((np.abs(x) < Lx/2.+1), (np.abs(x) > Lx/2.-1))*1.
	Nx = 40; Ny = 40; T = 21.

	x = np.linspace(0, Lx, Nx+1)  # mesh points in x dir
	y = np.linspace(0, Ly, Ny+1)  # mesh points in y dir
	dx = x[1] - x[0]
	dy = y[1] - y[0]
	dt = dx
	X,Y = np.meshgrid(x,y)
	u_exact = I(X,Y)
	for version in ['scalar','vectorized']:
		dt,time, T = bw.solver(I, None, None, c, 0.0, Lx, Ly, Nx, Ny, dt, T, version=version)
		n = int(T/dt - 5)
		fname = 'solution_%06d.txt' % n
		u = np.loadtxt(fname)
		diff = np.abs(u_exact - u[:,0]).max()
		mcrtmv(n, dt,Lx,Ly,Nx,Ny,savemovie=True, mvname=version)
		#nt.assert_almost_equal(diff, 0, delta=1E-14)

def test_quadratic():
	"""test for convergence rate in manufactured solution"""
	Lx = 4.
	Ly = 4.
	def c(x,y): 
		s = (np.size(x),np.size(y))
		c = np.zeros(s)
		c = c + 1
		return c
	mx = 4.; my = 4.; b = 0.01; omega=np.pi;
	def u_exact(x,y,t):
		return np.exp(-b*t)*np.cos(mx*x*np.pi/Lx)*np.cos(my*y*np.pi/Ly)*np.cos(omega*t)

	def f(x,y,t):
		return np.exp(-b*t)*np.cos(mx*np.pi*x/Lx)*np.cos(my*np.pi*y/Ly)*(b*omega*np.sin(omega*t) + \
			( - omega*omega + (mx*np.pi/Lx)**2 + (my*np.pi/Ly)**2 )*np.cos(omega*t))

	def I(x, y):
		return np.cos(mx*np.pi*x/Lx)*np.cos(my*np.pi*y/Ly)

	def V(x,y):
		return -b*np.cos(mx*np.pi*x/Lx)*np.cos(my*np.pi*y/Ly)
	Nx = 40; Ny = 40; T = 1.0
	dt = 0.01
	fac = 1.
	for version in ['scalar','vectorized']:
		errorlist=[0,0,0,0]
		for i in range(3):
			
			x = np.linspace(0, Lx, Nx+1); dx = x[1]-x[0]  # mesh points in x dir
			y = np.linspace(0, Ly, Ny+1); dy = y[1]-y[0]  # mesh points in y dir
			dt,time, T2 = bw.solver(I, V, f, c, b, Lx, Ly, Nx, Ny, dt, T, version=version)
			n = int(T/dt)
			fname = 'solution_%06d.txt' % n
			u = np.loadtxt(fname)
			X,Y = np.meshgrid(x,y)
			u_ex = u_exact(X,Y,T2)
			error = np.copy((u_ex-u)**2)
			error = np.sum(error)
			error = np.sqrt(error*dx*dy)
			errorlist[i] = error
			Nx = Nx*2;
			Ny = Ny*2;
			dt = dt/2

		ratio1 = errorlist[1]/errorlist[0]
		ratio2 = errorlist[2]/errorlist[1]
		nt.assert_almost_equal(ratio1, 0.25, delta=0.06)
		nt.assert_almost_equal(ratio2, 0.25, delta=0.06)
	

		#diff = np.abs(u_exact - u).max()
		#nt.assert_almost_equal(diff, 0, delta=1E-14)
		#mcrtmv(n, dt,Lx,Ly,Nx,Ny,savemovie=False, mvname=version)

# no need for any main
	