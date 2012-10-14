import nose.tools as nt
import sys, os
sys.path.insert(0, os.pardir)
import betterwave as bw
import numpy as np
import matplotlib.pyplot as plt
from mcrtmv import mcrtmv

def notest_constant_solution():
	"""compare to the constant solution u = 1.2"""
	Lx = 10
	Ly = 10
	def c(x,y): 
		s = (np.size(x),np.size(y))
		c = np.zeros(s)
		c = c + 1
		return c

	def I(x, y):
		"""ballemor"""
		"""constant solution u = 1.2"""
		return 1.2


	Nx = 40; Ny = 40; T = 20
	u_exact = 1.2*np.ones((Nx+1,Ny+1))
	for version in ['vectorized']:
		dt,time = bw.solver(I, None, None, c, 0.0, Lx, Ly, Nx, Ny, 0.1, T, version=version)
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
		"""constant solution u = 1.2"""
		return np.logical_and((np.abs(x) < 5), (np.abs(x) > 3))*1.
	Nx = 40; Ny = 40; T = 20
	u_exact = 1.2*np.ones((Nx+1,Ny+1))
	for version in ['scalar']:
		dt,time = bw.solver(I, None, None, c, 0.0, Lx, Ly, Nx, Ny, 0, T, version=version)
		n = int(T/dt - 1)
		fname = 'solution_%06d.txt' % n
		u = np.loadtxt(fname)
		diff = np.abs(u_exact - u).max()
		#nt.assert_almost_equal(diff, 0, delta=1E-14)
		mcrtmv(n, dt,Lx,Ly,Nx,Ny,savemovie=False, mvname=version)

def test_quadratic():
	"""test for convergence rate in manufactured solution"""
	


'''
def test_exact_solution_linear():
	"""compare to exact solution in linear case, should be equal to machine precission."""
	v0, rho, T, t0, rho_b, mu, body, diameter, volume, drag_coeff, cross_section, makeplot, dt_values = vt.read_command_line()
	rho = 0.0
	mass = rho_b*volume
	dt = dt_values[0]
	waterdrop = vt.Body(volume, mass, diameter, cross_section, drag_coeff, v0)
	air = vt.System(rho,mu)
	waterdrop.time_solution(air,t0,T,dt)
	exact_solution = np.zeros(len(waterdrop.t))  #= v0 - waterdrop.t*air.g
	exact_solution[0] = v0
	for i in range(len(waterdrop.t)-1):	
		exact_solution[i+1] = exact_solution[i] - air.g*dt
	diff = np.abs(waterdrop.v - exact_solution).max()	
	nt.assert_almost_equal(diff,0,delta=1E-14)

def test_convergence():
	v0, rho, T, t0, rho_b, mu, body, diameter, volume, drag_coeff, cross_section, makeplot, dt_values = vt.read_command_line()
	T = 40.
	mass = rho_b*volume
	dt = dt_values[0]
	waterdrop = vt.Body(volume, mass, diameter, cross_section, drag_coeff, v0)
	air = vt.System(rho,mu)
	waterdrop.time_solution(air,t0,T,dt)
	a = 0.5*waterdrop.drag_coeff*air.rho*waterdrop.cross_section/(waterdrop.rho_b*waterdrop.volume)
	b = air.g*(air.rho/waterdrop.rho_b -1)
	convergence_value = -np.sqrt(abs(b)/a)
	diff = abs(convergence_value - waterdrop.v[-1])
	nt.assert_almost_equal(diff, 0, delta=1E-14)
'''
# no need for any main