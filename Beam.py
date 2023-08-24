#------------------------------------------------------------------------#
#                               DESCRIPTION                              #
# The Python Class Beam implement the Rayleigh-Ritz approach to evaluate #
# flexural modes based on a Timoshenko approach                          #
#------------------------------------------------------------------------#
# Authors: Cédric GIRY, Arnaud MONTABERT, Jade LEPINE                    #
#------------------------------------------------------------------------#
# Reference: Montabert, A.; Giry, C.; Limoge Schraen, C.; Lépine, J.;    #
#           Choueiri, C. Mercerat, E. D.; Guéguen, P. 
#           An Open Database to Evaluate the Fundamental Frequency of    #
#           Historical Masonry Towers through Empirical and Physics-based#
#           Formulations. Buildings 2023                                 #
#------------------------------------------------------------------------#
import numpy as np
from   scipy import *
import sympy as sym
from   sympy import lambdify, symbols
import scipy.integrate as integrate
from   scipy.integrate import quad
import scipy.linalg as la
import itertools
import matplotlib.pyplot as plt
import shutil
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class ModelBellTower:
	"""Class Model approximated over polynomial function"""
	def __init__(self, param):
		"""Initialisation avec param"""
		# Param model
		paramm    = param["Model"]
		self.H    = paramm["H"]
		self.beam = paramm["Beam_model"]

		# Param tower
		if ("Tower" in param):
			paramT   = param["Tower"]
			self.S   = paramT["S"]
			self.Iz  = paramT["Iz"]
			self.rho = paramT["rho"]
			self.E   = paramT["E"]
			# Timoshenko
			if self.beam == 'Timoshenko':
				self.G        = paramT["G"]
				self.kshear   = paramT["kshear"] # see Cowper (1966) for the implementation of formulas
			self.Kt  = self.E*self.Iz*self.H
			self.Mt  = self.rho*self.S*self.H
		else:
			self.Kt  = 0.
			self.Mt  = 0. 

		# Param polynomial basis
		if ("Poly" in param):
			paramPoly = param["Poly"]
			self.npol = paramPoly["npol"]
			self.type = paramPoly["type"]
			if self.beam == 'Timoshenko':
				if 'npolt' in paramPoly:
					self.npolt = paramPoly["npolt"]
				else:
					self.npolt = self.npol
		else:
			self.npol  = 6
			self.type  = 'Simple_poly' 

		# Param bell	
		if ("Bell" in param):
			paramb     = param["Bell"]
			self.hb    = paramb["hb"]
			self.Mb    = paramb["Mb"]
			if self.Mt ==0:
				print("Warning !!! Parameters for the tower are missing !!!")
				self.alpMb = 0.
				self.alphb = 0.
			else:
				self.alpMb = self.Mb/self.Mt
				self.alphb = self.hb/self.H
		else:
			self.alpMb = 0.
			self.alphb = 0.

		# Param TNI
		if ("TNI" in param):
			paramTNI   = param["TNI"]
			self.h     = paramTNI["h"]
			self.kN    = paramTNI["kN"]
			if self.Kt ==0:
				print("Warning !!! Parameters for the tower are missing !!!")
				self.alph = 0.
				self.alpkN = 0.
			else:
				self.alph  = self.h/self.H
				self.alpkN = self.kN/self.Kt
				#self.alpkN = self.h*self.kN/self.Kt
		else:
			self.alph  = 0.
			self.alpkN = 0.

		# Param SSI
		if ("SSI" in param):
			paramSSI   = param["SSI"]
			self.inter = 'SSI'
			self.ks    = paramSSI["ks"]
			self.kro   = paramSSI["kro"]
			self.I0t   = paramSSI["I0tow"]
			if self.Kt ==0:
				print("Warning !!! Parameters for the tower are missing !!!")
				self.alpks = 0.
			else: 
				self.alpks  = self.ks/self.Kt
				self.alpkro = self.kro/self.Kt
		else:
			self.inter  = 'No'
			self.alpks  = 0.
			self.alpkro = 0.

		self.poly_basis()
		self.nbmo = len(self.phi) #number of modes
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#                         MASS BEAM                                   
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def mij(self, phi_i, phi_j):
		"""Compute the matrix coefficient of the mass matrix
		phi_i, phi_j : shape functions for the approximation
		"""
		x = symbols('x')
		function = self.rho*self.S*phi_i*phi_j
		f = lambdify(x, function, modules=['numpy'])
		m_ij = quad(f,0,self.H)[0]
		return float(m_ij)

	def mij_timov(self, phi_i, phi_j):
		"""Compute the matrix coefficient of the mass matrix
		phi_i, phi_j : shape functions for the approximation
		"""
		x = symbols('x')
		function = self.rho*self.S*phi_i*phi_j
		f = lambdify(x, function, modules=['numpy'])
		m_ij = quad(f,0,self.H)[0]
		return float(m_ij)

	def mij_timot(self, phi_i, phi_j):
		"""Compute the matrix coefficient of the mass matrix
		phi_i, phi_j : shape functions for the approximation
		"""
		x = symbols('x')
		function = self.rho*self.Iz*phi_i*phi_j
		f = lambdify(x, function, modules=['numpy'])
		m_ij = quad(f,0,self.H)[0]
		return float(m_ij)

	def mass_beam(self):
		"""Computation of the mass matrix for a cantilever beam
		Phi: list of the mode shape function for the approximation 
		"""
		if self.beam=='Euler_Bernoulli':			
			self.M = [[self.mij(self.phi[i], self.phi[j]) for j in range(self.nbmo)] for i in range(self.nbmo)]
		elif self.beam=='Timoshenko':
			self.M = [[self.mij_timov(self.phi[i], self.phi[j]) if (i < self.nbmov and j < self.nbmov) else self.mij_timot(self.phi[i], self.phi[j]) if (i >= self.nbmov and j >= self.nbmov) else 0. for j in range(self.nbmo)] for i in range(self.nbmo)]
		return self.M

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#                         MASS BELL                                   
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def mij_b(self, phi_i, phi_j):
		"""Compute the matrix coefficient of the mass matrix for the bell
		phi_i, phi_j : shape functions for the approximation
		"""
		x = symbols('x')
		function = self.alpMb*self.Mt*phi_i*phi_j
		m_ij_b = function.subs(x, self.alphb*self.H)
		return m_ij_b

	def mass_bell(self):
		"""
		Computation of the mass matrix for the bell
		hb [meter]: altitude of the bell
		Mb [kg]: mass of the bell
		lmod: list of the function for the approximation 
		"""
		if self.beam=='Timoshenko':
			self.Massb = [[0. if (i >= self.nbmov or j >= self.nbmov) else self.mij_b(self.phi[i], self.phi[j]) for j in range(self.nbmo)] for i in range(self.nbmo)]
		else:
			self.Massb = [[self.mij_b(self.phi[i], self.phi[j]) for j in range(self.nbmo)] for i in range(self.nbmo)]
		return self.Massb

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#                         BEAM STIFFNESS                                   
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def kij(self, phi_i, phi_j):
		"""Compute the matrix coefficient of the stiffness matrix
		phi_i, phi_j : shape functions for the approximation
		dphi_i, dphi_j: first derivative of shape function with respect to x
		ddphi_i, ddphi_j: second derivative of shape function with respect to x
		"""
		x = symbols('x')
		dphi_i   = sym.diff(phi_i, x)
		dphi_j   = sym.diff(phi_j, x)
		ddphi_i  = sym.diff(dphi_i, x)
		ddphi_j  = sym.diff(dphi_j, x)
		function = self.E*self.Iz*ddphi_i*ddphi_j
		f    = lambdify(x, function, modules=['numpy'])
		k_ij = quad(f,0, self.H)[0]
		return k_ij

	def kij_timovv(self, phi_i, phi_j):
		"""Compute the matrix coefficient of the stiffness matrix
		phi_i, phi_j : shape functions for the approximation
		dphi_i, dphi_j: first derivative of shape function with respect to x
		"""
		x = symbols('x')
		dphi_i = sym.diff(phi_i, x)
		dphi_j = sym.diff(phi_j, x)
		function = self.kshear*self.G*self.S*dphi_i*dphi_j
		f = lambdify(x, function, modules=['numpy'])
		k_ij = quad(f,0, self.H)[0]
		return k_ij

	def kij_timott(self, phi_i, phi_j):
		"""Compute the matrix coefficient of the stiffness matrix
		phi_i, phi_j : shape functions for the approximation
		dphi_i, dphi_j: first derivative of shape function with respect to x
		"""
		x = symbols('x')
		dphi_i = sym.diff(phi_i, x)
		dphi_j = sym.diff(phi_j, x)
		function = self.E*self.Iz*dphi_i*dphi_j+self.kshear*self.G*self.S*phi_i*phi_j
		f = lambdify(x, function, modules=['numpy'])
		k_ij = quad(f,0, self.H)[0]
		return k_ij

	def kij_timovt(self, phi_i, phi_j):
		"""Compute the matrix coefficient of the stiffness matrix
		phi_i, phi_j : shape functions for the approximation
		dphi_i, dphi_j: first derivative of shape function with respect to x
		"""
		x = symbols('x')
		dphi_i = sym.diff(phi_i, x)
		function = -self.kshear*self.G*self.S*phi_j*dphi_i
		f = lambdify(x, function, modules=['numpy'])
		k_ij = quad(f,0, self.H)[0]
		return k_ij

	def rigi_beam(self):
		"""Computation of the stiffness matrix for a cantilever beam
		Phi: list of the mode shape function for the approximation 
		"""
		if self.beam=='Euler_Bernoulli':
			self.k = [[self.kij(self.phi[i], self.phi[j]) for j in range(self.nbmo)] for i in range(self.nbmo)]
		elif self.beam=='Timoshenko':
			self.k = [[self.kij_timovv(self.phi[i], self.phi[j]) if (i < self.nbmov and j < self.nbmov) else self.kij_timott(self.phi[i], self.phi[j]) if (i >= self.nbmov and j >= self.nbmov) else self.kij_timovt(self.phi[i], self.phi[j]) if (i < self.nbmov and j >= self.nbmov) else self.kij_timovt(self.phi[j], self.phi[i]) for j in range(self.nbmo)] for i in range(self.nbmo)]
		return self.k

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#                         ADJACENT BUILDING INTERACTION STIFFNESS                        
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def kij_nave(self, phi_i, phi_j):
		"""Compute the matrix coefficient of the stiffness matrix interacting with the nave
		phi_i, phi_j : shape functions for the approximation
		"""
		x = symbols('x')
		kn = self.alpkN*self.Kt
		function = kn*phi_i*phi_j
		f = lambdify(x, function, modules=['numpy'])
		hn = self.alph*self.H
		k_ij = quad(f,0,hn)[0]
		return k_ij

	def rigi_nave(self):
		if self.beam=='Timoshenko':
			self.k_nave = [[self.kij_nave(self.phi[i], self.phi[j]) if (j<self.nbmov and i<self.nbmov) else 0. for j in range(self.nbmo)] for i in range(self.nbmo)]
		else:
			self.k_nave = [[self.kij_nave(self.phi[i], self.phi[j]) for j in range(self.nbmo)] for i in range(self.nbmo)]
		return self.k_nave

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#                         SOIL STRUCTURE INTERACTION STIFFNESS                        
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
	def rij_soilt(self,phi_i, phi_j):
		"""Compute the matrix coefficient of the stiffness matrix interacting with the soil
		phi_i, phi_j : shape functions for the approximation
		"""
		x = symbols('x')
		function  = self.alpks*self.Kt*phi_i*phi_j
		k_ij_soil = function.subs(x, 0.)
		return k_ij_soil

	def rij_soilr(self,phi_i, phi_j):
		"""Compute the matrix coefficient of the stiffness matrix interacting with the soil
		phi_i, phi_j : shape functions for the approximation
		"""
		x = symbols('x')
		dphi_i = sym.diff(phi_i, x)
		dphi_j = sym.diff(phi_j, x)
		if self.beam=='Timoshenko':
			function  = self.alpkro*self.Kt*phi_i*phi_j
		else:
			function  = self.alpkro*self.Kt*dphi_i*dphi_j
		k_ij_soil = function.subs(x, 0.)
		return k_ij_soil

	def rigi_soil(self):
		if self.beam=='Timoshenko':
			self.k_soil = [[self.rij_soilt(self.phi[i], self.phi[j]) if (i<self.nbmov and j<self.nbmov) else self.rij_soilr(self.phi[i], self.phi[j]) for j in range(self.nbmo)] for i in range(self.nbmo)]
		else:
			self.k_soil = [[(self.rij_soilt(self.phi[i], self.phi[j])+self.rij_soilr(self.phi[i], self.phi[j])) for j in range(self.nbmo)] for i in range(self.nbmo)]
		return self.k_soil

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#                         TOTAL STIFFNESS                        
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def rigi_tot(self):
		"""
		Computation of the total stiffness matrix over the approximated basis (polynomial basis)
		"""
		self.Ktot = np.array(self.rigi_beam(), dtype = float) + np.array(self.rigi_nave(), dtype = float) + np.array(self.rigi_soil(), dtype = float)
		return self.Ktot

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#                         TOTAL MASS                        
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def mass_tot(self):
		"""
		Computation of the total mass matrix over the 
		approximated basis (polynomial basis)
		"""
		self.Mtot = np.array(self.mass_beam(), dtype = float) + np.array(self.mass_bell(), dtype = float)
		return self.Mtot

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#                         Polynomial basis                        
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	def poly_basis(self):
		"""
		Definition of the approximated basis
		"""
		x = sym.symbols('x')
		if self.type == 'Euler_poly': # Mazanoglu (2015)
			if self.inter == 'SSI':
				self.phi = [(1-(x/self.H))**i for i in range(self.npol)]
			else:
				self.phi = [(x/self.H)**2*(1-(x/self.H))**i for i in range(self.npol)]
		if self.type == 'Timo_poly':   # Mazanoglu (2015)
			if self.inter == 'SSI':
				self.phi = [(1-(x/self.H))**i for i in range(self.npol)]
				self.phi+= [(1-(x/self.H))**i for i in range(self.npolt)]
				#self.phi = [(x/self.H)**2*(1-(x/self.H))**i for i in range(self.npol)]
				#self.phi+= [(x/self.H)*(1-(x/self.H))**i for i in range(self.npolt)]
				self.nbmov = self.npol
			else:
				self.phi = [(x/self.H)**2*(1-(x/self.H))**i for i in range(self.npol)]
				self.phi+= [(x/self.H)*(1-(x/self.H))**i for i in range(self.npolt)]
				self.nbmov = self.npol
		return self.phi

	def modal_analysis(self):
		"""
		Modal analysis over the approximated basis
		"""
		Ksym = self.rigi_tot()
		Msym = self.mass_tot() 
		om2, self.factorpoly= la.eigh(Ksym,Msym, eigvals_only=False)

		#normalizing eigenvector
		#for i in range(len(self.factorpoly)):
		#	self.factorpoly[:, i] *= (1./np.linalg.norm(self.factorpoly[:, i]))
		self.factorpoly = self.factorpoly.T

		#filtering real eigenvalues and associated eigen_vectors
		selectors = [item for item in om2.real]
		self.freq = (1./(2*np.pi))*np.sqrt(list(itertools.compress(om2, selectors)))

		self.factorpoly= list(itertools.compress(self.factorpoly, selectors))
		self.freq, self.factorpoly = (list(item) for item in zip(*sorted(zip(self.freq, self.factorpoly))))

		#building mode shape
		if self.beam=='Timoshenko':
			defov = [np.sum([self.factorpoly[j][i]*self.phi[i] if i<self.nbmov else 0.*self.phi[i] for i in range(len(self.phi))]) for j in range(self.npol)]
			defor = [np.sum([self.factorpoly[j][i]*self.phi[i] if i>=self.nbmov else 0.*self.phi[i] for i in range(len(self.phi))]) for j in range(self.npol)]
			self.defo = [defov,defor]
		else:
			self.defo = [np.sum([self.factorpoly[j][i]*self.phi[i] for i in range(len(self.phi))]) for j in range(np.shape(self.factorpoly)[0])]

		return self.freq, self.factorpoly, self.defo

if __name__ == '__main__':
	
	# Case study (cas_study)
	# ----------------------
	# 1. Cantilever beam
	# 2. Cantilever beam with interaction with nave
	# 3. Cantilever beam with interaction with nave and soil
	cas_study = 3
	
	# - Bell tower characteristics (paramT)
	# -------------------------------------
	# S   [m^2]    : section of the tower
	# Iz  [m^4]    : second moment of area for bending axis 
	# rho [kg/m^3] : density of masonry
	# E   [N/m^2]  : Young modulus of masonry
	# Mt  [kg]     : mass of the tower
	# Kt  [N/m]    : rigidity of the tower
	# H   [meter]  : height of the tower
	H        = 15
	model_ty = 'Timoshenko'
	#model_ty = 'Euler_Bernoulli'
	paramM = {"H":H, "Beam_model":model_ty}     #Model
	ep  = 0.9 ; B =  3.2; W = 3.;
	S   = B*W - ((B-2*ep)*(W-2*ep));
	Iz  = (B*(W**3))/12 - ((B-2*ep)*((W-2*ep)**3))/12
	rho = 2200. ; E   = 2.5e9 ; nu = 0.3 ; G = E/(2*(1+nu));
	m = ((B-2*ep)*ep)/(W*ep)
	n = (B-2*ep)/W
	# Cowper 1966
	kshear = (10*(1+nu)*(1+3*m)**2)/((12+72*m+150*m**2+90*m**3)+nu*(11+66*m+135*m**2+90*m**3)+10*n**2*((3+nu)*m+3*m**2))
	#paramT = {"S":S,"Iz":Iz,"rho":rho,"E":E}   #Bell Tower
	paramT = {"S":S,"Iz":Iz,"rho":rho,"E":E,"kshear":kshear,"G": G}   #Bell Tower

	# - Bell characteristics
	# -----
	# hb [meter]: altitude of the bell
	# Mb [kg]: mass of the bell
	hb = H ; Mb = 1000. ;
	paramb = {"hb":hb,"Mb":Mb}                 #Bell system

	# SSI characteristics
	# ---
	# ks  [N/m]: translational stiffness of soil/structure
	# kro [N.m]: rotational stiffness of soil/structure	
	# I0tow [kg.m2] : Inertia of the tower wrt to the z axis at the soil at midwidth
	#ks    = 1.E10 ;
	#kro   = 1.E7 ;
	ks    = 1.E7 ;
	kro   = 1.E8 ;
	I0tow = rho*H*((B*H**2*W/3+W**2/12)-((B-2*ep)*H**2*(W-2*ep)/3+(W-2*ep)**2/12))
	paramssi={"ks":ks,"kro":kro,"I0tow":I0tow}                         #SSI

	# Adjacent building interaction
	# ---
	# h [meter]: height of the adjacent building in interaction with the tower
	# kN [N/m2]: rigidity of the ineraction between the adjacent building and the tower
	#
	h  = 10.; kN = 1.e9/h;
	paramtni = {"kN":kN,"h":h}                 #Adjacent building interaction
	
	# - Polynomial basis
	# -----
	# npol : number of polynoms in the basis
	# type : type of polynoms ('Euler_poly', 'Timo_poly' )
	#npol = 10 ; typep = 'Euler_poly' ;
	#paramPoly = {"npol":npol,"type":typep}
	npol = 10 ; typep = 'Timo_poly' ;
	paramPoly = {"npol":npol,"type":typep,"npolt":npol}

	# Compilation of the parameters depeding of the case study
	# --------------------------------------------------------
	# Euler_Bernoulli
	# 1.simple tower
	# 2.simple tower + Nave interaction
	# 3.simple tower + Bell
	# 4.simple tower + SSI
	# 5.full
	#param={"Model":paramM,"Tower":paramT,"Poly":paramPoly}
	#param={"Model":paramM,"Tower":paramT,"TNI":paramtni,"Poly": paramPoly}
	#param={"Model":paramM,"Tower":paramT,"Bell":paramb, "Poly": paramPoly}
	#param={"Model":paramM,"Tower":paramT,"SSI":paramssi,"Poly": paramPoly}
	#param={"Model":paramM,"Tower":paramT,"TNI":paramtni,"Bell":paramb,"SSI":paramssi,"Poly": paramPoly}

	# Timoshenko
	# 1.simple tower
	# 2.simple tower + Nave interaction
	# 3.simple tower + SSI
	# 4.simple tower + Bell
	# 5.full
	if cas_study==1:
		param={"Model":paramM,"Tower":paramT,"Poly": paramPoly}
	elif cas_study==2:
		param={"Model":paramM,"Tower":paramT,"TNI":paramtni,"Poly": paramPoly}
	elif cas_study==3:
		param={"Model":paramM,"Tower":paramT,"TNI":paramtni,"SSI":paramssi,"Poly": paramPoly}
	else:
		disp('Case not implemented')
	#param={"Model":paramM,"Tower":paramT,"Bell":paramb,"Poly": paramPoly}
	#param={"Model":paramM,"Tower":paramT,"SSI":paramssi,"Poly": paramPoly}
	#param={"Model":paramM,"Tower":paramT,"TNI":paramtni,"SSI":paramssi,"Poly": paramPoly}
	#param={"Model":paramM,"Tower":paramT,"TNI":paramtni,"Bell":paramb,"SSI":paramssi,"Poly": paramPoly}

	# DEFINITION OF THE MODEL
	# -----------------------
	Model = ModelBellTower(param)

	# MODAL ANALYSIS
	# --------------
	freq, shape, defo = Model.modal_analysis()

	# POST-TREATMENT
	# --------------
	li_z = np.linspace(0, H, 40)
	x = sym.symbols('x')
	test_figure = True
	if test_figure:
		plt.figure()
		for i in range(3):
			text = 'Mode '+str(i+1)+' - f: '+str(round(freq[i],2))+' Hz'
			if Model.beam == 'Timoshenko':
				mode = [defo[0][i].subs(x, value) for value in li_z]
			else:
				mode = [defo[i].subs(x, value) for value in li_z]
			plt.plot([item*1./np.max(np.abs(mode)) for item in mode], li_z,'o-',label = text)
		plt.legend()
		nfic = 'Case_'+str(cas_study)+'modes.png'
		plt.savefig(nfic,dpi=1000)
		#plt.show()
