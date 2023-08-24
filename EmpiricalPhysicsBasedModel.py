#Class to deal with empirical relation
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy import stats
import ipdb as db
import sys
from figure_rule import *

#default Parameters to plot
mpl.rcParams['lines.linewidth'] = 2.0
mpl.rcParams.update({'font.size': 14})


class Formulation():
	"""Attributing empirical and physic based model to a single slender structure"""
	def __init__(self, slender_st):
		"""
		single attribut class : structures 
		"""
		self.slender_st        = slender_st
		self.shape             = self.slender_st['shape']
		self.height            = float(self.slender_st['H'])
		self.effective_height  = float(self.slender_st['Heff'])
		self.breadth           = float(self.slender_st['breadth'])
		self.length            = float(self.slender_st['length'])
		self.max_thickness     = float(self.slender_st['max_wall_thickness'])
		self.f0                = float(self.slender_st['f0'])
		try:
			self.young = 1e9*float(self.slender_st['E'])
		except:
			self.young   = np.max([1e9*float(item) for item in self.slender_st['E'].split(',')])

		try:	
			self.density = 100*float(self.slender_st['density'])
		except:
			self.density = np.max([100*float(item) for item in self.slender_st['density'].split(',')])


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# ~~~~~  Parameters for empirical models  ~~~~~
		# ~~~~~   cv Table ? in paper doi:        ~~~~~
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# A test is processed to check if parameters exist
		self.exist_height            = self.height != -1 
		self.exist_effective_height  = self.effective_height != -1 
		self.exist_breadth           = self.breadth != -1
		self.exist_length           = self.length != -1
		self.exist_max_thickness    = self.max_thickness != -1
		self.exist_E                 = self.young       != -1
		self.exist_rho               = self.density     != -1

		self.emp_model     = {'model_1_1' :{'a1' : [20.],        'b1' : [-3./4.], 'a1s' : [0.], 'b1s' : [0.],     'b2s' : [0.],  'b1e': [0.], 'exist': [self.exist_height]},
		                  'model_1_2' :{'a1' : [1./0.0187],  'b1' : [-1.],    'a1s' : [0.], 'b1s' : [0.],     'b2s' : [0.],  'b1e': [0.],    'exist': [self.exist_height]},
		                  'model_1_3' :{'a1' : [1./0.01137], 'b1' : [-1.138], 'a1s' : [0.], 'b1s' : [0.],     'b2s' : [0.],  'b1e': [0.]   , 'exist': [self.exist_height]},
		                  'model_1_4' :{'a1' : [1./0.0151],  'b1' : [-1.08],  'a1s' : [0.], 'b1s' : [0.],     'b2s' : [0.],  'b1e': [0.]   , 'exist': [self.exist_height]},
		                  'model_1_5' :{'a1' : [28.35],      'b1' : [-0.83],  'a1s' : [0.], 'b1s' : [0.],     'b2s' : [0.],  'b1e': [0.]   , 'exist': self.exist_height},
		                  'model_1_6' :{'a1' : [135.343],    'b1' : [-1.32],  'a1s' : [0.], 'b1s' : [0.],     'b2s' : [0.],  'b1e': [0.]   , 'exist': [self.exist_height]},
		                  'model_2' :{'a1' : [3.58],       'b1' : [0.],     'a1s' : [0.], 'b1s' : [0.57],   'b2s' : [0.],  'b1e': [0.]   , 'exist': [self.exist_height, self.exist_breadth]},
		                  'model_3' :{'a1' : [208.54],     'b1' : [-1.18],  'a1s' : [0.], 'b1s' : [0.55],   'b2s' : [0.],  'b1e': [0.]   , 'exist': [self.exist_breadth]},
		                  'model_4_1' :{'a1' : [1./0.06],    'b1' : [-0.5],   'a1s' : [2.], 'b1s' : [0.5],    'b2s' : [0.5], 'b1e': [0.]    , 'exist': [self.exist_height, self.exist_breadth]},
		                  'model_4_2':{'a1' : [1./0.03],    'b1' : [-0.83],  'a1s' : [1.], 'b1s' : [0.17],   'b2s' : [0.5], 'b1e': [0.]   , 'exist': [self.exist_height, self.exist_breadth]},
		                  'model_5':{'a1' : [1./0.0117],    'b1' : [0],  'a1s' : [-9.632], 'b1s' : [3],   'b2s' : [-1], 'a2s' : [94.786], 'a3s' : [144.461] ,'b1e': [0.]   , 'exist': [self.exist_height, self.exist_breadth]},
		                  'model_6':{'a1' : [12.96],      'b1' : [-0.686], 'a1s' : [0.], 'b1s' : [0.],     'b2s' : [0.],  'b1e': [-0.686], 'exist': [self.exist_height, self.exist_breadth, self.exist_effective_height]},
		                  'model_7':{'a1' : [14.61],      'b1' : [-0.811], 'a1s' : [0.], 'b1s' : [-0.254], 'b2s' : [0.],  'b1e': [-0.341], 'exist': [self.exist_height, self.exist_breadth, self.exist_effective_height]},
		                 }
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# ~~~~~  Parameters for physics based models  ~~~~~
		# ~~~~~       cv Table ? in paper doi:        ~~~~~
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		self.pb_model  = {'model_1' : {'C1': np.sqrt(1.375),  'C2' : 0, 'C3' : 1, 'exist' :[self.exist_breadth, self.exist_length, self.exist_max_thickness, self.exist_height, self.exist_E, self.exist_rho]},
	                          'model_2' : {'C1': 0.8       ,  'C2' : 1, 'C3' : 1, 'exist' :[self.exist_breadth, self.exist_length, self.exist_max_thickness, self.exist_height, self.exist_E, self.exist_rho, self.exist_effective_height]},
	                          'model_3' : {'C1': 0.8       ,  'C2' : 0, 'C3' : 1, 'exist' :[self.exist_breadth, self.exist_length, self.exist_max_thickness, self.exist_height, self.exist_E, self.exist_rho]},
	                          'model_4' : {'C1': 800       ,  'C2' : 0, 'C3' : 0, 'exist' :[self.exist_breadth, self.exist_length, self.exist_max_thickness, self.exist_height]}
				}


	def f0_emp(self, model_name):
		""" 
		General formulation of the fundamental frequency for empirical models
		model_name : Name of the physic based model to use
		H [m]: height of the tower
		ls[m]: breadth (characteristic size of the section, minimum size)
		hn[m]: interaction height between the tower and the adjacent structure
		"""

		param = self.emp_model[model_name]
		
		if False not in param['exist']:
		
			# Parameters of the empirical model
			# ---------------------------------
			param = self.emp_model[model_name]
			a1    = param['a1']
			b1    = param['b1']
			a1s   = param['a1s']
			b1s   = param['b1s']
			b2s   = param['b2s']
			b1e   = param['b1e']

			alphal = np.divide(self.breadth,self.height)
			alphah = np.divide(self.height - self.effective_height, self.height)
			#print('Empirical frequency succesfully computed from model: %s' %model_name)

			return np.multiply(a1,np.power(self.height,b1)*np.power(alphal,b1s)*np.power(np.add(1.,np.multiply(a1s,alphal)),b2s)*np.power(np.subtract(1,alphah),b1e))

		else:
			print("Warning: Empirical formulation cannot be used. A feature is missing.")
			return [None]
			pass

	
	def f0_phy(self, model_name, theta):
		""" 
		General formulation of the fundamental frequency for physics based models
		model_name : Name of the physic based model to use
		theta: rotation angle to compute the radius of inertia 
		"""
		# Parameters of the physics based model
		# ---------------------------------
		param = self.pb_model[model_name]
		c1    = param['C1']
		c2    = param['C2']
		c3    = param['C3']

		self.theta = theta

		if False not in param['exist']:
		
			self.alpha_t = self.max_thickness/self.breadth
			self.alpha_L = self.length/self.breadth
			self.alphah = np.divide(self.height - self.effective_height,self.height)
		
			# Compute the surface and Second moment of inertia
			if (self.shape == 'REC') or (self.shape == 'SQ'):
				self.alpha_shs = 1
				self.alpha_shi = 1/12
			else:
				self.alpha_shs = np.pi/4
				self.alpha_shi = np.pi/64
				
			S = 2*self.alpha_shs*(self.breadth**2)*(self.alpha_L + 1 - 2*self.alpha_t)
			self.alpha_S = 2*self.alpha_shs*(self.alpha_L + 1 - 2*self.alpha_t)
			self.alpha_Ix = self.alpha_shi*(self.alpha_L - (self.alpha_L - 2*self.alpha_t)*((1-2*self.alpha_t)**3))
			self.alpha_Iy = self.alpha_shi*(self.alpha_L**3 - ((self.alpha_L - 2*self.alpha_t)**3)*(1-2*self.alpha_t))
				

			# Compute radius of inertia
			if model_name == 'model_1':
				r = self.breadth*np.sqrt(np.divide(1,self.alpha_S))*np.sqrt(np.divide(self.alpha_Ix + self.alpha_Iy, 2) + (np.divide(self.alpha_Ix - self.alpha_Iy, 2)* np.cos(2*self.theta)))

			elif model_name == 'model_2':
				r = np.divide(self.breadth, np.sqrt(12))*1.5*(1-self.alpha_t)

			elif model_name == 'model_3':
				r = np.divide(self.breadth, np.sqrt(12))*1.125
			
			elif model_name == 'model_4':
				r = np.divide(self.breadth, np.sqrt(12))*1.125

			self.alpha_h = np.divide(self.height-self.effective_height,self.height)

			return c1*np.divide(1.875**2, 2*np.pi)*np.divide(r,np.power(self.height,2))*np.power(np.divide(1, 1 - self.alpha_h),c2)*np.power(np.sqrt(np.divide(self.young, self.density)),c3)

		else:
			print("Warning: Physics based formulation cannot be used. A feature is missing.")
			return None
			pass




class Regression():
	"""Attributing empirical and physic based model to a single slender structure"""
	def __init__(self, database):
		self.database = database

	def model_1_emp(self, x, a1, b1):
		"""
		Empirical model 1: Compute the regression analysis between the height (input) and the fundamental frequency (output)
		Correspond to model 1_1 -> model 1_6 in Montabert et al, 2023
		References: Eurocode 8;  Faccio et al, 2010; Rainieri et al, 2012; Shakya et al, 2016, Diafeiro et al, 2018.
		"""
		return np.multiply(a1, np.power(x, b1))
		
	def model_2_emp(self, x, a1, b1s):
		"""
		Empirical model 2: Compute the regression analysis between the height, the breadth (input) and the fundamental frequency (output)
		Correspond to model 2 in Montabert et al, 2023
		References: Shakya et al, 2016.
		"""
		return np.multiply(a1, np.power(x, b1s))

	def model_3_emp(self, x, a1, b1, b1s):
		"""
		Empirical model 3: Compute the regression analysis between the height, the breadth (input) and the fundamental frequency (output)
		Correspond to model 3 in Montabert et al, 2023
		References: Diafeiro et al, 2018.
		"""
		return np.multiply(np.multiply(a1, np.power(x[:,0], b1)), np.power(x[:,1], b1s))

	def model_4_emp(self, x, a1, b1, b1s, a1s, b2s):
		"""
		Empirical model 4_1 and model 4_2: Compute the regression analysis between the height, the breadth (input) and the fundamental frequency (output)
		Correspond to model 9 in Montabert et al, 2023
		References: NCSR-02
		"""
		return np.multiply(np.multiply(a1, np.power(x[:,0], b1)), np.multiply(np.power(x[:,1], b1s), np.power(1+(np.multiply(a1s,x[:,1])), b2s)))
		#return a1*(x[:,0]**b1)*(x[:,1]**b1s)*((1+(a1s*x[:,1]))**b2s)

	def model_5_emp(self, x, a1, b1s, a1s, b2s, a2s, a3s):
		"""
		Empirical model 5: Compute the regression analysis between the height, the breadth (input) and the fundamental frequency (output)
		Correspond to model 5 in Montabert et al, 2023
		References: Formisano et al, 2017.
		"""
		return np.multiply(a1, np.multiply(np.power(x[:,1], b1s), np.power(1+(a1s*x[:,1]) + np.multiply(a2s, np.power(x[:,1],2)) + np.multiply(a3s, np.power(x[:,1],3)), b2s)))
		#return a1*(x[:,1]**b1s)*((1+(a1s*x[:,1]) + (a2s*x[:,1])**2 + (a3s*x[:,1])**3)**b2s)

	def model_6_emp(self, x, a1, b1, b1e):
		"""
		Empirical model 6 : Compute the regression analysis between the height, the breadth (input) and the fundamental frequency (output)
		Correspond to model 6 in Montabert et al, 2023
		References: Diafeiro et al, 2018.
		"""
		#return np.multiply(np.multiply(a1, np.power(x[:,0], b1)), np.power( 1 - x[:,1], b1e))
		return a1*(x[:,0]**b1)*((1-x[:,1]))**b1e

	def model_7_emp(self, x, a1, b1, b1s, b1e):
		"""
		Empirical model 13 : Compute the regression analysis between the height, the breadth (input) and the fundamental frequency (output)
		Correspond to model 13 in Montabert et al, 2023
		References: Diafeiro et al, 2018.
		"""
		return np.multiply(np.multiply(a1, np.power(x[:,0], b1)), np.multiply(np.power(x[:,1], b1s) ,np.power( 1 - x[:,2], b1e)))


	def R_squarred(self, model, experiment):
		"""Compute R-squarred for a distribution"""
		self.SSR = np.sum((experiment - model)**2)
		self.SST = np.sum((experiment - np.mean(experiment))**2)

		return 1 - np.divide(self.SSR, self.SST)
	
	def power_law(self, x, a, b):
		""" Power law function for curve fitting"""
		return np.divide(1, a*np.power(x, b))

	def fit_power_model(self, model):
		"""Fit relation"""

		ydata = self.database['f0']
		if model == 'model_1':
			# ~~ Step1: clean database ~~
			self.database = self.database[(self.database['H'] != -1) & (self.database['f0'] != -1)]
			self.output = np.asarray(self.database['f0'])
			self.input = np.asarray(self.database['H'])
			self.f = self.model_1_emp

		elif model == 'model_2':
			# ~~ Step1: clean database ~~
			self.database = self.database[(self.database['H'] != -1) & (self.database['breadth'] != -1) & (self.database['f0'] != -1)]
			self.output = np.asarray(self.database['f0'])
			self.input =  np.asarray(np.divide(self.database['breadth'], self.database['H']))
			self.f = self.model_2_emp
			
		elif model == 'model_3':
			# ~~ Step1: clean database ~~
			self.database = self.database[(self.database['H'] != -1) & (self.database['breadth'] != -1) & (self.database['f0'] != -1)]
			self.output = np.asarray(self.database['f0'])
			self.input = np.vstack((np.asarray(self.database['H']), np.asarray(np.divide(self.database['breadth'], self.database['H'])))).T
			self.f = self.model_3_emp

		elif model == 'model_4':
			# ~~ Step1: clean database ~~
			self.database = self.database[(self.database['H'] != -1) & (self.database['breadth'] != -1) & (self.database['f0'] != -1)]
			self.output = np.asarray(self.database['f0'])
			self.input = np.vstack((np.asarray(self.database['H']), np.asarray(np.divide(self.database['breadth'], self.database['H'])))).T
			self.f = self.model_4_emp

		elif model == 'model_5':
			# ~~ Step1: clean database ~~
			self.database = self.database[(self.database['H'] != -1) & (self.database['breadth'] != -1) & (self.database['f0'] != -1)]
			self.output = np.asarray(self.database['f0'])
			self.input = np.vstack((np.asarray(self.database['H']), np.asarray(np.divide(self.database['breadth'], self.database['H'])))).T
			self.f = self.model_5_emp

		elif model == 'model_6':
			# ~~ Step1: clean database ~~
			self.database = self.database[(self.database['H'] != -1) & (self.database['Heff'] != -1) & (self.database['Heff'] != 0) & (self.database['f0'] != -1)]
			self.output = np.asarray(self.database['f0'])
			self.input = np.vstack((np.asarray(self.database['H']), np.asarray(np.divide(self.database['H'] - self.database['Heff'], self.database['H'])))).T
			self.f = self.model_6_emp

		elif model == 'model_7':
			# ~~ Step1: clean database ~~
			self.database = self.database[(self.database['H'] != -1) & (self.database['Heff'] != -1) & (self.database['Heff'] != 0) & (self.database['f0'] != -1) & (self.database['breadth'] != -1)]
			self.output = np.asarray(self.database['f0'])
			self.input = np.vstack((np.asarray(self.database['H']), np.asarray(np.divide(self.database['breadth'] , self.database['H'])), np.asarray(np.divide(self.database['H'] - self.database['Heff'], self.database['H'])))).T

			self.f = self.model_7_emp


		#Curve fit
		self.pars, self.cov = curve_fit(f=self.f, xdata= self.input, ydata= self.output, maxfev = 100000000)

		#standard error of parameters
		self.stderr = np.sqrt(np.diag(self.cov))

		
		#compute residuals
		self.res = self.database['f0'] - self.f(self.input, *self.pars)

		#Compute R-squarred
		self.R2 = self.R_squarred(self.f(self.input, *self.pars), self.database['f0'])
		
		return self.pars, self.cov, self.stderr, self.res, self.R2


	def IC(self, x, func):
		"""Compute interval confidence
		x : dedicated interval
		
		"""
		#correlated values
		a, b = unc.correlated_values(self.pars, self.cov)

		#distributed model
		cy   = func(x, a, b)

		#nominal value and standard deviation
		self.nom = unp.nominal_values(cy)
		self.std = unp.std_devs(cy)

		return self.nom, self.std



if __name__ == "__main__":

	# Load the masonry database
	path_data = 'TURRIS.xlsx'
	database = pd.read_excel(path_data)

	"""# TEST 1: Empirical and physics based formulation to evaluate the fundamental frequency
	# Choose your model
	model_emp = 'model_1'
	model_pb  = 'model_2'

	# Read one tower at a time
	for ii in range(len(database)):
		print("~~~~~~~~~~~~  Tower # %s  ~~~~~~~~~~~" % str(ii))
		tower = Formulation(database.iloc[ii]) # Build Empirical formulation object
		print("Tower name: %s. References: %s" % (tower.slender_st['building_name'], tower.slender_st['references']))
		# ~~~~ Empirical formulation of f0 ~~~~
		f0_model_emp = tower.f0_emp(model_emp)
		if f0_model_emp[0] != None:
			print('Empirical model: %s, frequency: %s Hz' % (model_emp, f0_model_emp[0]))
			print('Relative error between exp. f0 and emp. f0: %s' % np.divide(tower.f0 - f0_model_emp, tower.f0)[0])
		# ~~~~ Physics based formulation of f0 ~~~~
		f0_model_pb = tower.f0_phy(model_pb, 45)
		if f0_model_pb != None:
			print('Physics based model: %s, frequency: %s Hz' % (model_pb, f0_model_pb))
			print('Relative error between exp. f0 and emp. f0: %s' % np.divide(tower.f0 - f0_model_pb, tower.f0))
		print(' ')"""
 

	# TEST 2: Regression analysis
	# ~ You can create your own database or use the default one

	database = pd.read_excel(path_data) #default database
	database.fillna(-1, inplace=True)#fill Nan

	# Test the database excluding outliers
	#database = database[database['f0'] < 8]


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~             MODEL 1 --> 6                  ~~
	# ~~      Table 2 in Montabert et al:, 2023     ~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~ Build the Regression object~
	Reg1 = Regression(database)

	print('Regression analysis to update first empirical model (model_1 --> model_6)')
	pars, cov, ste, res, R2_tower = Reg1.fit_power_model("model_1")
	print('~~ Results ~~')
	print('a1 : ' , np.round(pars[0], 3), ' std : ', np.round(ste[0], 3) )
	print('b1 : ' , np.round(pars[1], 3), ' std : ', np.round(ste[1], 3) )
	print('R2 : ', np.round(R2_tower, 2))

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~                   MODEL 7                  ~~
	# ~~      Table 2 in Montabert et al:, 2023     ~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~ Build the Regression object~
	Reg2 = Regression(database)
	# fit power law: model 7 in Table 2
	print('Regression analysis to update first empirical model (model_7)')
	pars, cov, ste, res, R2_tower = Reg2.fit_power_model("model_7")
	print('~~ Results ~~')
	print('a1 : ' , np.round(pars[0], 3), ' std : ', np.round(ste[0], 3) )
	print('bs1 : ' , np.round(pars[1], 3), ' std : ', np.round(ste[1], 3) )
	print('R2 : ', np.round(R2_tower, 2))

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~                 MODEL 8                    ~~
	# ~~      Table 2 in Montabert et al:, 2023     ~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~ Build the Regression object~
	Reg3 = Regression(database)
	# fit power law: model 8 in Table 2
	print('Regression analysis to update first empirical model (model_8)')
	pars, cov, ste, res, R2_tower = Reg3.fit_power_model("model_8")
	print('~~ Results ~~')
	print('a1 : ' , np.round(pars[0], 3), ' std : ', np.round(ste[0], 3) )
	print('b1 : ' ,  np.round(pars[1], 3), ' std : ', np.round(ste[1], 3) )
	print('bs1 : ' ,  np.round(pars[2], 3), ' std : ', np.round(ste[2], 3) )
	print('R2 : ', np.round(R2_tower, 2))

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~             MODEL 9 --> 10                 ~~
	# ~~      Table 2 in Montabert et al:, 2023     ~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~ Build the Regression object~
	Reg4 = Regression(database)
	# fit power law: model 9 and 10 in Table 2
	print('Regression analysis to update first empirical model (model_9 and 10)')
	pars, cov, ste, res, R2_tower = Reg4.fit_power_model("model_9")
	print('~~ Results ~~')
	print('a1 : ' , np.round(pars[0], 3), ' std : ', np.round(ste[0], 3) )
	print('b1 : ' ,  np.round(pars[1], 3), ' std : ', np.round(ste[1], 3) )
	print('bs1 : ' ,  np.round(pars[2], 3), ' std : ', np.round(ste[2], 3) )
	print('a1s : ' ,  np.round(pars[3], 3), ' std : ', np.round(ste[3], 3) )
	print('b2s : ' ,  np.round(pars[4], 3), ' std : ', np.round(ste[4], 3) )
	print('R2 : ', np.round(R2_tower, 2))

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~                  MODEL 11                  ~~
	# ~~      Table 2 in Montabert et al:, 2023     ~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~ Build the Regression object~
	Reg5 = Regression(database)
	# fit power law: model 11 in Table 2
	print('Regression analysis to update first empirical model (model_11)')
	pars, cov, ste, res, R2_tower = Reg5.fit_power_model("model_11")
	print('~~ Results ~~')
	print('a1 : ' , np.round(pars[0], 3), ' std : ', np.round(ste[0], 3) )
	print('bs1 : ' ,  np.round(pars[1], 3), ' std : ', np.round(ste[1], 3) )
	print('a1s : ' ,  np.round(pars[2], 3), ' std : ', np.round(ste[2], 3) )
	print('b2s : ' ,  np.round(pars[3], 3), ' std : ', np.round(ste[3], 3) )
	print('a2s : ' ,  np.round(pars[4], 3), ' std : ', np.round(ste[4], 3) )
	print('a3s : ' ,  np.round(pars[5], 3), ' std : ', np.round(ste[5], 3) )
	print('R2 : ', np.round(R2_tower, 2))

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~                  MODEL 12                  ~~
	# ~~      Table 2 in Montabert et al:, 2023     ~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~ Build the Regression object~
	Reg6 = Regression(database)
	# fit power law: model 12 in Table 2
	print('Regression analysis to update first empirical model (model_12)')
	pars, cov, ste, res, R2_tower = Reg6.fit_power_model("model_12")
	print('~~ Results ~~')
	print('a1 : ' , np.round(pars[0], 3), ' std : ', np.round(ste[0], 3) )
	print('b1 : ' ,  np.round(pars[1], 3), ' std : ', np.round(ste[1], 3) )
	print('b2e : ' ,  np.round(pars[2], 3), ' std : ', np.round(ste[2], 3) )
	print('R2 : ', np.round(R2_tower, 2))

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~                  MODEL 13                  ~~
	# ~~      Table 2 in Montabert et al:, 2023     ~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~ Build the Regression object~
	Reg7 = Regression(database)
	# fit power law: model 13 in Table 2
	print('Regression analysis to update first empirical model (model_13)')
	pars, cov, ste, res, R2_tower = Reg7.fit_power_model("model_13")
	print('~~ Results ~~')
	print('a1 : ' , np.round(pars[0], 3), ' std : ', np.round(ste[0], 3) )
	print('b1 : ' ,  np.round(pars[1], 3), ' std : ', np.round(ste[1], 3) )
	print('b1s : ' ,  np.round(pars[2], 3), ' std : ', np.round(ste[2], 3) )
	print('b2e : ' ,  np.round(pars[3], 3), ' std : ', np.round(ste[3], 3) )
	print('R2 : ', np.round(R2_tower, 2))




	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~      PLOT FIGURE FOR EMPIRICAL LAW         ~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	fig = plt.figure(1, figsize = (10, 10))
	spec = gridspec.GridSpec(figure = fig, ncols=3, nrows=3, top = 0.95, bottom = 0.06, left = 0.1, right = 0.97, wspace = 0.1, hspace = 0.3)

	# ~ prepare axis
	ax11  = fig.add_subplot(spec[0, 0])
	ax12  = fig.add_subplot(spec[0, 1])
	ax13  = fig.add_subplot(spec[0, 2])
	ax21  = fig.add_subplot(spec[1, 0])
	ax22  = fig.add_subplot(spec[1, 1], sharex = ax12)
	ax23  = fig.add_subplot(spec[1, 2], sharex = ax13)
	ax31  = fig.add_subplot(spec[2, 0], sharex = ax21)
	"""ax32  = fig.add_subplot(spec[2, 1], sharex = ax22,sharey = ax31)
	ax33  = fig.add_subplot(spec[2, 2], sharex = ax22,sharey = ax31)
	ax41  = fig.add_subplot(spec[3, 0], sharex = ax22,sharey = ax31)
	ax42  = fig.add_subplot(spec[3, 1], sharex = ax22,sharey = ax31)"""

	# ~~ default color 
	mfc = 'dimgray'
	mec = 'k'

	# ~~ model 1 --> model 6
	prediction1  = Reg1.model_1_emp(Reg1.input, *(Reg1.pars))
	ax11.plot(Reg1.database['H'], Reg1.database['f0'], linewidth = 0, marker = 'o' , markersize = 6,markerfacecolor = mfc, markeredgecolor = mec)
	ax11.plot(Reg1.database['H'], prediction1, 'r+', ms = 6.,          label = 'model 1 -> 6')

	# ~~ model 7
	prediction2  = Reg2.model_7_emp(Reg2.input, *(Reg2.pars))
	ax12.plot(Reg2.database['H'], Reg2.database['f0'], linewidth = 0, marker = 'o' , markersize = 6,markerfacecolor = mfc, markeredgecolor = mec)
	ax12.plot(Reg2.database['H'], prediction2, 'r+',         label = 'model 7')

	# ~~ model 8
	prediction3  = Reg3.model_8_emp(Reg3.input, *(Reg3.pars))
	ax13.plot(Reg3.database['H'], Reg3.database['f0'], linewidth = 0, marker = 'o' , markersize = 6,markerfacecolor = mfc, markeredgecolor = mec)
	ax13.plot(Reg3.database['H'], prediction3, 'r+',         label = 'model 8')

	# ~~ model 9 & 10
	prediction4  = Reg4.model_9_emp(Reg4.input, *(Reg4.pars))
	ax21.plot(Reg4.database['H'], Reg4.database['f0'], linewidth = 0, marker = 'o' , markersize = 6,markerfacecolor = mfc, markeredgecolor = mec)
	ax21.plot(Reg4.database['H'], prediction4, 'r+',        label = 'model 9 & 10')

	# ~~ model 11
	prediction5  = Reg5.model_11_emp(Reg5.input, *(Reg5.pars))
	ax22.plot(Reg5.database['H'], Reg5.database['f0'], linewidth = 0, marker = 'o' , markersize = 6,markerfacecolor = mfc, markeredgecolor = mec)
	ax22.plot(Reg5.database['H'], prediction4, 'r+',         label = 'model 11')

	# ~~ model 12
	prediction6  = Reg6.model_12_emp(Reg6.input, *(Reg6.pars))
	ax23.plot(Reg6.database['H'], Reg6.database['f0'], linewidth = 0, marker = 'o' , markersize = 6,markerfacecolor = mfc, markeredgecolor = mec)
	ax23.plot(Reg6.database['H'], prediction6, 'r+',         label = 'model 12')

	# ~~ model 13
	prediction7  = Reg7.model_13_emp(Reg7.input, *(Reg7.pars))
	ax31.plot(Reg7.database['H'], Reg7.database['f0'], linewidth = 0, marker = 'o' , markersize = 6,markerfacecolor = mfc, markeredgecolor = mec, label = 'Measured ' +r'$f_0$')
	ax31.plot(Reg7.database['H'], prediction7, 'r+', label = 'Computed ' +r'$f_0$')

	# ~~ Physics based model
	data = database[(database['H'] != -1) & (database['Heff'] != -1) & (database['breadth'] != -1) & (database['length'] != -1) &  (database['max_wall_thickness'] != -1) & (database['f0'] != -1) & (database['density'] != -1) & (database['E'] != -1)]

	li_H = data['H'][0:-1]
	li_f0    = data['f0'][0:-1]
	li_mod_pb1 = []
	li_mod_pb2 = []
	li_mod_pb3 = []
	li_mod_pb4 = []
	th = 90
	for ii in data.index[0:-1]:
		print("~~~~~~~~~~~~  Tower # %s  ~~~~~~~~~~~" % str(ii))
		tower = Formulation(database.iloc[ii]) # Build Empirical formulation object
		f0_model_pb1 = tower.f0_phy('model_1', theta = th)
		f0_model_pb2 = tower.f0_phy('model_2', theta = th)
		f0_model_pb3 = tower.f0_phy('model_3', theta = th)
		f0_model_pb4 = tower.f0_phy('model_4', theta = th)

		# save data
		li_mod_pb1.append(f0_model_pb1)
		li_mod_pb2.append(f0_model_pb2)
		li_mod_pb3.append(f0_model_pb3)
		li_mod_pb4.append(f0_model_pb4)

	# ~~ Compute R2 for empirical model
	SST = np.sum((li_f0 - np.mean(li_f0))**2)
	SSR_pb1 = np.sum((li_f0 - li_mod_pb1)**2)
	SSR_pb2 = np.sum((li_f0 - li_mod_pb2)**2)
	SSR_pb3 = np.sum((li_f0 - li_mod_pb3)**2)
	SSR_pb4 = np.sum((li_f0 - li_mod_pb4)**2)

	R2_pb1 = 1 - np.divide(SSR_pb1, SST)
	R2_pb2 = 1 - np.divide(SSR_pb2, SST)
	R2_pb3 = 1 - np.divide(SSR_pb3, SST)
	R2_pb4 = 1 - np.divide(SSR_pb4, SST)



	# ~~ axis process ~~
	li_R2 = [Reg1.R2, Reg2.R2, Reg3.R2, Reg4.R2, Reg5.R2, Reg6.R2, Reg7.R2, R2_pb1, R2_pb2, R2_pb3, R2_pb4]
	li_title = ['Emp. model 1', 'Emp. model 2', 'Emp. model 3', 'Emp. model 4', 'Emp. model 5', 'Emp. model 6', 'Emp. model 7']
	for aa, ax in enumerate([ax11, ax12, ax13, ax21, ax22, ax23, ax31]):
		ax_settings_plot(ax)
		ax.text(90, 6, r'$R^2 = %s$'%np.round(li_R2[aa], 2))
		ax.text(60, 12, li_title[aa], bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'), fontsize = 12)
	
	"""for aa, ax in enumerate([ax11, ax12, ax13, ax21, ax22, ax23, ax31]):
		ax.text(90, 6, r'$R^2 = %s$'%np.round(li_R2[aa], 2))"""


	for ax in [ax11, ax21, ax31]:
		ax.set_ylabel(r'$f_0\,[Hz]$')

	for ax in [ax22, ax23, ax31]:
		ax.set_xlabel('H [m]')


	# ~~ LEGEND
	ax31.legend(bbox_to_anchor=(2,0.75), frameon=False)

	# ~~ PLOT AND SAVE
	plt.savefig('Empirical_Relation.png', dpi = 500)






	#plt.show()


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~      PLOT FIGURE FOR EMPIRICAL LAW         ~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	fig = plt.figure(2, figsize = (10, 5))
	spec = gridspec.GridSpec(figure = fig, ncols=3, nrows=2, top = 0.9, bottom = 0.12, left = 0.1, right = 0.97, wspace = 0.1, hspace = 0.3)

	ax32  = fig.add_subplot(spec[0, 0], sharex = ax22,sharey = ax31)
	ax33  = fig.add_subplot(spec[0, 1], sharex = ax22,sharey = ax31)
	ax41  = fig.add_subplot(spec[0, 2], sharex = ax22,sharey = ax31)
	ax42  = fig.add_subplot(spec[1, 0], sharex = ax22,sharey = ax31)



	# ~~ Physics based model 1
	ax32.plot(li_H, li_f0, linewidth = 0, marker = 'o' , markersize = 6,markerfacecolor = mfc, markeredgecolor = mec)
	ax32.plot(li_H, li_mod_pb1, 'r+')

	# ~~ Physics based model 2
	ax33.plot(li_H, li_f0, linewidth = 0, marker = 'o' , markersize = 6,markerfacecolor = mfc, markeredgecolor = mec)
	ax33.plot(li_H, li_mod_pb2, 'r+')

	# ~~ Physics based model 3
	ax41.plot(li_H, li_f0, linewidth = 0, marker = 'o' , markersize = 6,markerfacecolor = mfc, markeredgecolor = mec)
	ax41.plot(li_H, li_mod_pb3, 'r+')

	# ~~ Physics based model 2
	ax42.plot(li_H, li_f0, linewidth = 0, marker = 'o' , markersize = 6,markerfacecolor = mfc, markeredgecolor = mec, label = 'Measured ' +r'$f_0$')
	ax42.plot(li_H, li_mod_pb4, 'r+', label = 'Computed ' +r'$f_0$')
	# ~~ axis process ~~
	li_R2 = [R2_pb1, R2_pb2, R2_pb3, R2_pb4]
	li_title = ['Phys. based model 1', 'Phys. based model 2', 'Phys. based model 3', 'Phys. based model 4']
	for aa, ax in enumerate([ax32, ax33, ax41, ax42]):
		ax_settings_plot(ax)
		ax.text(90, 6, r'$R^2 = %s$'%np.round(li_R2[aa], 2))
		if ax in [ax32, ax33, ax41, ax42]:
			ax.text(30, 12, li_title[aa], bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'), fontsize = 12)
		else:
			ax.text(60, 12, li_title[aa], bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'), fontsize = 12)
	
	"""for aa, ax in enumerate([ax11, ax12, ax13, ax21, ax22, ax23, ax31]):
		ax.text(90, 6, r'$R^2 = %s$'%np.round(li_R2[aa], 2))"""


	for ax in [ax32, ax42]:
		ax.set_ylabel(r'$f_0\,[Hz]$')

	for ax in [ax42, ax33, ax41]:
		ax.set_xlabel('H [m]')


	# ~~ LEGEND
	ax42.legend(bbox_to_anchor=(2,0.75), frameon=False)

	# ~~ PLOT AND SAVE
	plt.savefig('Physics_Based.png', dpi = 500)

	plt.show()


