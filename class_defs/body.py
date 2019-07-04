'''
Author:	Thorsten Tepper Garcia
'''

import sys
sys.path.insert(0,"../.")							# include top directory in module search path
from utils import funcs

class Body():
	"""General body class representing a body"""
	def __init__(self,mass=None,pot=None,r_vec=None,v_vec=None):
		"""
		NAME:

			__init__

		PURPOSE:

			Initialize a Body instance

		INPUT:

			mass - total mass of body (scalar)

			pot	- potential of body (function of r=[x,y,z])

			r_vec - initial position vector of body (list)

			v_vec - initial velocity vector of body (list)

		OUTPUT:

			Instance          

		HISTORY:

			2019-07-04 - Written - TTG

		"""
		if mass is None:
			raise ValueError("mass is a required parameter of Body object.")
		else:
			self.mass = mass
		if pot is None:
			raise ValueError("pot is a required parameter of Body object.")
		else:
			self.potential = pot
		if r_vec is None:
			raise ValueError("r_vec is a required parameter of Body object.")
		elif v_vec is None:
			raise ValueError("v_vec is a required parameter of Body object.")
		elif len(r_vec) != len(v_vec):
			raise ValueError("r_vec and v_vec have differen dimensionality.")
		else:
			self.x = r_vec[0]
			self.y = r_vec[1]
			self.z = r_vec[2]
			self.vx = v_vec[0]
			self.vy = v_vec[1]
			self.vz = v_vec[2]

		# set cumulative mass function self-consistently
		if pot.__name__ == "Kepler_Pot": 
			self.mass_cum = funcs.Kepler_Mass(mass)
		elif pot.__name__ == "Plummer_Pot":
			_a = pot.__getattribute__('_a')
			self.mass_cum = funcs.Plummer_Mass(mass,_a)
		elif pot.__name__ == "Hernquist_Pot":
			_a = pot.__getattribute__('_a')
			self.mass_cum = funcs.Hernquist_Mass(mass,_a)
		elif pot.__name__ == "NFW_Pot":
			_rho0 = pot.__getattribute__('_rho0')
			_rs = pot.__getattribute__('_rs')
			self.mass_cum = funcs.NFW_Mass(_rho0,_rs)
		else:
			raise ValueError("No cumulative mass function defined for Body object.")

	def pos(self):
		'''Absolute position vector of Body object'''
		return [self.x,self.y,self.z]

	def vel(self):
		'''Absolute velocity vector of Body object'''
		return [self.vx,self.vy,self.vz]

	def dist(self):
		'''Absolute distance of Body object'''
		return funcs.norm(*self.pos())

	def speed(self):
		'''Absolute speed of Body object'''
		return funcs.norm(*self.vel())

	def pos_rel(self,b):
		'''Relative position vector of Body object relative to another Body object'''
		return [self.x-b.x,self.y-b.y,self.z-b.z]

	def vel_rel(self,b):
		'''Relative velocity vector of Body object relative to another Body object'''
		return [self.vx-b.vx,self.vy-b.vy,self.vz-b.vz]

	def dist_rel(self,b):
		'''Relative distance of Body object relative to another Body object'''
		return funcs.norm(*self.pos_rel(b))

	def speed_rel(self,b):
		'''Relative speed of Body object relative to another Body object'''
		return funcs.norm(*self.vel_rel(b))

