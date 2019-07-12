'''
Author:	Thorsten Tepper Garcia
'''

import sys
sys.path.insert(0,"../.")							# include top directory in module search path
from utils import funcs
import config.phys_consts as pc

class Body():
	"""General body class representing a body"""
	def __init__(self,mass=None,pot=None,r_vec=None,v_vec=None, df = None):
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
			
			df = dynamical friction switch (boolean)

		OUTPUT:

			Instance          

		HISTORY:

			2019-07-04 - Written - TTG

		"""
		if mass is None:
			raise ValueError("mass is a required parameter of Body instance.")
		else:
			self.mass = mass
		if r_vec is None:
			raise ValueError("r_vec is a required parameter of Body instance.")
		elif v_vec is None:
			raise ValueError("v_vec is a required parameter of Body instance.")
		elif len(r_vec) != len(v_vec):
			raise ValueError("r_vec and v_vec have differen dimensionality.")
		else:
			self.x = r_vec[0]
			self.y = r_vec[1]
			self.z = r_vec[2]
			self.vx = v_vec[0]
			self.vy = v_vec[1]
			self.vz = v_vec[2]
			self.pos = [self.x,self.y,self.z]
			self.vel = [self.vx,self.vy,self.vz]
			self.dist = funcs.norm(*self.pos)
			self.speed = funcs.norm(*self.vel)

			if pot is None:
				print("\nWARNING: No potential specified for Body object.")
				print("Assuming a Kepler (point-like) potential.\n")
				self.potential = funcs.Kepler_Potential(mass=self.mass)
			else:
				self.potential = pot

			# set cumulative mass function self-consistently
			pot_name = self.potential.__name__
			if pot_name == "Kepler_Pot": 
				self.mass_cum = funcs.Kepler_Mass(self.mass)
			elif pot_name == "Plummer_Pot":
				_a = self.potential.__getattribute__('_a')
				_mass = self.potential.__getattribute__('_mass')
				if _mass != self.mass:
					raise ValueError("Non-matching mass in Body's potential {}".format(pot_name))
				else:
					self.mass_cum = funcs.Plummer_Mass(self.mass,_a)
			elif pot_name == "Hernquist_Pot":
				_a = self.potential.__getattribute__('_a')
				_mass = self.potential.__getattribute__('_mass')
				if _mass != self.mass:
					raise ValueError("Non-matching mass in Body's potential {}".format(pot_name))
				else:
					self.mass_cum = funcs.Hernquist_Mass(self.mass,_a)
			elif pot_name == "NFW_Pot":
				_rho0 = self.potential.__getattribute__('_rho0')
				_rs = self.potential.__getattribute__('_rs')
				self.mass_cum = funcs.NFW_Mass(_rho0,_rs)
			else:
				raise ValueError("No cumulative mass function available for Body instance with pot = {}.".\
				format(pot_name))

			# set density profile required for dynamical friction calculation
			# self-consistently
			self.dynamical_friction = df
			if self.dynamical_friction is not None:
				if pot_name == "Plummer_Pot":
					_a = self.potential.__getattribute__('_a')
					self.dens = funcs.Plummer_Density(mass,_a)
				elif pot_name == "Hernquist_Pot":
					_a = self.potential.__getattribute__('_a')
					self.dens = funcs.Hernquist_Density(mass,_a)
				elif pot_name == "NFW_Pot":
					_rho0 = self.potential.__getattribute__('_rho0')
					_rs = self.potential.__getattribute__('_rs')
					self.dens = funcs.NFW_Density(_rho0,_rs)
				else:
					raise ValueError("No density profile available for Body instance with pot = {}.".\
					format(pot_name))


	def pos_rel(self,b):
		'''Relative position vector of Body instance relative to another Body instance'''
		return [self.x-b.x,self.y-b.y,self.z-b.z]

	def vel_rel(self,b):
		'''Relative velocity vector of Body instance relative to another Body instance'''
		return [self.vx-b.vx,self.vy-b.vy,self.vz-b.vz]

	def dist_rel(self,b):
		'''Relative distance of Body instance relative to another Body instance'''
		return funcs.norm(*self.pos_rel(b))

	def speed_rel(self,b):
		'''Relative speed of Body instance relative to another Body instance'''
		return funcs.norm(*self.vel_rel(b))

