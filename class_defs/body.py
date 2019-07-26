'''
Author:	Thorsten Tepper Garcia
'''

import sys
from utils import funcs
import config.phys_consts as pc

class Body():
	"""General body class representing a body"""
	def __init__(self,mass=None,pot=None,r_vec=None,v_vec=None, df = None, massevol = None):
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
			
			df - dynamical friction function

			massevol - mass evolution function (if df is switched on)

		OUTPUT:

			Instance          

		HISTORY:

			2019-07-04 - Written - TTG

		"""
		if mass is None:
			raise ValueError("mass is a required parameter of Body instance.")
		else:
			self.mass_scale = mass
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
				raise AttributeError("No potential specified for Body object.")
			else:
				self.potential = pot

			# set cumulative mass function self-consistently
			# IMPORTANT: Must include parameter checks for NFW and PITS
			pot_name = self.potential.__name__
			if pot_name == "Kepler_Pot": 
				self.mass_cum = funcs.Kepler_Mass(self.mass_scale)
			elif pot_name == "Plummer_Pot":
				_a = self.potential.__getattribute__('_a')
				_mass = self.potential.__getattribute__('_mass')
				if _mass != self.mass_scale:
					raise ValueError("Non-matching mass in Body's potential {}".format(pot_name))
				else:
					self.mass_cum = funcs.Plummer_Mass(self.mass_scale,_a)
			elif pot_name == "Hernquist_Pot":
				_a = self.potential.__getattribute__('_a')
				_mass = self.potential.__getattribute__('_mass')
				if _mass != self.mass_scale:
					raise ValueError("Non-matching mass in Body's potential {}".format(pot_name))
				else:
					self.mass_cum = funcs.Hernquist_Mass(self.mass_scale,_a)
			elif pot_name == "NFW_Pot":
				_mass = self.potential.__getattribute__('_mass')
				_rs = self.potential.__getattribute__('_rs')
				if _mass != self.mass_scale:
					raise ValueError("Non-matching mass in Body's potential {}".format(pot_name))
				else:
					self.mass_cum = funcs.NFW_Mass(_mass,_rs)
			elif pot_name == "PITS_Pot":
				_mass = self.potential.__getattribute__('_mass')
				_a = self.potential.__getattribute__('_a')
				if _mass != self.mass_scale:
					raise ValueError("Non-matching mass in Body's potential {}".format(pot_name))
				else:
					self.mass_cum = funcs.PITS_Mass(_mass,_a)
			else:
				raise ValueError("No cumulative mass function available for Body instance with pot = {}.".\
				format(pot_name))

			# set density profile and velocity dispersion self-consistently
			# for dynamical friction calculation
			# IMPORTANT: Must include parameter checks for NFW and PITS
			self.dynamical_friction = df
			if self.dynamical_friction is not None:
				if pot_name == "Plummer_Pot":
					_mass = self.potential.__getattribute__('_mass')
					_a = self.potential.__getattribute__('_a')
					if _mass != self.mass_scale:
						raise ValueError("Non-matching mass in Body's potential {}".format(pot_name))
					else:
						self.dens = funcs.Plummer_Density(self.mass_scale,_a)
						self.vel_disp = funcs.Plummer_VelDisp(self.mass_scale, _a)
				elif pot_name == "Hernquist_Pot":
					_mass = self.potential.__getattribute__('_mass')
					_a = self.potential.__getattribute__('_a')
					if _mass != self.mass_scale:
						raise ValueError("Non-matching mass in Body's potential {}".format(pot_name))
					else:
						self.dens = funcs.Hernquist_Density(self.mass_scale,_a)
						self.vel_disp = funcs.Hernquist_VelDisp(self.mass_scale,_a)
				elif pot_name == "NFW_Pot":
					_mass = self.potential.__getattribute__('_mass')
					_rs = self.potential.__getattribute__('_rs')
					self.dens = funcs.NFW_Density(_mass,_rs)
					self.vel_disp = funcs.NFW_VelDisp(_mass,_rs)
				elif pot_name == "PITS_Pot":
					_mass = self.potential.__getattribute__('_mass')
					_a = self.potential.__getattribute__('_a')
					self.dens = funcs.PITS_Density(_mass,_a)
					self.vel_disp = funcs.PITS_VelDisp(_mass,_a)
				else:
					raise \
					ValueError("No density profile / velocity dispersion available for Body instance with pot = {}.".\
						format(pot_name))

			# set mass evolution function (if not set by input parameter)
			# only relevant for dynamical friction calculation (for now)
			# IMPORTANT: mass evolution = mass loss (for now)
			if massevol is None:			# no mass evolution
				def _mass_evol(t,*r):
					return self.mass_scale
			else:								# mass evolution as set by input paramameter
				def _mass_evol(t,*r):	# for now: mass evolution = mass loss
					mb = massevol(t,*r)
					if mb < self.mass_scale:
						self.mass_scale = mb
					return self.mass_scale
			self.mass_evol = _mass_evol



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

