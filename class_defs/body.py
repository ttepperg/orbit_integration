'''
Author:	Thorsten Tepper Garcia
'''

import sys
from utils import funcs
import config.phys_consts as pc
import math

class Body():
	"""General body class representing a body"""
	def __init__(self,mass=None,pot=None,r_vec=None,v_vec=None, df = None, massevol = None, rt = None, mmin = None):
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
			
			rt - truncation radius; required if massevol is set

			mmin - present-day mass; required if massevol is set and backwards integration is on

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
				if rt is None:
					raise ValueError("NFW model is mass-divergent and requires to set a truncation radius.")
				_mass = self.potential.__getattribute__('_mass')
				_rs = self.potential.__getattribute__('_rs')
				if _mass != self.mass_scale:
					raise ValueError("Non-matching mass in Body's potential {}".format(pot_name))
				else:
					self.mass_cum = funcs.NFW_Mass(_mass,_rs)
			elif pot_name == "PITS_Pot":
				if rt is None:
					raise ValueError("PITS model is mass-divergent and requires to set a truncation radius.")
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
			self.dynamical_friction = df
			if self.dynamical_friction is not None:
				if pot_name == "Plummer_Pot":
					_mass = self.potential.__getattribute__('_mass')
					_a = self.potential.__getattribute__('_a')
					self.dens = funcs.Plummer_Density(self.mass_scale,_a)
					self.vel_disp = funcs.Plummer_VelDisp(self.mass_scale, _a)
				elif pot_name == "Hernquist_Pot":
					_mass = self.potential.__getattribute__('_mass')
					_a = self.potential.__getattribute__('_a')
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
			# IMPORTANT:
			# mass gain and mass loss are not (yet) fully symmetric. this
			# has to do with the fact that calculating the mass increase at
			# pericentre when integrating backwards requires knowledge of the
			# earlier pericentric distance, which is impossible unless the
			# orbit is known in advance.

			if rt is not None:
				self.rtrunc = rt
			else:
				self.rtrunc = 1.e10

			self.mass_bound = self.mass_cum(self.rtrunc)

			if massevol is None:							# no mass evolution (trivial)

				def _mass_evol_x(t,r,v,s):
					return self.mass_bound

				def _mass_evol_y(t,r,v,s):
					return self.mass_bound

				def _mass_evol_z(t,r,v,s):
					return self.mass_bound

			else:											# mass evolution as set by input parameter

				self.vrsave_x = -1
				def _mass_evol_x(t,r,v,s):					# s - backwards integration switch (+/-1)
					_mb = massevol(t,r)
					_vr = funcs.dot_prod(r,v)
					_r = funcs.norm(*r)
					if _vr>0 and self.vrsave_x<0:
						_massb = _mb
					else:
						_massb = 0.
					self.vrsave_x = math.copysign(1,_vr)
					if s == 1:								# forward integration: mass evolution = mass loss
						if _mb < self.mass_bound:
							self.mass_bound = _mb
						if _massb > 0.:
							print("x(+1): {} {} {:E}".format(t,_r,self.mass_bound))
					else:									# backwards integration: mass evolution = mass 'gain'
						if t == 0:
							self.mass_bound = mmin			# = to the present-day mass
						else:
							if _massb > 0.:
								self.mass_bound = _massb
						if _massb > 0.:
							print("x(-1): {} {} {:E}".format(t,_r,self.mass_bound))
					return self.mass_bound

				self.vrsave_y = -1
				def _mass_evol_y(t,r,v,s):					# s - backwards integration switch (+/-1)
					_mb = massevol(t,r)
					_vr = funcs.dot_prod(r,v)
					_r = funcs.norm(*r)
					if _vr>0 and self.vrsave_y<0:
						_massb = _mb
					else:
						_massb = 0.
					self.vrsave_y = math.copysign(1,_vr)
					if s == 1:								# forward integration: mass evolution = mass loss
						if _mb < self.mass_bound:
							self.mass_bound = _mb
						if _massb > 0.:
							print("y(+1): {} {} {:E}".format(t,_r,self.mass_bound))
					else:									# backwards integration: mass evolution = mass 'gain'
						if t == 0:
							self.mass_bound = mmin			# = to the present-day mass
						else:
							if _massb > 0.:
								self.mass_bound = _massb
						if _massb > 0.:
							print("y(-1): {} {} {:E}".format(t,_r,self.mass_bound))
					return self.mass_bound

				self.vrsave_z = -1
				def _mass_evol_z(t,r,v,s):					# s - backwards integration switch (+/-1)
					_mb = massevol(t,r)
					_vr = funcs.dot_prod(r,v)
					_r = funcs.norm(*r)
					if _vr>0 and self.vrsave_z<0:
						_massb = _mb
					else:
						_massb = 0.
					self.vrsave_z = math.copysign(1,_vr)
					if s == 1:								# forward integration: mass evolution = mass loss
						if _mb < self.mass_bound:
							self.mass_bound = _mb
						if _massb > 0.:
							print("z(+1): {} {} {:E}".format(t,_r,self.mass_bound))
					else:									# backwards integration: mass evolution = mass 'gain'
						if t == 0:
							self.mass_bound = mmin			# = to the present-day mass
						else:
							if _massb > 0.:
								self.mass_bound = _massb
						if _massb > 0.:
							print("z(-1): {} {} {:E}".format(t,_r,self.mass_bound))
					return self.mass_bound

			self.mass_evol_x = _mass_evol_x
			self.mass_evol_y = _mass_evol_y
			self.mass_evol_z = _mass_evol_z


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

