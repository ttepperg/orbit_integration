'''
Author:	Thorsten Tepper Garcia
'''

import math
import sys
sys.path.insert(0,"../.")							# include top directory in module search path
from utils import funcs
import config.phys_consts as pc

class Orbit():
	"""General Orbit class representing an orbit"""
	def __init__(self, b1 = None, b2 = None):
		"""
		NAME:

			__init__

		PURPOSE:

			Initialize an Orbit instance

		INPUT:

			b1 - a Body instance
			b2 - a Body instance

		OUTPUT:

			Instance          

		HISTORY:

			2019-07-05 - Written - TTG

		"""
		if b1 is None or b2 is None:
			raise ValueError("Two Body instances are required when initialising an Orbit object.")
		else:
			self.b1 = b1
			self.b2 = b2


	def pos(self):
		'''position of body 2 relative to body 1'''
		return self.b2.pos_rel(self.b1)


	def vel(self):
		'''velocity of body 2 relative to body 1'''
		return self.b2.vel_rel(self.b1)


	def dist(self):
		'''distance of body 2 relative to body 1'''
		return self.b2.dist_rel(self.b1)


	def speed(self):
		'''speed of body 2 relative to body 1'''
		return self.b2.speed_rel(self.b1)


	def mtot(self):
		'''total mass'''
		return self.b1.mass + self.b2.mass


	def mred(self):
		'''reduced mass'''
		return self.b1.mass*self.b2.mass/self.mtot()

	def grav_param(self):
		'''gravitational parameter'''
		_r = self.pos()
		_mtot_r = self.b1.mass_cum(*_r)+self.b2.mass_cum(*_r)
		return pc.Grav*_mtot_r


	def ang_mom_vec(self):
		'''Angular momentum of body 2 relative to body 1'''
		return funcs.cross_prod(self.pos(),self.vel())


	def ang_mom(self):
		'''Magnitude ofngular momentum of body 2 relative to body 1'''
		return funcs.norm(*self.ang_mom_vec())


	def ang_mom_vec_norm(self):
		'''Normalised angular momentum of body 2 relative to body 1'''
		_L = self.ang_mom()
		if _L > 0:
			return [ l/_L for l in self.ang_mom_vec()]
		else:
			raise ValueError("Initial conditions imply a purely radial orbit (vanishing angular momentum).")


	def v_rad(self):
		'''radial velocity'''
		_r = self.pos()
		_d = self.dist()
		_v = self.vel()
		if _d > 0:
			return funcs.dot_prod(_r,_v) / _d
		else:
			raise ValueError("Non-positive distance in v_rad()")


	def v_tan(self):
		'''tangential velocity'''
		_r = self.pos()
		_d = self.dist()
		_v = self.vel()
		if _d > 0:
			return funcs.norm(*funcs.cross_prod(_r,_v)) / _d
		else:
			raise ValueError("Non-positive distance in v_tan()")


	def energy_pot(self):
		'''relative potential energy'''
		_r = self.pos()
		return self.mred()*(self.b1.potential(*_r)+self.b2.potential(*_r))


	def energy_kin(self):
		'''relative kinetic energy'''
		_v = self.speed()
		return 0.5*self.mred()*_v**2


	def eccentricity_vec(self):
		'''calculates the specific Laplace-Runge-Lenz vector a.k.a. 'eccentricity vector'
		'specific' means it is normalised by (G Mtot)*(M_reduced)
		the direction of the LRL vector lies along the symmetry axis of the
		conic section and points from the center of force toward the periapsis,
		i.e., the point of closest approach.
		See: https://en.wikipedia.org/wiki/Laplace–Runge–Lenz_vector'''
		_vxh = funcs.cross_prod(self.vel(),self.ang_mom_vec())
		_r = self.pos()
		_d = self.dist()
		_mu = self.grav_param()
		return [a/_mu - b/_d for a, b in zip(_vxh, _r)]


	def eccentricity(self):
		'''scalar eccentricity'''
		return funcs.norm(*self.eccentricity_vec())


	def semi_latus_rec(self):
		'''semi-latus rectum'''
		return (self.ang_mom())**2/self.grav_param()


	def semimajor(self):
		'''semimajor axis; undefined for e = 1'''
		e = self.eccentricity()
		if e < 1 or e > 1:
			return self.semi_latus_rec() / (1. - e**2)
		else:
			print("\nWARNING: semimajor axis not defined for e = 1.\n")
			return float('NaN')


	def semiminor(self):
		'''semiminor axis; undefined for e = 1'''
		e = self.eccentricity()
		if e < 1 or e > 1:
			return self.semi_latus_rec() / math.sqrt(abs(1. - e**2))
		else:
			print("\nWARNING: semiminor axis not defined for e = 1.\n")
			return float('NaN')


	def apocenter(self):
		'''apocentric distance; undefined for e >= 1'''
		e = self.eccentricity()
		if e < 1:
			return	self.semi_latus_rec() / (1. - e)
		else:
			print("\nWARNING: apocentre not defined for e >= 1.\n")
			return float('NaN')


	def pericenter(self):
		'''pericentric distance'''
		return	self.semi_latus_rec() / (1. + self.eccentricity())


	def v_apo(self):
		'''speed at apocenter'''
		return self.ang_mom() / self.apocenter()


	def v_peri(self):
		'''speed at pericenter'''
		return self.ang_mom() / self.pericenter()


	def circumference(self):
		''' Returns the approximate cirumference of an orbit with 0 < e < 1,
		exact for e = 0, and undefined for e >= 1.'''
		a = self.semimajor()
		b = self.semiminor()
		if a>=b:
			# correction factor
			corr = 0.25 * ((a-b)/(a+b))**2
			# approximate circumference to order 4 accuracy;
			# higher orders expansion exist!
			return math.pi * (a+b) * (1.+corr+(corr/4.)**4)
		else:
			print("\nWARNING: orbital circumference not defined for e >= 1.\n")
			return float('NaN')


	def period(self):
		'''orbital period; undefined for e >= 1'''
		a = self.semimajor()
		if a > 0:
			return 2.*math.pi*a*math.sqrt(a/self.grav_param())
		else:
			print("\nWARNING: orbital period not defined for a < 0.\n")
			return float('NaN')


	def pericentric_freq(self):
		'''orbital frequency at pericenter; undefined for e >= 1'''
		return self.v_peri() / self.circumference()


	def orbital_plane(self, vec = None):
		'''Determine transformation that maps state vectors onto
		orbital plane; makes use of Rodrigues' rotation formula
		See: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
		The transformation is such that the transformed angular
		momentum aligns with the positive z-axis.
		If vec is not provided, the function returns the inclination
		angle (radian) and the reference unit vector.'''
		_zaxis = [0.,0.,1.]
		_l = self.ang_mom_vec_norm()
		_cos_theta = funcs.dot_prod(_zaxis,_l)
		incl = math.acos(_cos_theta)
		# order of cross producto very important!
		_k = funcs.cross_prod(_l,_zaxis)
		k_norm = funcs.norm(*_k)
		if k_norm > 0:
			k_vec = [k/k_norm for k in _k]
		else:
			k_vec = _k
		if vec is None:
			return incl, k_vec
		else:
			return funcs.rodrigues_rot(vec,k_vec,incl)


	def apsidal_angle(self):
		'''Calcualates the angle between the periapsis -- given by
		the eccentricity vector -- and the x-axis on the orbital plane.'''
		_evec = self.eccentricity_vec()
		_evecrot = self.orbital_plane(_evec)
		_enorm = funcs.norm(*_evecrot)
		if _enorm > 1.e-10:
			return math.atan2(_evecrot[1],_evecrot[0])
		else:
			return 0.


	def asc_node_vec(self):
		'''vector of the ascending node;
		the reference vector is taken to be the z-axis;
		the reference plane, the xy-plane.'''
		_k = [0.,0.,1.]
		_h = self.ang_mom_vec()
		return funcs.cross_prod(_k,_h)


	def long_asc_node(self):
		'''Calculate the longitude of the ascending node in radians.
		The reference vector is taken to be the z-axis; the reference
		plane, the xy-plane.
		See https://en.wikipedia.org/wiki/Longitude_of_the_ascending_node'''
		_vec = self.asc_node_vec()
		_norm = funcs.norm(*_vec)
		if _norm > 0:
			_n = [n/_norm for n in _vec]
			if _n[1] >= 0:
				return math.acos(_n[0])
			else:
				return 2.*math.pi - math.acos(_n[0])
		else:
# 			print("\nWARNING: longitude of ascending node undefined. Set to 0 by convention.")
			return 0.


	def arg_periapsis(self):
		'''Calculates the argument of periapsis in radians.
		The reference vector is taken to be the z-axis;
		the reference plane, the xy-plane.'''
		_n = self.asc_node_vec()
		_nnorm = funcs.norm(*_n)
		_e = self.eccentricity_vec()
		_enorm = self.eccentricity()
		if _nnorm*_enorm > 0:
			_ndote = funcs.dot_prod(_n,_e)
			if _e[2] >= 0:
				return math.acos(_ndote/_nnorm/_enorm)
			else:
				return 2.*math.pi - math.acos(_ndote/_nnorm/_enorm)
		elif _enorm > 1.e-10:
			_l = self.ang_mom_vec()
			if _l[2] >= 0:
				return math.atan2(_e[1],_e[0])
			else:
				return 2.*math.pi - math.atan2(_e[1],_e[0])
		else:
			return 0.


	def kepler_orbit_polar(self,theta):
		'''Evaluates the analytic kepler orbit:
			r(theta) = p / 1 + e cos(theta-theta0)
			at theta (given in radians).'''
		_p = self.semi_latus_rec()
		_e = self.eccentricity()
		_theta0 = self.apsidal_angle()
		return _p / (1.+ _e * math.cos(theta-_theta0))


	def kepler_orbit_cartesian(self,theta):
		'''Returns the cartesian x,y coordinates of a Kepler orbit on the orbital
		plane'''
		_r = self.kepler_orbit_polar(theta)
		return _r*math.cos(theta), _r*math.sin(theta)



