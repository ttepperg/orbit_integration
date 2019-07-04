'''
Author:	Thorsten Tepper Garcia
'''

import sys
sys.path.insert(0,"../.")							# include top directory in module search path
from utils import funcs

class Body():
	"""General body class representing a body"""
	def __init__(self,mass=None,pot=None,r_vec=None,v_vec=None,mass_cum=None):
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

		OPTIONAL INPUT:

		mass_cum - cumulative mass function of body (function of r=[x,y,z])

		OUTPUT:

		Instance          

		HISTORY:

		2019-07-04 - Written - TTG

		TO DO:

		- add cumulative mass function corresponding to potential
		e.g. if Phi1.__name__ == "Kepler_Pot": ...

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

# NEED to do this in a consistent way; the problem is that mass functions such as
# Plummer_Mass depend on additional parameters, as does the corresponding potential
# the trick is to get the value of this parameters from the potential definition.
# 		if pot.__name__ == "Kepler_Pot": 
# 			self.mass_cum = funcs.Kepler_Mass(mass)
# 		elif pot.__name__ == "Plummer_Pot": 
# 			self.mass_cum = funcs.Plummer_Mass(mass)
# 		else:
		if mass_cum is None:
			print("\nWARNING: No cumulative mass function defined for Body object.")
			print("Assuming a point-like mass distribution.\n")
			self.mass_cum = None
		else:
			self.mass_cum = mass_cum

	def position(self):
		return [self.x,self.y,self.z]

	def velocity(self):
		return [self.vx,self.vy,self.vz]




