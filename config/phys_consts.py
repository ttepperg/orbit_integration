'''
Author:	Thorsten Tepper Garcia
Date:	03/07/2019
'''
# from astropy import constants
# 
# Interesting note:
# The use of Grav obtained from astropy and normalised as follows
# leads to a significant increase in the integration time, presumably
# because Grav carries a lot of structure with (it is a Quantity object)
# it and the function that depends on it is called several 1000 times in a loop.
# 
# # Gravitational constant (in astrophysical units):
# myUnit = units.Unit("kpc km^2 / Msun s^2")
# myGrav = G.to(myUnit)
# 
# # Make it dimensionless (but recall that it is still a Quantity object):
# myGrav = myGrav * (units.Msun*units.s**2/units.kpc/units.km**2)
# # Grav = myGrav

# Therefore, it is better to define it directly as a simple float:

Grav = 4.300917270069976e-06	# units: kpc km^2 / Msun s^2

# Note that the adopted units fix the units of mass, length, and
# velocity; the last two fix the unit of time. See units.py.
