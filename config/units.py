'''
Author:	Thorsten Tepper Garcia
Date:	02/07/2019
'''

from astropy import units


MASS = 1. * units.M_sun
LENGTH = 1. * units.kpc
VELOCITY = (1.*units.km/units.s)
TIME = (1.*units.kpc / (units.km / units.s)).decompose().to(units.Gyr)
