'''
Author:	Thorsten Tepper Garcia
Date:	02/07/2019
'''

from astropy import units


# Unit definition
MASS = 1. * units.M_sun
LENGTH = 1. * units.kpc
VELOCITY = (1.*units.km/units.s)
TIME = (1.*units.kpc / (units.km / units.s)).decompose().to(units.Gyr)

# Unit name
TIME_UNIT_STR = str(TIME.unit)
MASS_UNIT_STR = str(MASS.unit)
LENGTH_UNIT_STR =  str(LENGTH.unit)
VEL_UNIT_STR = str(VELOCITY.unit).replace(' / ','/')
ANG_MOM_UNIT_STR = str(LENGTH.unit*VELOCITY.unit).replace(' / ','/')
ENERGY_UNIT_STR = str(VELOCITY.unit**2).replace(' / ','/')