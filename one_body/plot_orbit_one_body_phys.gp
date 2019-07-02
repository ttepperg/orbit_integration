#!/opt/local/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.2 patchlevel 6    last modified 2019-01-01 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2018
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')


# Author: Thorsten Tepper-Garcia
# Date: 02/07/2019

#----------------------------------------------------------------------------------------
# Note: It is possible to plot 3D orbits projected onto 2D, by setting dim=2.

# Invoke example:

# $> gnuplot -e "dim=2; dataFile='./output/orbit_int_kepler_2d_leap_phys.dat'" plot_orbit_one_body.gp

#----------------------------------------------------------------------------------------
# User settings

if(!exists('dim')){
	dim = 3
}

if(!exists('dataFile')){
	dataFile = "./output/orbit_int_kepler_3d_leap_phys.dat"
}

lengthUnitName = "kpc"
lengthUnit = 1.
velUnitName = "kms^{-1}"
velUnit = 1.
velVectorScale = 0.1	# for nice arrows on plot
timeUnitName = "Gyr"
timeUnit = 0.978
timeFreq = 10
timeStep = 0.001
pauseStep = 0.01

#----------------------------------------------------------------------------------------
# Internal settings

unit(s) = sprintf(" [%s]",s)
max(a,b) = a>b ? a : b
length(x,y,z) = sqrt(x**2+y**2+z**2)
ang_mom(x,vx,y,vy,z,vz) = length(y*vz-z*vy,-x*vz+z*vx,x*vy-y*vx)

# Energy conservation
set terminal x11 0 persist title "Energy" size 600,600 font "Times-Roman,14" enhanced solid
set xlabel "Time".unit(timeUnitName)
set ylabel "|{/Symbol D} E(t)/E(0) |"
pl dataFile u ($1):(abs(($2+$3)/(columnhead(2)+columnhead(3))-1.)) w l not


# Angular momentum conservation
set terminal x11 1 persist title "Angular momentum" size 600,600 font "Times-Roman,14" enhanced solid
set xlabel "Time".unit(timeUnitName)
set ylabel "|{/Symbol D} h(t)/h(0)|"
pl dataFile u ($1):(abs(($4)/columnhead(4)-1.)) w l not


# Orbit
set terminal x11 2 persist title "Orbit" size 600,600 font "Times-Roman,14" enhanced solid

set size square

stats dataFile u 5 nooutput name 'X'
stats dataFile u 7 nooutput name 'Y'
Z_min = 0.
Z_max = 0.
if(dim > 2){
	stats dataFile u 9 nooutput name 'Z'
}

X = max(abs(X_min),abs(X_max))
Y = max(abs(Y_min),abs(Y_max))
Z = max(abs(Z_min),abs(Z_max))

plotRange = 1.1*max(abs(Z),max(abs(X),abs(Y)))

print "Plot range: ", plotRange


set xrange [-plotRange:plotRange]
set yrange [-plotRange:plotRange]
if(dim > 2){
	set zrange [-plotRange:plotRange]
}

set xlabel "x".unit(lengthUnitName)
set ylabel "y".unit(lengthUnitName)
if(dim > 2){
	set zlabel "z".unit(lengthUnitName)
}

if(dim > 2){

	do for [i=0:X_records-1]{ \
		unset label
		set label sprintf("T = %5.3f %s", i*timeUnit*timeStep*timeFreq,timeUnitName) \
		at screen 0.1,0.9
		spl dataFile u ($5):($0 <= i ? $7 : 1/0):($9) w l t "orbit", \
		'' u ($5):($0 == i ? $7 : 1/0):($9) w p pt 7 ps 2 t "body", \
		'' u ($5):($0 == i ? $7 : 1/0):($9):(velVectorScale*$6):(velVectorScale*$8):(velVectorScale*$10) \
		w vectors t "velocity", \
		'+' u (0):(0):(0) w p pt 1 ps 2 lw 2 t "potential centre (fixed)"; \
		pause pauseStep \
	} # do

} else {

	do for [i=0:X_records-1]{ \
		unset label
		set label sprintf("T = %5.3f %s", i*timeUnit*timeStep*timeFreq,timeUnitName) \
		at first -0.9*plotRange,0.9*plotRange
		pl dataFile u ($5):($0 <= i ? $7 : 1/0) w l t "orbit", \
		'' u ($5):($0 == i ? $7 : 1/0) w p pt 7 ps 2 t "body", \
		'' u ($5):($0 == i ? $7 : 1/0):(velVectorScale*$6):(velVectorScale*$8) w vectors t "velocity", \
		'+' u (0):(0) w p pt 1 ps 2 lw 2 t "potential centre (fixed)"; \
		pause pauseStep \
	} # do

} # if-else

#EOF