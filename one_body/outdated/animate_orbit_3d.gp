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

dataDir = "./output/"
dataFile = "orbit_int_kepler_3d_leap_phys.dat"


lengthUnitName = "kpc"
lengthUnit = 1.
velUnitName = "kms^{-1}"
velUnit = 1.
velVectorScale = 0.1	# for nice arrows on plot
timeUnitName = "Gyr"
timeUnit = 0.978
timeStep = 0.001

plotRange = 1.2


# Settings
unit(s) = sprintf(" [%s]",s)

set terminal x11 0 persist title "Energy" size 600,600 font "Times-Roman,14" enhanced solid
set xlabel "Time".unit(timeUnitName)
set ylabel "E(t)/E(0)"
pl dataDir.dataFile u 1:((abs(($2+$3)/(columnhead(2)+columnhead(3))-1.))) w l not


set terminal x11 1 persist title "Orbit" size 600,600 font "Times-Roman,14" enhanced solid

set size square
set xrange [-plotRange:plotRange]
set yrange [-plotRange:plotRange]
set zrange [-plotRange:plotRange]

set xlabel "x".unit(lengthUnitName)
set ylabel "y".unit(lengthUnitName)
set zlabel "z".unit(lengthUnitName)

stats dataDir.dataFile u 4 nooutput
do for [i=0:STATS_records-1]{ \
unset label
set label sprintf("T = %5.3f %s", i*timeUnit*timeStep,timeUnitName) at screen 0.1,0.9
spl dataDir.dataFile u ($4):($0 <= i ? $6 : 1/0):($8) w l t "orbit", \
'' u ($4):($0 == i ? $6 : 1/0):($8) w p pt 7 ps 2 t "body", \
'' u ($4):($0 == i ? $6 : 1/0):($8):(velVectorScale*$5):(velVectorScale*$7):(velVectorScale*$9) w vectors t "velocity", \
'+' u (0):(0):(0) w p pt 1 ps 2 lw 2 t "potential centre (fixed)"; \
pause 0.01 \
}