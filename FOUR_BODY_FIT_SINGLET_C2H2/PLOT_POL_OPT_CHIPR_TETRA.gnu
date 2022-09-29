#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.6 patchlevel 3    last modified 2013-04-12 
#    	Build System: Linux x86_64
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2013
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal wxt 0
# set output

set terminal postscript landscape enhanced color "Helvetica" 16
set output "TETRA_CHIPR.ps"

set samples 5000

################################################################################
################################################################################

set xrange [-0:15.00]
set yrange [-0.2:0.4]
set format y '%.2f'
#set multiplot 
#set xtics 1,1,10.0
#set ytics -0.6,0.1,0.7 offset 0.5,0
set ylabel "{/Helvetica-Bold=20  Energy/E_h}"
set xlabel "{/Helvetica-Bold=20  R/a_0}"
set xzeroaxis lt 1 lc 7 lw 3 #lw 2 lc 6 lw 2

set key at graph 0.95,0.9
set key width 5

set label "{/Helvetica-Bold=20 BASIS SET}" at graph 0.4,0.25

################################################################################
################################################################################

#plot "plot_basis_abinitio.res" index 0 u 1:2 w p pt 7 lc 1 t "ab initio points Y1",\
"plot_basis.res" index 0 u 1:2 smooth csplines lt 1 lw 5 lc rgb "gray" t "old Y1 basis",\
"plot_nodes.res" index 0 u 1:2 w p pt 7 lc 7 t "origins of each old contracted Y1 basis",\
"plot_basis.res" index 0 u 1:3 smooth csplines lt 1 lw 5 lc rgb "green" t "new Y1 basis",\
"plot_nodes.res" index 0 u 3:4 w p pt 13 ps 1.2 lc 7 t "origin of each new contracted Y1 basis",\
"plot_basis_abinitio.res" index 1 u 1:2 w p pt 7 lc 4 t "ab initio points Y2",\
"plot_basis.res" index 1 u 1:2 smooth csplines lt 1 lw 5 lc rgb "blue" t "old Y2 basis",\
"plot_nodes.res" index 1 u 1:2 w p pt 7 lc 7 t "origin of each old contracted Y2 basis",\
"plot_basis.res" index 1 u 1:3 smooth csplines lt 1 lw 5 lc rgb "olive" t "new Y2 basis",\
"plot_nodes.res" index 1 u 3:4 w p pt 13 ps 1.2 lc 7 t "origin of each new contracted Y2 basis",\
"plot_basis_abinitio.res" index 2 u 1:2 w p pt 7 lc rgb "brown" t "ab initio points Y3",\
"plot_basis.res" index 2 u 1:2 smooth csplines lt 1 lw 5 lc rgb "cyan" t "old Y3 basis",\
"plot_nodes.res" index 2 u 1:2 w p pt 7 lc 7 t "origin of each old contracted Y3 basis",\
"plot_basis.res" index 2 u 1:3 smooth csplines lt 1 lw 5 lc rgb "orange" t "new Y3 basis",\
"plot_nodes.res" index 2 u 3:4 w p pt 13 ps 1.2 lc 7 t "origin of each new contracted Y3 basis"


################################################################################
################################################################################

reset
#set xrange [0.00]
#set yrange [-0.3:0.6]
set format y '%.1f'
#set multiplot 
#set xtics 1,1,10.0
#set ytics -0.6,0.1,0.7 offset 0.5,0
set ylabel "{/Helvetica-Bold=20  abs. deviation/cm^{-1}}"
set xlabel "{/Helvetica-Bold=20  point number}"
set xzeroaxis lt 1 lc 7 lw 3 #lw 2 lc 6 lw 2

set style rect fc lt -1 fs solid 0.1 noborder back
set obj rect from 0, -350 to 90000, 350 behind

set boxwidth 0.5 #absolute
set style fill solid 0.5 noborder #lt -1
set style histogram clustered

################################################################################
################################################################################

plot "dev_abinitio_points.res" u 7:8 with boxes t ""

################################################################################
################################################################################

reset
set xrange [0.00:250000]
#set yrange [-0.3:0.6]
set format y '%.1f'
#set multiplot 
#set xtics 1,1,10.0
#set ytics -0.6,0.1,0.7 offset 0.5,0
set ylabel "{/Helvetica-Bold=20  rmsd/cm^{-1}}"
set xlabel "{/Helvetica-Bold=20  stratum/cm^{-1}}"
set xzeroaxis lt 1 lc 7 lw 3 #lw 2 lc 6 lw 2

set style rect fc lt -1 fs transparent solid 0.1 noborder back
#set obj rect from 27809.1914,0 to 1000000,1000000 front behind

#set boxwidth 0.15 #absolute
set style fill transparent solid 0.25 border lt -1
#set style histogram clustered

################################################################################
################################################################################

plot "strat_rmsd.res" u 1:4 with boxes t ""

set ylabel "{/Helvetica-Bold=20  # of points}"
set xlabel "{/Helvetica-Bold=20  stratum/cm^{-1}}"
#set xrange [0.00:250000]

plot "strat_rmsd.res" u 1:2 with boxes t ""

################################################################################
################################################################################

reset

ymax=6500

cmtokjmol=0.0119626592E+00

n=50 #number of intervals

max=25.3 #max value

min=-25.3 #min value

width=(max-min)/n #interval width

#function used to map a value to the intervals

hist(x,width)=width*floor(x/width)+width/2.0

#set table "check.data"
#plot "dev_abinitio_points.res" u (hist(($5*cmtokjmol),width)):(1.0) smooth freq w boxes lc rgb "magenta" notitle
#unset table

#a=1#30000
#sigma=1#1.72824872032598
mean=0.0

FIT_LIMIT=1.0E-7

#gauss(x)=(a/(sqrt(2*pi*sigma**2)))*exp(-(x-mean)**2/(2*sigma**2)) # NORMAL DISTRIBUTION

#fit gauss(x) "check.data" u 1:2 via a,sigma,mean

chemacc=4.184E+00

set arrow from -chemacc+mean,0 to -chemacc+mean,ymax nohead lw 5 lt 1 lc 7 front
set arrow from chemacc+mean,0 to chemacc+mean,ymax nohead lw 5 lt 1 lc 7 front

chemacc2=2*4.184E+00

set arrow from -chemacc2+mean,0 to -chemacc2+mean,ymax nohead lw 5 lt 4 lc 7 front
set arrow from chemacc2+mean,0 to chemacc2+mean,ymax nohead lw 5 lt 4 lc 7 front

max=50 #max value

min=-50 #min value

set xrange [min:max]
set yrange [0:ymax]

set xlabel "{/Helvetica-Bold=28 deviation/kJ mol^{-1}}" offset 0.0,-1.0
set ylabel "{/Helvetica-Bold=28 # of points}" offset -3.0,0
set xtics min,(max-min)/5,max nomirror font "Helvetica,22" out 
set ytics 0,1000,ymax nomirror font "Helvetica,22" in
set format x "%.1f"
set format y "%.0f
set mxtics 2
set mytics 2

set key at graph 0.31,0.95 #0.995,0.245
set key maxrow 5
set key right center Left
set key font "Helvetica,24"
set key samplen 3
set key spacing 4
set key width 3.0

set boxwidth width*0.9
set style fill solid 0.25 #fillstyle
set style histogram clustered

set samples 1000

##################################################################################################################
##################################################################################################################
 
plot "dev_abinitio_points.res" u (hist(($8*cmtokjmol),width)):(1.0) smooth freq w boxes lc rgb "magenta" notitle

#,\
gauss(x) w l lt 8 lc rgb "navy" lw 8 t "normal dist." #"{/Symbol \116}({/Symbol \155},{/Symbol \163})" 

##################################################################################################################
##################################################################################################################

exit
