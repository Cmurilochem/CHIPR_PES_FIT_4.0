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
set output "TRIAT_CHIPR.ps"

set samples 5000

################################################################################
################################################################################

set xrange [-0:15.00]
set yrange [-0.8:5.5]
set format y '%.2f'
#set multiplot 
#set xtics 1,1,10.0
#set ytics -0.6,0.1,0.7 offset 0.5,0
set ylabel "{/Helvetica-Bold=20  Energy/E_h}"
set xlabel "{/Helvetica-Bold=20  R/a_0}"
set xzeroaxis dt 2 lc rgb "dark-gray" lw 1 #lw 2 lc 6 lw 2

#set key at graph 0.95,0.9
#set key width 5

set label "{/Helvetica-Bold=20 BASIS SET}" at graph 0.4,0.6

################################################################################
################################################################################

plot "plot_basis.res" index 0 u 1:3 smooth csplines lt 1 lw 5 lc rgb "magenta" t "new Y1 basis",\
"plot_nodes.res" index 0 u 3:4 w p pt 13 ps 1.2 lc 7 t "origin of each new contracted Y1 basis",\
"plot_basis.res" index 1 u 1:3 smooth csplines lt 1 lw 5 lc rgb "green" t "new Y2 basis",\
"plot_nodes.res" index 1 u 3:4 w p pt 13 ps 1.2 lc 7 t "origin of each new contracted Y2 basis"

#"plot_basis_abinitio.res" index 0 u 1:2 w p pt 7 lc 1 t "ab initio points Y1",\
 "plot_basis.res" index 0 u 1:2 smooth csplines lt 2 lw 5 lc rgb "gray" t "old Y1 basis",\
 "plot_nodes.res" index 0 u 1:2 w p pt 7 lc 7 t "origin of each old contracted Y1 basis",\
 "plot_basis_abinitio.res" index 1 u 1:2 w p pt 7 lc 1 t "ab initio points Y2",\
"plot_basis.res" index 1 u 1:2 smooth csplines lt 2 lw 5 lc rgb "brown" t "old Y2 basis",\
"plot_nodes.res" index 1 u 1:2 w p pt 7 lc 7 t "origin of each old contracted Y2 basis",\

################################################################################
################################################################################

unset label

ASYMP = -364.51919784E+00

set xrange [1.0:15.00]
set yrange [-0.6:0.7]
set format y '%.1f'
set xtics 1,1,15.0
set ytics -0.6,0.2,0.7 offset 0.5,0
set ylabel " {/Helvetica-Bold=20  Energy/E_h}"
set xlabel " {/Helvetica-Bold=20  R_{CC-Si}/a_0}"

set label "{/Helvetica-Bold=18 C_{2}----Si Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=90 degs}" at graph 0.4, graph 0.625

set key at graph 1.01,0.97
set key vertical maxrows 5

TSKOPUT=0.008119

set arrow from 1,(-364.99163597-ASYMP) to 15,(-364.99163597-ASYMP) nohead lw 2 lc rgb "black" dt 2

set arrow from 1,((-364.99163597-ASYMP)+TSKOPUT) to 15,((-364.99163597-ASYMP)+TSKOPUT) nohead lw 2 lc rgb "red" dt 3

#set arrow from 1,(-113.79673256-ASYMP) to 10,(-113.79673256-ASYMP) nohead lw 2 lc 3 lt 4

#PLOT 90

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/C2_Si_90.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.00000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_90.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.10000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_90.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=2.30000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_90.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=2.47913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_90.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=2.72913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_90.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=2.97913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_90.txt" index 6 u 2:($4-ASYMP) w p pt 7 lc 7 t "r=3.22913152 a_0",\
     "fort.2001" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.2001" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.2001" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.2001" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.2001" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.2001" index 5 u 2:4 w l dt 1 lc 6 t "",\
     "fort.2001" index 6 u 2:4 w l dt 1 lc 7 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 C_{2}----Si Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=82.5 degs}" at graph 0.4, graph 0.625

#PLOT 82.5

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/C2_Si_82.5.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.00000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_82.5.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.10000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_82.5.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=2.30000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_82.5.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=2.47913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_82.5.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=2.72913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_82.5.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=2.97913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_82.5.txt" index 6 u 2:($4-ASYMP) w p pt 7 lc 7 t "r=3.22913152 a_0",\
     "fort.2009" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.2009" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.2009" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.2009" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.2009" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.2009" index 5 u 2:4 w l dt 1 lc 6 t "",\
     "fort.2009" index 6 u 2:4 w l dt 1 lc 7 t ""
     
#
     
################################################################################
################################################################################  

unset label
set label "{/Helvetica-Bold=18 C_{2}----Si Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=75 degs}" at graph 0.4, graph 0.625

#PLOT 75

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/C2_Si_75.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.00000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_75.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.10000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_75.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=2.30000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_75.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=2.47913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_75.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=2.72913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_75.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=2.97913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_75.txt" index 6 u 2:($4-ASYMP) w p pt 7 lc 7 t "r=3.22913152 a_0",\
     "fort.2002" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.2002" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.2002" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.2002" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.2002" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.2002" index 5 u 2:4 w l dt 1 lc 6 t "",\
     "fort.2002" index 6 u 2:4 w l dt 1 lc 7 t ""
     
################################################################################
################################################################################     
     
unset label
set label "{/Helvetica-Bold=18 C_{2}----Si Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=67.5 degs}" at graph 0.4, graph 0.625

#PLOT 67.5

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/C2_Si_67.5.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.00000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_67.5.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.10000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_67.5.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=2.30000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_67.5.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=2.47913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_67.5.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=2.72913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_67.5.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=2.97913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_67.5.txt" index 6 u 2:($4-ASYMP) w p pt 7 lc 7 t "r=3.22913152 a_0",\
     "fort.2008" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.2008" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.2008" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.2008" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.2008" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.2008" index 5 u 2:4 w l dt 1 lc 6 t "",\
     "fort.2008" index 6 u 2:4 w l dt 1 lc 7 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 C_{2}----Si Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=60 degs}" at graph 0.4, graph 0.625

#PLOT 60

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/C2_Si_60.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.00000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_60.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.10000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_60.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=2.30000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_60.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=2.47913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_60.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=2.72913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_60.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=2.97913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_60.txt" index 6 u 2:($4-ASYMP) w p pt 7 lc 7 t "r=3.22913152 a_0",\
     "fort.2003" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.2003" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.2003" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.2003" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.2003" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.2003" index 5 u 2:4 w l dt 1 lc 6 t "",\
     "fort.2003" index 6 u 2:4 w l dt 1 lc 7 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 C_{2}----Si Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=45 degs}" at graph 0.4, graph 0.625

#PLOT 45

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/C2_Si_45.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.00000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_45.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.10000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_45.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=2.30000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_45.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=2.47913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_45.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=2.72913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_45.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=2.97913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_45.txt" index 6 u 2:($4-ASYMP) w p pt 7 lc 7 t "r=3.22913152 a_0",\
     "fort.2004" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.2004" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.2004" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.2004" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.2004" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.2004" index 5 u 2:4 w l dt 1 lc 6 t "",\
     "fort.2004" index 6 u 2:4 w l dt 1 lc 7 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 C_{2}----Si Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=30 degs}" at graph 0.4, graph 0.625

#PLOT 30

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/C2_Si_30.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.00000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_30.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.10000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_30.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=2.30000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_30.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=2.47913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_30.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=2.72913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_30.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=2.97913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_30.txt" index 6 u 2:($4-ASYMP) w p pt 7 lc 7 t "r=3.22913152 a_0",\
     "fort.2005" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.2005" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.2005" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.2005" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.2005" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.2005" index 5 u 2:4 w l dt 1 lc 6 t "",\
     "fort.2005" index 6 u 2:4 w l dt 1 lc 7 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 C_{2}----Si Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=15 degs}" at graph 0.4, graph 0.625

#PLOT 15

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/C2_Si_15.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.00000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_15.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.10000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_15.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=2.30000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_15.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=2.47913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_15.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=2.72913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_15.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=2.97913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_15.txt" index 6 u 2:($4-ASYMP) w p pt 7 lc 7 t "r=3.22913152 a_0",\
     "fort.2006" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.2006" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.2006" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.2006" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.2006" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.2006" index 5 u 2:4 w l dt 1 lc 6 t "",\
     "fort.2006" index 6 u 2:4 w l dt 1 lc 7 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 C_{2}----Si Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=0 degs}" at graph 0.4, graph 0.625

#PLOT 0

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/C2_Si_0.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.00000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_0.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.10000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_0.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=2.30000000 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_0.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=2.47913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_0.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=2.72913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_0.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=2.97913152 a_0",\
     "ABINTIO_DATA_POINTS/C2_Si_0.txt" index 6 u 2:($4-ASYMP) w p pt 7 lc 7 t "r=3.22913152 a_0",\
     "fort.2007" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.2007" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.2007" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.2007" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.2007" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.2007" index 5 u 2:4 w l dt 1 lc 6 t "",\
     "fort.2007" index 6 u 2:4 w l dt 1 lc 7 t ""     

################################################################################
################################################################################

unset label
set ylabel " {/Helvetica-Bold=20  Energy/E_h}"
set xlabel " {/Helvetica-Bold=20  R_{SiC-C}/a_0}"
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=180 degs}" at graph 0.4, graph 0.625

#PLOT 180

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/SiC_C_180.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_180.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_180.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_180.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_180.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_180.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     "fort.3001" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3001" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3001" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3001" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3001" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3001" index 5 u 2:4 w l dt 1 lc 6 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=165 degs}" at graph 0.4, graph 0.625

#PLOT 165

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/SiC_C_165.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_165.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_165.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_165.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_165.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_165.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     "fort.3002" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3002" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3002" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3002" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3002" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3002" index 5 u 2:4 w l dt 1 lc 6 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=150 degs}" at graph 0.4, graph 0.625

#PLOT 150

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/SiC_C_150.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_150.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_150.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_150.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_150.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_150.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     "fort.3003" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3003" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3003" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3003" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3003" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3003" index 5 u 2:4 w l dt 1 lc 6 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=135 degs}" at graph 0.4, graph 0.625

#PLOT 135

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/SiC_C_135.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_135.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_135.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_135.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_135.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_135.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     "fort.3004" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3004" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3004" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3004" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3004" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3004" index 5 u 2:4 w l dt 1 lc 6 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=127.5 degs}" at graph 0.4, graph 0.625

#PLOT 127.5

################################################################################
################################################################################

plot "fort.3015" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3015" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3015" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3015" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3015" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3015" index 5 u 2:4 w l dt 1 lc 6 t "",\
     "ABINTIO_DATA_POINTS/SiC_C_127.5.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_127.5.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_127.5.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_127.5.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_127.5.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_127.5.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     
################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=120 degs}" at graph 0.4, graph 0.625

#PLOT 120

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/SiC_C_120.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_120.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_120.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_120.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_120.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_120.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     "fort.3005" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3005" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3005" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3005" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3005" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3005" index 5 u 2:4 w l dt 1 lc 6 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=112.5 degs}" at graph 0.4, graph 0.625

#PLOT 112.5

################################################################################
################################################################################

plot "fort.3014" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3014" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3014" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3014" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3014" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3014" index 5 u 2:4 w l dt 1 lc 6 t "",\
     "ABINTIO_DATA_POINTS/SiC_C_112.5.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_112.5.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_112.5.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_112.5.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_112.5.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_112.5.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0"     

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=105 degs}" at graph 0.4, graph 0.625

#PLOT 105

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/SiC_C_105.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_105.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_105.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_105.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_105.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_105.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     "fort.3006" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3006" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3006" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3006" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3006" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3006" index 5 u 2:4 w l dt 1 lc 6 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=90 degs}" at graph 0.4, graph 0.625

#PLOT 90

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/SiC_C_90.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_90.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_90.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_90.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_90.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_90.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     "fort.3007" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3007" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3007" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3007" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3007" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3007" index 5 u 2:4 w l dt 1 lc 6 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=75 degs}" at graph 0.4, graph 0.625

#PLOT 75

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/SiC_C_75.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_75.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_75.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_75.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_75.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_75.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     "fort.3008" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3008" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3008" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3008" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3008" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3008" index 5 u 2:4 w l dt 1 lc 6 t ""



################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=60 degs}" at graph 0.4, graph 0.625

#PLOT 60

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/SiC_C_60.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_60.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_60.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_60.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_60.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_60.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     "fort.3009" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3009" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3009" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3009" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3009" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3009" index 5 u 2:4 w l dt 1 lc 6 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=45 degs}" at graph 0.4, graph 0.625

#PLOT 45

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/SiC_C_45.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_45.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_45.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_45.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_45.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_45.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     "fort.3010" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3010" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3010" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3010" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3010" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3010" index 5 u 2:4 w l dt 1 lc 6 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=30 degs}" at graph 0.4, graph 0.625

#PLOT 30

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/SiC_C_30.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_30.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_30.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_30.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_30.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_30.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     "fort.3011" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3011" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3011" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3011" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3011" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3011" index 5 u 2:4 w l dt 1 lc 6 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=15 degs}" at graph 0.4, graph 0.625

#PLOT 15

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/SiC_C_15.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_15.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_15.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_15.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_15.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_15.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     "fort.3012" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3012" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3012" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3012" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3012" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3012" index 5 u 2:4 w l dt 1 lc 6 t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 SiC----C Channel}" at graph 0.365, graph 0.67
set label "{/Helvetica-Bold=18 {/Symbol \121}=0 degs}" at graph 0.4, graph 0.625

#PLOT 0

################################################################################
################################################################################

plot "ABINTIO_DATA_POINTS/SiC_C_0.txt" index 0 u 2:($4-ASYMP) w p pt 7 lc 1 t "r=2.80000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_0.txt" index 1 u 2:($4-ASYMP) w p pt 7 lc 2 t "r=2.96000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_0.txt" index 2 u 2:($4-ASYMP) w p pt 7 lc 3 t "r=3.26000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_0.txt" index 3 u 2:($4-ASYMP) w p pt 7 lc 4 t "r=3.56000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_0.txt" index 4 u 2:($4-ASYMP) w p pt 7 lc 5 t "r=3.86000000 a_0",\
     "ABINTIO_DATA_POINTS/SiC_C_0.txt" index 5 u 2:($4-ASYMP) w p pt 7 lc 6 t "r=4.15000000 a_0",\
     "fort.3013" index 0 u 2:4 w l dt 1 lc 1 t "",\
     "fort.3013" index 1 u 2:4 w l dt 1 lc 2 t "",\
     "fort.3013" index 2 u 2:4 w l dt 1 lc 3 t "",\
     "fort.3013" index 3 u 2:4 w l dt 1 lc 4 t "",\
     "fort.3013" index 4 u 2:4 w l dt 1 lc 5 t "",\
     "fort.3013" index 5 u 2:4 w l dt 1 lc 6 t ""     


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
set obj rect from 0, -350 to 6000, 350 behind

set boxwidth 0.5 #absolute
set style fill solid 0.5 noborder #lt -1
set style histogram clustered

################################################################################
################################################################################

plot "dev_abinitio_points.res" u 4:5 with boxes t ""

################################################################################
################################################################################

reset
set xrange [0.00:300000]
#set yrange [-0.3:0.6]
set format y '%.1f'
#set multiplot 
#set xtics 1,1,10.0
#set ytics -0.6,0.1,0.7 offset 0.5,0
set ylabel "{/Helvetica-Bold=20  rmsd/cm^{-1}}"
set xlabel "{/Helvetica-Bold=20  stratum}"
set xzeroaxis lt 1 lc 7 lw 3 #lw 2 lc 6 lw 2

set style rect fc lt -1 fs transparent solid 0.1 noborder back
set obj rect from 104713.63289864504,0 to 1000000,1000000 front behind

#set boxwidth 0.15 #absolute
set style fill transparent solid 0.25 border lt -1
#set style histogram clustered

################################################################################
################################################################################

plot "strat_rmsd.res" u 1:4 with boxes t ""

set ylabel "{/Helvetica-Bold=20  # of points}"
set xlabel "{/Helvetica-Bold=20  stratum/cm^{-1}}"
#set xrange [0.00:300000]

plot "strat_rmsd.res" u 1:2 with boxes t ""

################################################################################
################################################################################

reset

ymax=600

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


set xrange [min:max]
set yrange [0:ymax]

set xlabel "{/Helvetica-Bold=28 deviation/kJ mol^{-1}}" offset 0.0,-1.0
set ylabel "{/Helvetica-Bold=28 # of points}" offset -3.0,0
set xtics min,(max-min)/6,max nomirror font "Helvetica,22" out 
set ytics 0,300,ymax nomirror font "Helvetica,22" in
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
 
plot "dev_abinitio_points.res" u (hist(($5*cmtokjmol),width)):(1.0) smooth freq w boxes lc rgb "magenta" notitle


exit
