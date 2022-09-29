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
set yrange [-0.3:1.5]
set format y '%.2f'
#set multiplot 
#set xtics 1,1,10.0
#set ytics -0.6,0.1,0.7 offset 0.5,0
set ylabel "{/Helvetica-Bold=20  Energy/E_h}"
set xlabel "{/Helvetica-Bold=20  R/a_0}"
set xzeroaxis lt 1 lc 7 lw 3 #lw 2 lc 6 lw 2

#set key at graph 0.95,0.9
#set key width 5

set label "{/Helvetica-Bold=20 BASIS SET}" at graph 0.4,0.6

################################################################################
################################################################################

plot "plot_basis_abinitio.res" index 0 u 1:2 w p pt 7 lc 1 t "ab initio points Y1",\
"plot_basis.res" index 0 u 1:2 smooth csplines lt 2 lw 5 lc rgb "gray" t "old Y1 basis",\
"plot_nodes.res" index 0 u 1:2 w p pt 7 lc 7 t "origin of each old contracted Y1 basis",\
"plot_basis.res" index 0 u 1:3 smooth csplines lt 1 lw 5 lc rgb "gray" t "new Y1 basis",\
"plot_nodes.res" index 0 u 3:4 w p pt 13 ps 1.2 lc 7 t "origin of each new contracted Y1 basis"

################################################################################
################################################################################

ASYMP = -113.36656732E+00

unset label
set xrange [-0:15.00]
set yrange [-0.6:1.0]
set format y '%.2f'

set arrow from 1,(-113.84367774-ASYMP) to 10,(-113.84367774-ASYMP) nohead lw 2 lc 1 lt 2

################################################################################
################################################################################

#plot "plot_pot.res" u 1:4 smooth csplines lt 2 lw 5 lc rgb "gray" t "",\
"plot_abinitio.res" u 1:4 w p pt 7 lc 1 t ""

################################################################################
################################################################################

unset label
set xrange [1.0:10.00]
set yrange [-0.6:0.7]
set format y '%.1f'
set xtics 1,1,10.0
set ytics -0.6,0.2,0.7 offset 0.5,0
set ylabel " {/Helvetica-Bold=20  Energy/E_h}"
set xlabel " {/Helvetica-Bold=20  R_{2}/a_0}"

set label "{/Helvetica-Bold=18 {/Symbol \320}CCC=60 degs}" at graph 0.4, graph 0.625

set key at graph 1.01,0.97
set key vertical maxrows 5

set arrow from 1,(-113.84367774-ASYMP) to 10,(-113.84367774-ASYMP) nohead lw 2 lc 1 lt 2

set arrow from 1,(-113.73778197-ASYMP) to 10,(-113.73778197-ASYMP) nohead lw 2 lc 2 lt 3

set arrow from 1,(-113.79673256-ASYMP) to 10,(-113.79673256-ASYMP) nohead lw 2 lc 3 lt 4

# PLOT 60

################################################################################
################################################################################

plot "fort.101" u 2:4 smooth csplines lt 1 lw 4 lc rgb "black" t "fit",\
"ABINTIO_DATA_POINTS/60.txt" index 0 u 2:($4-ASYMP) pt 7 lc 1 t "R_{1}=1.8 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 1 u 2:($4-ASYMP) pt 7 lc 2 t "R_{1}=1.9 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 2 u 2:($4-ASYMP) pt 7 lc 3 t "R_{1}=2.0 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 3 u 2:($4-ASYMP) pt 7 lc 4 t "R_{1}=2.1 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 4 u 2:($4-ASYMP) pt 7 lc 5 t "R_{1}=2.2 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 5 u 2:($4-ASYMP) pt 7 lc 6 t "R_{1}=2.3 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 6 u 2:($4-ASYMP) pt 7 lc 7 t "R_{1}=2.4 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 7 u 2:($4-ASYMP) pt 7 lc 8 t "R_{1}=2.5 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 8 u 2:($4-ASYMP) pt 7 lc 9 t "R_{1}=2.6 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 9 u 2:($4-ASYMP) pt 7 lc 10 t "R_{1}=2.7 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 10 u 2:($4-ASYMP) pt 7 lc 11 t "R_{1}=2.8 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 11 u 2:($4-ASYMP) pt 7 lc 12 t "R_{1}=2.9 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 12 u 2:($4-ASYMP) pt 7 lc 13 t "R_{1}=3.0 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 13 u 2:($4-ASYMP) pt 7 lc 14 t "R_{1}=3.1 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 14 u 2:($4-ASYMP) pt 7 lc 15 t "R_{1}=3.2 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 15 u 2:($4-ASYMP) pt 7 lc 16 t "R_{1}=3.3 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 16 u 2:($4-ASYMP) pt 7 lc 17 t "R_{1}=3.4 a_0",\
"ABINTIO_DATA_POINTS/60.txt" index 17 u 2:($4-ASYMP) pt 7 lc 18 t "R_{1}=3.5 a_0",\
"ABINTIO_DATA_POINTS/d3h.txt" u 1:($7-ASYMP) pt 13 ps 1.5 lc 6 t "" 
#,\
"doll-crop-1.jpg" binary filetype=jpg origin=(7.7,0.14) dx=2./1000 dy=2./4500 with rgbimage t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 {/Symbol \320}CCC=70 degs}" at graph 0.4, graph 0.625

# PLOT 70

################################################################################
################################################################################

plot "fort.102" u 2:4 smooth csplines lt 1 lw 4 lc rgb "black" t "fit",\
"ABINTIO_DATA_POINTS/70.txt" index 0 u 2:($4-ASYMP) pt 7 lc 1 t "R_{1}=1.8 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 1 u 2:($4-ASYMP) pt 7 lc 2 t "R_{1}=1.9 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 2 u 2:($4-ASYMP) pt 7 lc 3 t "R_{1}=2.0 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 3 u 2:($4-ASYMP) pt 7 lc 4 t "R_{1}=2.1 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 4 u 2:($4-ASYMP) pt 7 lc 5 t "R_{1}=2.2 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 5 u 2:($4-ASYMP) pt 7 lc 6 t "R_{1}=2.3 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 6 u 2:($4-ASYMP) pt 7 lc 7 t "R_{1}=2.4 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 7 u 2:($4-ASYMP) pt 7 lc 8 t "R_{1}=2.5 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 8 u 2:($4-ASYMP) pt 7 lc 9 t "R_{1}=2.6 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 9 u 2:($4-ASYMP) pt 7 lc 10 t "R_{1}=2.7 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 10 u 2:($4-ASYMP) pt 7 lc 11 t "R_{1}=2.8 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 11 u 2:($4-ASYMP) pt 7 lc 12 t "R_{1}=2.9 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 12 u 2:($4-ASYMP) pt 7 lc 13 t "R_{1}=3.0 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 13 u 2:($4-ASYMP) pt 7 lc 14 t "R_{1}=3.1 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 14 u 2:($4-ASYMP) pt 7 lc 15 t "R_{1}=3.2 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 15 u 2:($4-ASYMP) pt 7 lc 16 t "R_{1}=3.3 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 16 u 2:($4-ASYMP) pt 7 lc 17 t "R_{1}=3.4 a_0",\
"ABINTIO_DATA_POINTS/70.txt" index 17 u 2:($4-ASYMP) pt 7 lc 18 t "R_{1}=3.5 a_0"
#,\
"doll-crop-1.jpg" binary filetype=jpg origin=(7.7,0.14) dx=2./1000 dy=2./4500 with rgbimage t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 {/Symbol \320}CCC=80 degs}" at graph 0.4, graph 0.625

# PLOT 80

################################################################################
################################################################################

plot "fort.103" u 2:4 smooth csplines lt 1 lw 4 lc rgb "black" t "fit",\
"ABINTIO_DATA_POINTS/80.txt" index 0 u 2:($4-ASYMP) pt 7 lc 1 t "R_{1}=1.8 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 1 u 2:($4-ASYMP) pt 7 lc 2 t "R_{1}=1.9 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 2 u 2:($4-ASYMP) pt 7 lc 3 t "R_{1}=2.0 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 3 u 2:($4-ASYMP) pt 7 lc 4 t "R_{1}=2.1 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 4 u 2:($4-ASYMP) pt 7 lc 5 t "R_{1}=2.2 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 5 u 2:($4-ASYMP) pt 7 lc 6 t "R_{1}=2.3 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 6 u 2:($4-ASYMP) pt 7 lc 7 t "R_{1}=2.4 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 7 u 2:($4-ASYMP) pt 7 lc 8 t "R_{1}=2.5 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 8 u 2:($4-ASYMP) pt 7 lc 9 t "R_{1}=2.6 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 9 u 2:($4-ASYMP) pt 7 lc 10 t "R_{1}=2.7 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 10 u 2:($4-ASYMP) pt 7 lc 11 t "R_{1}=2.8 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 11 u 2:($4-ASYMP) pt 7 lc 12 t "R_{1}=2.9 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 12 u 2:($4-ASYMP) pt 7 lc 13 t "R_{1}=3.0 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 13 u 2:($4-ASYMP) pt 7 lc 14 t "R_{1}=3.1 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 14 u 2:($4-ASYMP) pt 7 lc 15 t "R_{1}=3.2 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 15 u 2:($4-ASYMP) pt 7 lc 16 t "R_{1}=3.3 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 16 u 2:($4-ASYMP) pt 7 lc 17 t "R_{1}=3.4 a_0",\
"ABINTIO_DATA_POINTS/80.txt" index 17 u 2:($4-ASYMP) pt 7 lc 18 t "R_{1}=3.5 a_0"
#,\
"doll-crop-1.jpg" binary filetype=jpg origin=(7.7,0.14) dx=2./1000 dy=2./4500 with rgbimage t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 {/Symbol \320}CCC=90 degs}" at graph 0.4, graph 0.625

# PLOT 90

################################################################################
################################################################################

plot "fort.104" u 2:4 smooth csplines lt 1 lw 4 lc rgb "black" t "fit",\
"ABINTIO_DATA_POINTS/90.txt" index 0 u 2:($4-ASYMP) pt 7 lc 1 t "R_{1}=1.8 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 1 u 2:($4-ASYMP) pt 7 lc 2 t "R_{1}=1.9 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 2 u 2:($4-ASYMP) pt 7 lc 3 t "R_{1}=2.0 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 3 u 2:($4-ASYMP) pt 7 lc 4 t "R_{1}=2.1 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 4 u 2:($4-ASYMP) pt 7 lc 5 t "R_{1}=2.2 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 5 u 2:($4-ASYMP) pt 7 lc 6 t "R_{1}=2.3 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 6 u 2:($4-ASYMP) pt 7 lc 7 t "R_{1}=2.4 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 7 u 2:($4-ASYMP) pt 7 lc 8 t "R_{1}=2.5 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 8 u 2:($4-ASYMP) pt 7 lc 9 t "R_{1}=2.6 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 9 u 2:($4-ASYMP) pt 7 lc 10 t "R_{1}=2.7 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 10 u 2:($4-ASYMP) pt 7 lc 11 t "R_{1}=2.8 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 11 u 2:($4-ASYMP) pt 7 lc 12 t "R_{1}=2.9 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 12 u 2:($4-ASYMP) pt 7 lc 13 t "R_{1}=3.0 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 13 u 2:($4-ASYMP) pt 7 lc 14 t "R_{1}=3.1 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 14 u 2:($4-ASYMP) pt 7 lc 15 t "R_{1}=3.2 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 15 u 2:($4-ASYMP) pt 7 lc 16 t "R_{1}=3.3 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 16 u 2:($4-ASYMP) pt 7 lc 17 t "R_{1}=3.4 a_0",\
"ABINTIO_DATA_POINTS/90.txt" index 17 u 2:($4-ASYMP) pt 7 lc 18 t "R_{1}=3.5 a_0"
#,\
"doll-crop-1.jpg" binary filetype=jpg origin=(7.7,0.14) dx=2./1000 dy=2./4500 with rgbimage t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 {/Symbol \320}CCC=100 degs}" at graph 0.4, graph 0.625

# PLOT 100

################################################################################
################################################################################

plot "fort.105" u 2:4 smooth csplines lt 1 lw 4 lc rgb "black" t "fit",\
"ABINTIO_DATA_POINTS/100.txt" index 0 u 2:($4-ASYMP) pt 7 lc 1 t "R_{1}=1.8 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 1 u 2:($4-ASYMP) pt 7 lc 2 t "R_{1}=1.9 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 2 u 2:($4-ASYMP) pt 7 lc 3 t "R_{1}=2.0 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 3 u 2:($4-ASYMP) pt 7 lc 4 t "R_{1}=2.1 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 4 u 2:($4-ASYMP) pt 7 lc 5 t "R_{1}=2.2 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 5 u 2:($4-ASYMP) pt 7 lc 6 t "R_{1}=2.3 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 6 u 2:($4-ASYMP) pt 7 lc 7 t "R_{1}=2.4 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 7 u 2:($4-ASYMP) pt 7 lc 8 t "R_{1}=2.5 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 8 u 2:($4-ASYMP) pt 7 lc 9 t "R_{1}=2.6 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 9 u 2:($4-ASYMP) pt 7 lc 10 t "R_{1}=2.7 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 10 u 2:($4-ASYMP) pt 7 lc 11 t "R_{1}=2.8 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 11 u 2:($4-ASYMP) pt 7 lc 12 t "R_{1}=2.9 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 12 u 2:($4-ASYMP) pt 7 lc 13 t "R_{1}=3.0 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 13 u 2:($4-ASYMP) pt 7 lc 14 t "R_{1}=3.1 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 14 u 2:($4-ASYMP) pt 7 lc 15 t "R_{1}=3.2 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 15 u 2:($4-ASYMP) pt 7 lc 16 t "R_{1}=3.3 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 16 u 2:($4-ASYMP) pt 7 lc 17 t "R_{1}=3.4 a_0",\
"ABINTIO_DATA_POINTS/100.txt" index 17 u 2:($4-ASYMP) pt 7 lc 18 t "R_{1}=3.5 a_0"
#,\
"doll-crop-1.jpg" binary filetype=jpg origin=(7.7,0.14) dx=2./1000 dy=2./4500 with rgbimage t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 {/Symbol \320}CCC=110 degs}" at graph 0.4, graph 0.625

# PLOT 110

################################################################################
################################################################################

plot "fort.106" u 2:4 smooth csplines lt 1 lw 4 lc rgb "black" t "fit",\
"ABINTIO_DATA_POINTS/110.txt" index 0 u 2:($4-ASYMP) pt 7 lc 1 t "R_{1}=1.8 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 1 u 2:($4-ASYMP) pt 7 lc 2 t "R_{1}=1.9 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 2 u 2:($4-ASYMP) pt 7 lc 3 t "R_{1}=2.0 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 3 u 2:($4-ASYMP) pt 7 lc 4 t "R_{1}=2.1 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 4 u 2:($4-ASYMP) pt 7 lc 5 t "R_{1}=2.2 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 5 u 2:($4-ASYMP) pt 7 lc 6 t "R_{1}=2.3 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 6 u 2:($4-ASYMP) pt 7 lc 7 t "R_{1}=2.4 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 7 u 2:($4-ASYMP) pt 7 lc 8 t "R_{1}=2.5 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 8 u 2:($4-ASYMP) pt 7 lc 9 t "R_{1}=2.6 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 9 u 2:($4-ASYMP) pt 7 lc 10 t "R_{1}=2.7 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 10 u 2:($4-ASYMP) pt 7 lc 11 t "R_{1}=2.8 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 11 u 2:($4-ASYMP) pt 7 lc 12 t "R_{1}=2.9 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 12 u 2:($4-ASYMP) pt 7 lc 13 t "R_{1}=3.0 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 13 u 2:($4-ASYMP) pt 7 lc 14 t "R_{1}=3.1 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 14 u 2:($4-ASYMP) pt 7 lc 15 t "R_{1}=3.2 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 15 u 2:($4-ASYMP) pt 7 lc 16 t "R_{1}=3.3 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 16 u 2:($4-ASYMP) pt 7 lc 17 t "R_{1}=3.4 a_0",\
"ABINTIO_DATA_POINTS/110.txt" index 17 u 2:($4-ASYMP) pt 7 lc 18 t "R_{1}=3.5 a_0"
#,\
"doll-crop-1.jpg" binary filetype=jpg origin=(7.7,0.14) dx=2./1000 dy=2./4500 with rgbimage t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 {/Symbol \320}CCC=120 degs}" at graph 0.4, graph 0.625

# PLOT 120

################################################################################
################################################################################

plot "fort.107" u 2:4 smooth csplines lt 1 lw 4 lc rgb "black" t "fit",\
"ABINTIO_DATA_POINTS/120.txt" index 0 u 2:($4-ASYMP) pt 7 lc 1 t "R_{1}=1.8 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 1 u 2:($4-ASYMP) pt 7 lc 2 t "R_{1}=1.9 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 2 u 2:($4-ASYMP) pt 7 lc 3 t "R_{1}=2.0 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 3 u 2:($4-ASYMP) pt 7 lc 4 t "R_{1}=2.1 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 4 u 2:($4-ASYMP) pt 7 lc 5 t "R_{1}=2.2 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 5 u 2:($4-ASYMP) pt 7 lc 6 t "R_{1}=2.3 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 6 u 2:($4-ASYMP) pt 7 lc 7 t "R_{1}=2.4 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 7 u 2:($4-ASYMP) pt 7 lc 8 t "R_{1}=2.5 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 8 u 2:($4-ASYMP) pt 7 lc 9 t "R_{1}=2.6 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 9 u 2:($4-ASYMP) pt 7 lc 10 t "R_{1}=2.7 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 10 u 2:($4-ASYMP) pt 7 lc 11 t "R_{1}=2.8 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 11 u 2:($4-ASYMP) pt 7 lc 12 t "R_{1}=2.9 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 12 u 2:($4-ASYMP) pt 7 lc 13 t "R_{1}=3.0 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 13 u 2:($4-ASYMP) pt 7 lc 14 t "R_{1}=3.1 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 14 u 2:($4-ASYMP) pt 7 lc 15 t "R_{1}=3.2 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 15 u 2:($4-ASYMP) pt 7 lc 16 t "R_{1}=3.3 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 16 u 2:($4-ASYMP) pt 7 lc 17 t "R_{1}=3.4 a_0",\
"ABINTIO_DATA_POINTS/120.txt" index 17 u 2:($4-ASYMP) pt 7 lc 18 t "R_{1}=3.5 a_0"
#,\
"doll-crop-1.jpg" binary filetype=jpg origin=(7.7,0.14) dx=2./1000 dy=2./4500 with rgbimage t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 {/Symbol \320}CCC=130 degs}" at graph 0.4, graph 0.625

# PLOT 130

################################################################################
################################################################################

plot "fort.108" u 2:4 smooth csplines lt 1 lw 4 lc rgb "black" t "fit"
#,\
"doll-crop-1.jpg" binary filetype=jpg origin=(7.7,0.14) dx=2./1000 dy=2./4500 with rgbimage t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 {/Symbol \320}CCC=140 degs}" at graph 0.4, graph 0.625

# PLOT 140

################################################################################
################################################################################

plot "fort.109" u 2:4 smooth csplines lt 1 lw 4 lc rgb "black" t "fit"
#,\
"doll-crop-1.jpg" binary filetype=jpg origin=(7.7,0.14) dx=2./1000 dy=2./4500 with rgbimage t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 {/Symbol \320}CCC=150 degs}" at graph 0.4, graph 0.625

# PLOT 150

################################################################################
################################################################################

plot "fort.110" u 2:4 smooth csplines lt 1 lw 4 lc rgb "black" t "fit",\
"ABINTIO_DATA_POINTS/150.txt" index 0 u 2:($4-ASYMP) pt 7 lc 1 t "R_{1}=1.8 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 1 u 2:($4-ASYMP) pt 7 lc 2 t "R_{1}=1.9 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 2 u 2:($4-ASYMP) pt 7 lc 3 t "R_{1}=2.0 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 3 u 2:($4-ASYMP) pt 7 lc 4 t "R_{1}=2.1 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 4 u 2:($4-ASYMP) pt 7 lc 5 t "R_{1}=2.2 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 5 u 2:($4-ASYMP) pt 7 lc 6 t "R_{1}=2.3 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 6 u 2:($4-ASYMP) pt 7 lc 7 t "R_{1}=2.4 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 7 u 2:($4-ASYMP) pt 7 lc 8 t "R_{1}=2.5 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 8 u 2:($4-ASYMP) pt 7 lc 9 t "R_{1}=2.6 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 9 u 2:($4-ASYMP) pt 7 lc 10 t "R_{1}=2.7 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 10 u 2:($4-ASYMP) pt 7 lc 11 t "R_{1}=2.8 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 11 u 2:($4-ASYMP) pt 7 lc 12 t "R_{1}=2.9 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 12 u 2:($4-ASYMP) pt 7 lc 13 t "R_{1}=3.0 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 13 u 2:($4-ASYMP) pt 7 lc 14 t "R_{1}=3.1 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 14 u 2:($4-ASYMP) pt 7 lc 15 t "R_{1}=3.2 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 15 u 2:($4-ASYMP) pt 7 lc 16 t "R_{1}=3.3 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 16 u 2:($4-ASYMP) pt 7 lc 17 t "R_{1}=3.4 a_0",\
"ABINTIO_DATA_POINTS/150.txt" index 17 u 2:($4-ASYMP) pt 7 lc 18 t "R_{1}=3.5 a_0"
#,\
"doll-crop-1.jpg" binary filetype=jpg origin=(7.7,0.14) dx=2./1000 dy=2./4500 with rgbimage t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 {/Symbol \320}CCC=160 degs}" at graph 0.4, graph 0.625

# PLOT 160

################################################################################
################################################################################

plot "fort.111" u 2:4 smooth csplines lt 1 lw 4 lc rgb "black" t "fit",\
"ABINTIO_DATA_POINTS/160.txt" index 0 u 2:($4-ASYMP) pt 7 lc 1 t "R_{1}=1.8 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 1 u 2:($4-ASYMP) pt 7 lc 2 t "R_{1}=1.9 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 2 u 2:($4-ASYMP) pt 7 lc 3 t "R_{1}=2.0 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 3 u 2:($4-ASYMP) pt 7 lc 4 t "R_{1}=2.1 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 4 u 2:($4-ASYMP) pt 7 lc 5 t "R_{1}=2.2 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 5 u 2:($4-ASYMP) pt 7 lc 6 t "R_{1}=2.3 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 6 u 2:($4-ASYMP) pt 7 lc 7 t "R_{1}=2.4 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 7 u 2:($4-ASYMP) pt 7 lc 8 t "R_{1}=2.5 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 8 u 2:($4-ASYMP) pt 7 lc 9 t "R_{1}=2.6 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 9 u 2:($4-ASYMP) pt 7 lc 10 t "R_{1}=2.7 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 10 u 2:($4-ASYMP) pt 7 lc 11 t "R_{1}=2.8 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 11 u 2:($4-ASYMP) pt 7 lc 12 t "R_{1}=2.9 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 12 u 2:($4-ASYMP) pt 7 lc 13 t "R_{1}=3.0 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 13 u 2:($4-ASYMP) pt 7 lc 14 t "R_{1}=3.1 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 14 u 2:($4-ASYMP) pt 7 lc 15 t "R_{1}=3.2 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 15 u 2:($4-ASYMP) pt 7 lc 16 t "R_{1}=3.3 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 16 u 2:($4-ASYMP) pt 7 lc 17 t "R_{1}=3.4 a_0",\
"ABINTIO_DATA_POINTS/160.txt" index 17 u 2:($4-ASYMP) pt 7 lc 18 t "R_{1}=3.5 a_0"
#,\
"doll-crop-1.jpg" binary filetype=jpg origin=(7.7,0.14) dx=2./1000 dy=2./4500 with rgbimage t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 {/Symbol \320}CCC=170 degs}" at graph 0.4, graph 0.625

# PLOT 170

################################################################################
################################################################################

plot "fort.112" u 2:4 smooth csplines lt 1 lw 4 lc rgb "black" t "fit",\
"ABINTIO_DATA_POINTS/170.txt" index 0 u 2:($4-ASYMP) pt 7 lc 1 t "R_{1}=1.8 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 1 u 2:($4-ASYMP) pt 7 lc 2 t "R_{1}=1.9 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 2 u 2:($4-ASYMP) pt 7 lc 3 t "R_{1}=2.0 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 3 u 2:($4-ASYMP) pt 7 lc 4 t "R_{1}=2.1 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 4 u 2:($4-ASYMP) pt 7 lc 5 t "R_{1}=2.2 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 5 u 2:($4-ASYMP) pt 7 lc 6 t "R_{1}=2.3 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 6 u 2:($4-ASYMP) pt 7 lc 7 t "R_{1}=2.4 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 7 u 2:($4-ASYMP) pt 7 lc 8 t "R_{1}=2.5 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 8 u 2:($4-ASYMP) pt 7 lc 9 t "R_{1}=2.6 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 9 u 2:($4-ASYMP) pt 7 lc 10 t "R_{1}=2.7 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 10 u 2:($4-ASYMP) pt 7 lc 11 t "R_{1}=2.8 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 11 u 2:($4-ASYMP) pt 7 lc 12 t "R_{1}=2.9 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 12 u 2:($4-ASYMP) pt 7 lc 13 t "R_{1}=3.0 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 13 u 2:($4-ASYMP) pt 7 lc 14 t "R_{1}=3.1 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 14 u 2:($4-ASYMP) pt 7 lc 15 t "R_{1}=3.2 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 15 u 2:($4-ASYMP) pt 7 lc 16 t "R_{1}=3.3 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 16 u 2:($4-ASYMP) pt 7 lc 17 t "R_{1}=3.4 a_0",\
"ABINTIO_DATA_POINTS/170.txt" index 17 u 2:($4-ASYMP) pt 7 lc 18 t "R_{1}=3.5 a_0"
#,\
"doll-crop-1.jpg" binary filetype=jpg origin=(7.7,0.14) dx=2./1000 dy=2./4500 with rgbimage t ""

################################################################################
################################################################################

unset label
set label "{/Helvetica-Bold=18 {/Symbol \320}CCC=180 degs}" at graph 0.4, graph 0.625

# PLOT 180

################################################################################
################################################################################

plot "fort.113" u 2:4 smooth csplines lt 1 lw 4 lc rgb "black" t "fit",\
"ABINTIO_DATA_POINTS/180.txt" index 0 u 2:($4-ASYMP) pt 7 lc 1 t "R_{1}=1.8 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 1 u 2:($4-ASYMP) pt 7 lc 2 t "R_{1}=1.9 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 2 u 2:($4-ASYMP) pt 7 lc 3 t "R_{1}=2.0 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 3 u 2:($4-ASYMP) pt 7 lc 4 t "R_{1}=2.1 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 4 u 2:($4-ASYMP) pt 7 lc 5 t "R_{1}=2.2 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 5 u 2:($4-ASYMP) pt 7 lc 6 t "R_{1}=2.3 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 6 u 2:($4-ASYMP) pt 7 lc 7 t "R_{1}=2.4 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 7 u 2:($4-ASYMP) pt 7 lc 8 t "R_{1}=2.5 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 8 u 2:($4-ASYMP) pt 7 lc 9 t "R_{1}=2.6 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 9 u 2:($4-ASYMP) pt 7 lc 10 t "R_{1}=2.7 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 10 u 2:($4-ASYMP) pt 7 lc 11 t "R_{1}=2.8 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 11 u 2:($4-ASYMP) pt 7 lc 12 t "R_{1}=2.9 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 12 u 2:($4-ASYMP) pt 7 lc 13 t "R_{1}=3.0 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 13 u 2:($4-ASYMP) pt 7 lc 14 t "R_{1}=3.1 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 14 u 2:($4-ASYMP) pt 7 lc 15 t "R_{1}=3.2 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 15 u 2:($4-ASYMP) pt 7 lc 16 t "R_{1}=3.3 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 16 u 2:($4-ASYMP) pt 7 lc 17 t "R_{1}=3.4 a_0",\
"ABINTIO_DATA_POINTS/180.txt" index 17 u 2:($4-ASYMP) pt 7 lc 18 t "R_{1}=3.5 a_0",\
"ABINTIO_DATA_POINTS/ts_dinfh.txt" u 1:($7-ASYMP) pt 13 ps 1.5 lc 6 t "",\
"ABINTIO_DATA_POINTS/cinfv.txt" u 1:($7-ASYMP) pt 13 ps 1.5 lc 6 t ""
#,\
"doll-crop-1.jpg" binary filetype=jpg origin=(7.7,0.14) dx=2./1000 dy=2./4500 with rgbimage t ""

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

ymax=1800

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
