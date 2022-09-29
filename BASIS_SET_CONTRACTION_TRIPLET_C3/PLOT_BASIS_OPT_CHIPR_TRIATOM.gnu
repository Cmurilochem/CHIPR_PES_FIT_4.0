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
set yrange [-0.3:0.6]
set format y '%.2f'
#set multiplot 
#set xtics 1,1,10.0
#set ytics -0.6,0.1,0.7 offset 0.5,0
set ylabel "{/Helvetica-Bold=20  Energy/E_h}"
set xlabel "{/Helvetica-Bold=20  R/a_0}"
set xzeroaxis lt 1 lc 7 lw 3 #lw 2 lc 6 lw 2

set key at graph 0.95,0.9
set key width 5

set label "{/Helvetica-Bold=20 BASIS SET}" at graph 0.4,0.6

################################################################################
################################################################################

plot "plot_basis_abinitio.res" index 0 u 1:2 w p pt 7 lc 1 t "ab initio points Y1",\
"plot_basis.res" index 0 u 1:2 smooth csplines lt 1 lw 5 lc rgb "gray" t "fitted curve for coordinate Y1",\
"plot_nodes.res" index 0 u 1:2 w p pt 7 lc 7 t "origin of each contracted basis for Y1"

# FOR THE CASE OF AB2 AND ABC MOLECULES
#,\
"plot_basis_abinitio.res" index 1 u 1:2 w p pt 9 lc 2 t "ab initio points Y2",\
"plot_basis.res" index 1 u 1:2 smooth csplines lt 1 lw 5 lc rgb "blue" t "fitted curve for coordinate Y2",\
"plot_nodes.res" index 1 u 1:2 w p pt 7 lc 7 t "origin of each contracted basis for Y2",\
"plot_basis_abinitio.res" index 2 u 1:2 w p pt 11 lc 3 t "ab initio points Y3",\
"plot_basis.res" index 2 u 1:2 smooth csplines lt 1 lw 5 lc rgb "red" t "fitted curve for coordinate Y3",\
"plot_nodes.res" index 2 u 1:2 w p pt 7 lc 7 t "origin of each contracted basis for Y3"

################################################################################
################################################################################

exit
