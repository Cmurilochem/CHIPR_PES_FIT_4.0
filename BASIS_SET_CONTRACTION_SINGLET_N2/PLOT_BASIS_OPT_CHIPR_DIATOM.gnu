#!/usr/local/bin/gnuplot -persist
#
#
#       G N U P L O T
#       Version 4.3 patchlevel 0
#       last modified February 2007
#       System: Linux 2.6.20-1.2944.fc6
#
#       Copyright (C) 1986 - 1993, 1998, 2004, 2007
#       Thomas Williams, Colin Kelley and many others
#
#       Type `help` to access the on-line reference manual.
#       The gnuplot FAQ is available from
#               http://www.gnuplot.info/faq/
#
#       Send comments and help requests to  <gnuplot-beta@lists.sourceforge.net>
#       Send bug reports and suggestions to <gnuplot-beta@lists.sourceforge.net>
#
# set terminal x11
# set output

set terminal postscript landscape enhanced color "Helvetica" 16
set output "DIAT_CHIPR.ps"

################################################################################
################################################################################

set xrange [-0:15.00]
set yrange [-0.4:0.6]
set format y '%.2f'
#set multiplot 
#set xtics 1,1,10.0
#set ytics -0.6,0.1,0.7 offset 0.5,0
set ylabel "{/Helvetica-Bold=20  Energy/E_h}"
set xlabel "{/Helvetica-Bold=20  R/a_0}"
set xzeroaxis lt 1 lc 7 lw 3 #lw 2 lc 6 lw 2

set key at graph 0.95,0.9

set label "{/Helvetica-Bold=20 BASIS SET}" at graph 0.4,0.6

################################################################################
################################################################################

plot "plot_basis_abinitio.res" index 0 u 1:2 w p pt 7 lc 1 t "ab initio points Y1",\
"plot_basis.res" index 0 u 1:2 smooth csplines lt 1 lw 5 lc rgb "gray" t "fitted curve for coordinate Y1",\
"plot_nodes.res" index 0 u 1:2 w p pt 7 lc 7 t "origin of each contracted basis for Y1"

################################################################################
################################################################################

exit

