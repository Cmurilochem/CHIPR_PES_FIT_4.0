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
set yrange [-0.3:0.6]
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

#plot "plot_basis_abinitio.res" u 1:2 w p pt 7 lc 1 t "ab initio points",\
"plot_basis.res" u 1:2 smooth csplines lt 2 lw 5 lc rgb "gray" t "old basis",\
"plot_nodes.res" u 1:2 w p pt 7 lc 7 t "origin of each old contracted basis",\
"plot_basis.res" u 1:3 smooth csplines lt 1 lw 5 lc rgb "gray" t "new basis",\
"plot_nodes.res" u 3:4 w p pt 13 ps 1.2 lc 7 t "origin of each new contracted basis"

################################################################################
################################################################################

reset
set multiplot
set size nosquare 1
set xtics border nomirror norotate 
set ytics border nomirror norotate 
set lmargin 2
set bmargin 2
set mxtics 2
set mytics 2
x1=0.35
y1=0.0
y2=0.10
y3=0.25
y4=0.60
a=1.5
b=2.5
set xtics 0,2
set format x '%g'
set format y '%.2f'

################################################################################
################################################################################

unset label
set zeroaxis
set origin x1,y2
set size nosquare 0.4,0.2
set format x '%g'
set format y '%g'
set ytics -1000,1000,1000
set xtics border nomirror norotate
set ytics border nomirror norotate
set xrange [0:12]
set yrange [-1000:1000]
set xlabel "{/Helvetica-Bold=20  R/a_0}" offset 0,-0.5
set ylabel "{/Helvetica-Bold=20  Error/cm^{-1}}" offset -1.5,0.3

plot "dev_abinitio_points.res" u 1:($2) with imp lt 1 lw 3 lc rgb "red" notitle

################################################################################
################################################################################

unset label
set size 0.25,0.25
set origin 0.48,0.68
unset key
unset title
set xrange [2:3.25]
set yrange [-0.25:-0.15]
set format y '%.2f'
set xtics 2,0.5,3.0 offset 0,0.5
set xtics font "Helvetica, 10"
set ytics -0.25,0.05,-0.15 offset 0.5,0
set ytics font "Helvetica, 10"
set ylabel " {/Helvetica-Bold=10.5  V^{(2)}(R)/E_h}" offset 4.1,0
set xlabel " {/Helvetica-Bold=10.5  R/a_0}" offset -0.6,1.2

plot "plot_pot.res" t '' smooth csplines lt 1 lw 3 lc rgb "black",\
"plot_nodes.res" u 3:5 t '' w p pt 7 ps 1.0 lc 7,\
"plot_abinitio.res" u 1:2 t '' w p pt 6 ps 1.0 lc 1 

################################################################################
################################################################################

set size 0.25,0.555
set origin 0.755,0.095
#unset key
unset title
set xrange [0:25.0]
set yrange [-0.30:0.30]
set format y '%.1f'
set format x '%.1f'
set xtics 0,5,25 offset 0,0.5
set xtics font "Helvetica, 10"
set ytics -0.3,0.1,0.3 offset 0.5,0
set ytics font "Helvetica, 10"
set ylabel " {/Helvetica-Bold=10.5  Error/cm^{-1}}" offset 3.3,0
set xlabel " {/Helvetica-Bold=10.5  Level/10^{3}cm^{-1}}" offset -0.6,1.2

set key at graph 1.1,1.5

NAN(x)=100000.00

plot "plot_levels.res" u ($2/1000.00):($3) with imp lt 1 lw 3 lc rgb "blue" notitle,\
"plot_levels.res" u ($2/1000.00):($3) t '' w p pt 7 ps 1.0 lc rgb "blue",\
NAN(x) w p pt 6 ps 1.0 lc 1 t "ab initio data",\
NAN(x) w p pt 7 ps 1.0 lc 3 t "spectrosc. data",\
NAN(x) w l lt 1 lw 3 lc rgb "black" t "CHIPR",\
NAN(x) w p pt 7 ps 1.0 lc 7 t "nodes"

################################################################################
################################################################################

#set label 1 "{/Helvetica-Bold=18 C_{2}(a^3{/Symbol P}_{u})}" at graph 0.65, graph 0.65

set origin x1,y3
set xtics 0,2
set ytics -0.25,0.5
set format x ''
set format y '%.2f'
set zeroaxis
#set size nosquare 0.3333,0.3333
set size nosquare 0.4,0.4
set xrange [0:12]
set yrange [-0.25:0.05]
set xtics 0,2
set ytics -0.25,0.05,0.00
set xtics font "Helvetica, 16"
set ytics font "Helvetica, 16"
set xlabel " " 
set ylabel offset 1,4
set ylabel " {/Helvetica-Bold=20  V^{(2)}(R)/E_h}" 

plot "plot_abinitio.res" u 1:2 t '' w p pt 6 ps 1.0 lc 1,\
"plot_pot.res" u 1:2 t '' smooth csplines lt 1 lw 3 lc rgb "black",\
"plot_nodes.res" u 3:5 t '' w p pt 7 ps 1.0 lc 7

################################################################################
################################################################################

unset label
#set size nosquare 0.3333,0.23
set size nosquare 0.4,0.35
set xtics 0,2
set mytics 10
set ytics 1e-2,1e1
set xtics font "Helvetica, 16"
set ytics font "Helvetica, 16"
set origin x1,y4
set logscale y
set format y '10^{%T}'
set xrange [0:12]
set yrange [0.4E-1:1E3]
set xlabel " " 
set ylabel " "
set format x ''
#set encoding iso_8859_1
#set label 1 "{/Helvetica-Bold=18 C_{2}({/Helvetica-Italic-Bold a} ^3{/Symbol P}_{{/Helvetica-Italic-Bold u}})}" at graph 0.65, graph 0.85
set encoding default 
set key right top Right
set key spacing 1.25
set key font "Helvetica,12"
set key samplen 2
set key at graph 0.98, graph 0.60

plot "plot_pot.res" w l lt 1 lw 3 lc rgb "black" notitle,\
"plot_abinitio.res" u 1:2 t '' w p pt 6 ps 1.0 lc 1,\
"plot_nodes.res" u 3:5 t '' w p pt 7 ps 1.0 lc 7 


################################################################################
################################################################################

unset multiplot
exit

