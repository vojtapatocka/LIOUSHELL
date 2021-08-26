set term pngcairo enhanced size 1800,800 font ',17'
set o '../figs/profiles.png'

set key right top
set style line 1 lw 3

filc='../run/gridcentres.dat'
filed='../run/gridedges.dat'

set multiplot layout 1,4
set ylabel 'depth [1000 km]'
km=1.e6

set rmargin 4.
set xtics font ",14"
set xlabel 'log eta [Pa s]'
set format x "%.0f"
set xtics 0,2.0,40
p filed u (log10($3)):($1/km) w p ti 'viscosity'
set xtics auto

set ylabel ''
set format y ''
set lmargin 0.

#set format x "%.1f"
#set xtics 0,0.2,10
set xlabel 'G / 1e10 [Pa]'
p filed u (($2)/1.e10):($1/km) w p ti 'shear modulus'
set xtics auto

set xlabel 'density [g/cm3]'
p filc u ($3/1000):($1/km) w p ti 'density'

set xlabel 'gravity [m/s2]'
p filc u ($2):($1/km) w p ti 'gravity'

unset multiplot
set o 

#unset log x

