set term pngcairo enhanced size 1800,800 font ',17'
set o '../figs/profiles.png'

set key right top
set style line 1 lw 3

filc='../run/gridcentres.dat'
filed='../run/gridedges.dat'
sidday=system("grep sidday ../param.in | awk '{print $2}'")

set multiplot layout 1,4 title 'Spin period'.sidday.' s'
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

set format x "%.1f"
temperature=0
if ( temperature==1 ){ 
    set xlabel 'T [K]'
    p filed u ($4):($1/km) w p ti 'temperature'
} else {
    set xlabel 'G / 1e10 [Pa]'
    p filed u (($2)/1.e10):($1/km) w p ti 'shear modulus'
}

set xlabel 'density [g/cm3]'
#set xtics 0.05
p filc u ($3/1000):($1/km) w p ti 'mantle',\
 '< tail -n 1 '.filc u ($3/1000):($1/km) w p ps 2 pt 4 ti 'core'

set xtics auto
set xlabel 'gravity [m/s2]'
p filc u ($2):($1/km) w p ti 'gravity'

unset multiplot
set o 

#unset log x

