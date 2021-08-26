set term png size 1200,800 font ",20" linewidth 2

if ( id==1 ) set o '../figs/m1m2_'.compdir.'.png'
if ( id==2 ) set o '../figs/ur21_'.compdir.'.png'
if ( id==3 ) set o '../figs/Phi21_'.compdir.'.png'
if ( id==4 ) set o '../figs/lin_dLOD_'.compdir.'.png'
if ( id==5 ) set o '../figs/ur20_'.compdir.'.png'
if ( id==6 ) set o '../figs/Phi20_'.compdir.'.png'
if ( id==7 ) set o '../figs/id7_'.compdir.'.png'
if ( id==8 ) set o '../figs/3Dplot_'.compdir.'.png'
if ( id==9 ) set o '../figs/full_dLOD_'.compdir.'.png'
if ( id==10 ) set o '../figs/Colatitude_'.compdir.'.png'
if ( id==11 ) set o '../figs/Longitude_'.compdir.'.png'
if ( id==12 ) set o '../figs/NormAlog_'.compdir.'.png'

set key right center
set style line 1 lw 3
set xlabel 'time [Myr]'

sk=1
My=1.e6
ky=1.e3
#set xrange[0:0.07]

file='../run/icetest.dat'
fc='../../'.compdir.'/run/icetest.dat'
file2='../run/tpw_process.dat'
fc2='../../'.compdir.'/run/tpw_process.dat'
spf2='../../spada/veen70.out'
spf4='../../spada/veen1.out'

if ( id==1 ) {
 set xlabel 'time [kyr]'
 #set xrange [0:5]
 set ylabel 'omega(1:2)/omega0 [degree]'
 p file u ($1/ky):(($2)/pi*180) ev sk ti 'm1', '' u ($1/ky):(($3)/pi*180) ev sk ti 'm2',\
  fc u ($1/ky):(($2)/pi*180) ev sk w l ti 'm1 '.compdir, '' u ($1/ky):(($3)/pi*180) ev sk w l ti 'm2 '.compdir,\
  spf4 u 2:3 w l lw 2 ti 'LD m1', spf4 u 2:4 w l lw 2 ti 'LD m2'
}
if ( id==2 ) p file u ($1/My):6 ev sk ti 're(ur_2,1)', '' u ($1/My):7 ev sk ti 'im(ur_2,1)'
if ( id==3 ) p file u ($1/My):(-$8-$10) ev sk ti 're(gravity_2,1)', '' u ($1/My):(-$9-$11) ev sk ti 'im(gravity_2,1)'
if ( id==4 ) {
 set xlabel 'time [kyr]'
 p file u ($1/ky):($4*1000) ev sk ti 'dLOD [ms]', spf2 u 2:(-$3) w l lw 2 ti 'LD -dLOD',\
  fc u ($1/ky):($4*1000) ev sk w l ti 'comp dLOD [ms]'
}
if ( id==5 ) p file u ($1/My):5 ev sk ti 're(ur_2,0)'
if ( id==6 ) p file u ($1/My):(-$12-$14) ev sk ti 're(gravity_2,0)'
if ( id==7 ) p file u 2:3:1 pal z w l, '' u 16:17:1 pal z
if ( id==8 ) sp file2 u 2:3:4:1 ev 10 pal z w l, sqrt(1-x**2-y**2)

sidday = 86164.1
if ( id==9 ) {
p file2 u ($1/ky):( -1000.*sidday*(sqrt($2**2+$3**2+$4**2)-1) ) w l ti 'true dLOD (ms)',\
 fc2 u ($1/ky):( -1000.*sidday*(sqrt($2**2+$3**2+$4**2)-1) ) w l ti compdir,\
 spf2 u 2:(-$3) w l lw 2 ti 'LD -dLOD'
}
if ( id==10 ) {
 set multiplot layout 2,1
 p file2 u ($1/My):(acos($4/sqrt($2*$2+$3*$3+$4*$4))/pi*180) ev sk w l ti 'colatitude (deg)',\
  fc2 u ($1/My):(acos($4/sqrt($2*$2+$3*$3+$4*$4))/pi*180) ev sk w l ti compdir.', col'
 p file2 u ($1/My):(atan($3/$2)/pi*180) ev sk w l ti 'longitude (deg)',\
  fc2 u ($1/My):(atan($3/$2)/pi*180) ev sk w l ti compdir.', lon'
 unset multiplot
 }
if ( id==11 ) p file2 u ($1/My):(atan($3/$2)/pi*180) ev sk w l ti 'longitude (deg)'
if ( id==12 ) {
set multiplot layout 1,2
p '../run/venus.dat' u ($1/My):2 ev sk ti 'normal direction',\
 '../../'.compdir.'/run/venus.dat' u ($1/My):2 ev sk w l ti compdir.' normal'
p '../run/venus.dat' u ($1/My):3 ev sk ti 'along track',\
 '../../'.compdir.'/run/venus.dat' u ($1/My):3 ev sk w l ti compdir.' along'
unset multiplot
}

unset log x
set o 

