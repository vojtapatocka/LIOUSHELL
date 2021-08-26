reset
set term png size 1200,800 font ',17' lw 3
set o '../figs/balance_pr'.pr.'_'.compdir.'.png'

set log x 
#if ( pr==1 ) unset log x

set key left top
set style line 1 lw 3
set ylabel 'energy(t) - energy(0) [J]'

fileb='../run/urad_process.dat'
if ( pr==0 ) fileb='../run/urad_ini.dat'
fc='../../'.compdir.'/run/urad_process.dat'
if ( pr==0 ) fc='../../'.compdir.'/run/urad_ini.dat'

sk=1
ky=1.e3

if ( pr==0 ) set multiplot layout 2,1
set xlabel ''
set format x ''
#p fileb u ($1/ky):4 ev sk ti 'E: rotational', '' u ($1/ky):5 ev sk ti 'E: gravitational',\

p fileb u ($1/ky):($4+$5-Erot0-Egrav0) ev sk ti 'E: rot + grav',\
 '' u ($1/ky):($7) ev sk ti 'E: dissipation','' u ($1/ky):($6-Eel0) ev sk ti 'E: elastic',\
 '' u ($1/ky):(($4)+($5)+($6)+($7)-Erot0-Egrav0-Eel0) ev sk ti 'E: sum',\
 fc u ($1/ky):($4+$5-Erot0-Egrav0) ev sk w l ti compdir.' E: rot + grav',\
 '' u ($1/ky):($7) ev sk w l ti 'cE: dissipation','' u ($1/ky):($6-Eel0) ev sk w l ti 'cE: elastic',\
 '' u ($1/ky):(($4)+($5)+($6)+($7)-Erot0-Egrav0-Eel0) ev sk w l ti 'cE: sum',\

set ylabel 'displacement [m]'
set xlabel 'time [ky]'
set format x "%g"

if ( pr==0 ) {
  set key left bottom
  p fileb u ($1/ky):2 ti 'surface displacement', '' u ($1/ky):3 ti 'CMB displacement',\
   '' u ($1/ky):11 ti 'driving potential / g',\
   fc u ($1/ky):2 w l ti 'c: surface displacement', '' u ($1/ky):3 w l ti 'c: CMB displacement',\
   '' u ($1/ky):11 w l ti 'c: driving potential / g'
  unset multiplot
} 
# else { set ylabel 'relative error in H'
#  p '../run/tpw_process.dat' u ($1/ky):($6-1.0) ti 'H / H0 - 1' }

set o 
#unset log x
