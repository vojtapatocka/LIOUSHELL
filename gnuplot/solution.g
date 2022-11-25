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

set multiplot layout 2,1
set xlabel ''
set format x ''
#set xrange [1.e4:]
print 'solution.g is using ref values: Erot0 ', Erot0, ', Etid0 ', Etid0, ', Egrav ', Egrav0, ', Eel ', Eel0
Erefsum = Erot0 + Egrav0 + Eel0 + Etid0

#p fileb u ($1/ky):($4+$5-Erot0-Egrav0) ev sk ti 'E: rot + grav',\

if (compdir ne 'pure') {
 print 'gnuploting a comparison: ', compdir
 p fileb u ($1/ky):($4-Erot0) ev sk ti 'E: rotational', '' u ($1/ky):($5-Egrav0) ev sk ti 'E: gravitational',\
 '' u ($1/ky):($8-Etid0) ev sk ti 'E: tidal',\
 '' u ($1/ky):($7) ev sk ti 'E: dissipation','' u ($1/ky):($6-Eel0) ev sk ti 'E: elastic',\
 '' u ($1/ky):(($4)+($5)+($6)+($7)+($8)-Erefsum) ev sk ti 'E: sum', '' u ($1/ky):($9) ev sk w l ti 'E: ext. src',\
 fc u ($1/ky):($4-Erot0) ev sk ti 'cE: rotational', '' u ($1/ky):($5-Egrav0) ev sk ti 'cE: gravitational',\
  '' u ($1/ky):($8-Etid0) ev sk ti 'cE: tidal',\
  '' u ($1/ky):($7) ev sk ti 'cE: dissipation','' u ($1/ky):($6-Eel0) ev sk ti 'cE: elastic',\
  '' u ($1/ky):(($4)+($5)+($6)+($7)+($8)-Erefsum) ev sk ti 'cE: sum', '' u ($1/ky):($9) ev sk w l ti 'cE: ext. src'
} else {
 p fileb u ($1/ky):($4-Erot0) ev sk ti 'E: rotational', '' u ($1/ky):($5-Egrav0) ev sk ti 'E: gravitational',\
 '' u ($1/ky):($8-Etid0) ev sk ti 'E: tidal',\
 '' u ($1/ky):($7) ev sk ti 'E: dissipation','' u ($1/ky):($6-Eel0) ev sk ti 'E: elastic',\
 '' u ($1/ky):(($4)+($5)+($6)+($7)+($8)-Erefsum) ev sk ti 'E: sum', '' u ($1/ky):($9) ev sk w l ti 'E: ext. src'
}

set ylabel 'displacement [m]'
set xlabel 'time [ky]'
set format x "%g"
set key left bottom

if (compdir ne 'pure') {
 p fileb u ($1/ky):2 ti 'surface displacement', '' u ($1/ky):3 ti 'CMB displacement', '' u ($1/ky):13 ti 'surf (rot+sg pot) / g0',\
  fc u ($1/ky):2 w l ti 'c: surface displacement', '' u ($1/ky):3 w l ti 'c: CMB displacement', '' u ($1/ky):13 w l ti 'c: surf (rot+sg pot) / g0'
} else {
 p fileb u ($1/ky):2 ti 'surface displacement', '' u ($1/ky):3 ti 'CMB displacement', '' u ($1/ky):13 ti 'surf (rot+sg pot) / g0'
}

unset multiplot

# else { set ylabel 'relative error in H'
#  p '../run/tpw_process.dat' u ($1/ky):($6-1.0) ti 'H / H0 - 1' }

set o 
#unset log x
