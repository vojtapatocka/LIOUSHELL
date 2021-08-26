pr=0

set term png size 1200,800 font ',22'
if ( pr==1 ) set o '../figs/love_pr'.pr.'.png'
if ( pr==0 ) set o '../figs/love_pr'.pr.'.png'

if ( pr==0 ) set log x 
if ( pr==1 ) unset log x 
set key right top
set style line 1 lw 3
set ylabel 'love number k2t'

fileb='../run/urad_process.dat'
if ( pr==0 ) fileb='../run/urad_ini.dat'

sk=1
ky=1.
if ( pr==0 ) ky=1000.

set key right center
set xlabel 'time [ky]'
p fileb u ($1/ky):10 ev sk w l lw 3 ti 'k2t'

set o 

#unset log x

