pr=1
plot=1
if ( plot==1 ) set term pngcairo enhanced size 1600,1200 font ",22"
if ( plot==1 ) set o '../figs/terms_pr'.pr.'.png'

set log x
#if ( pr==1 ) unset log x
ky=1000.
set xlabel 'time [ky]'

set key left bottom
set style line 1 lw 3
filec='../run/enip_process.dat'
if ( pr==0 ) filec='../run/enip_ini.dat'

#set multiplot layout 2,1
sk=1
P3a='P: Surface traction'
E3a='E: Surface topography in ref. grav. potential'
P3b='P: CMB traction due to CMB topography'
E3b='E: CMB topography in ref. grav. potential'
P4='P: Selfg body force plus CMB pressure due to SG'
E4='E: Deformed body in SG potential'
P2='P: Centrifugal force plus CMB rotational pressure'
E2='E: Rotational energy'
P1el='P: Rate of change of elastic energy'
E1el='E: Elastic energy'
P3c='P: term (u grad(rho) g0)'
E3c='E: Energy due to grad(rho)'

p filec u ($1/ky):2 ev sk ti P3a, '' u ($1/ky):3 ev sk w l lw 2 ti E3a,\
'' u ($1/ky):4 ev sk ti P3b, '' u ($1/ky):5 ev sk w l lw 2 ti E3b,\
'' u ($1/ky):($6+$8) ev sk ti P4, '' u ($1/ky):18 ev sk w l lw 2 ti E4,\
'' u ($1/ky):10 ev sk ti P2, '' u ($1/ky):11 ev sk w l lw 2 ti E2,\
'' u ($1/ky):12 ev sk ti P1el, '' u ($1/ky):13 ev sk w l lw 2 ti E1el,\
'' u ($1/ky):($15) ev sk ti P3c, '' u ($1/ky):($16) ev sk w l ti E3c
#'' u ($1/ky):($2+$4+$6+$8+$10+$12+$14+$15) ev sk ti 'sum of all terms in P.R.','' u ($1/ky):($3+$5+$18+$11+$13+$14+$16) ev sk w l ti 'sum of all energies'

#p filec u ($1/ky):($2+$4+$6+$8+$10+$12+$15) ev sk ti 'sum of integrated powers', '' u ($1/ky):($3+$5+$18+$11+$13+$16) ev sk w l lw 2 ti 'sum of energies',\
#'' u ($1/ky):(-$14) ev sk w l lw 2 ti '-Dissipation'

#unset multiplot

if ( plot==1 ) set o 
#unset log x
