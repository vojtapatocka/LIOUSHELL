#pr=0
plot=1
set term pngcairo enhanced size 1800,1200 font ",22" linewidth 3
set o '../figs/terms_pr'.pr.'.png'
set pointsize 2
set title 'Term by term analysis (cf. Fig. A1, doi: 10.1093/gji/ggx469)'

matplotlibc = 1
if ( matplotlibc==1 ) {
    set linetype  1 lc rgb "dark-violet" lw 1
    set linetype  2 lc rgb "#009e73" lw 1
    set linetype  3 lc rgb "#56b4e9" lw 1
    set linetype  4 lc rgb "#e69f00" lw 1
    set linetype  5 lc rgb "#f0e442" lw 1
    set linetype  6 lc rgb "#0072b2" lw 1
    set linetype  7 lc rgb "#e51e10" lw 1
    set linetype  8 lc rgb "black"   lw 1
    set linetype  9 lc rgb "gray50"  lw 1
    set linetype cycle  9
}

set log x
set xrange [1:]
#if ( pr==1 ) unset log x
ky=1000.
set xlabel 'time [ky]' # offset 0,0.5
set ylabel 'Energy [J]'

set key below left font ",21" spacing 1.0 width -5 height +1 #title 'time [ky]' maxrows 4 #outside
set style line 1 lw 3
filec='../run/enip_process.dat'
if ( pr==0 ) filec='../run/enip_ini.dat'

sk=1
P3a='P: Surface traction'
E3a='E: Surf topo in ref. grav. pot.'
P3b='P: CMB traction due to CMB topo'
E3b='E: CMB topo in ref. grav. pot.'
P4='P: Selfg b.f. + CMB press. due to SG'
E4='E: Deformed body in SG potential'
P2='P: Centrifugal f. + CMB rot. press.'
P2b='P: Tidal f. + CMB tid. press.'
E2='E: Rotational energy'
E2b='E: Tidal energy'
P1el='P: Elastic power'
E1el='E: Elastic energy'
P3c='P: term (u grad(rho) g0)'
E3c='E: Energy due to grad(rho)'

model = system("grep model ../param.in | awk '{print $2}'")
mdl = model[2:2] + 0
print 'Radial structure: model ', mdl
hid(mdl) = ( mdl>1 ? 2 : 0 )

p filec u ($1/ky):2 ev sk ti P3a, '' u ($1/ky):3 ev sk w l lw 2 ti E3a,\
'' u ($1/ky):4 ev sk ti P3b, '' u ($1/ky):5 ev sk w l lw 2 ti E3b,\
'' u ($1/ky):($6+$8) ev sk ti P4, '' u ($1/ky):18 ev sk w l lw 2 ti E4,\
'' u ($1/ky):($10+$21) ev sk ti P2, '' u ($1/ky):11 ev sk w l lw 2 ti E2,\
'' u ($1/ky):($19+$22) ev sk ti P2b, '' u ($1/ky):20 ev sk w l lw 2 ti E2b,\
'' u ($1/ky):12 ev sk ti P1el, '' u ($1/ky):13 ev sk w l lw 2 ti E1el,\
'' u ($1/ky):($15) ev sk lw hid(mdl) ti P3c, '' u ($1/ky):($16) ev sk w l lw hid(mdl) ti E3c

#'' u ($1/ky):($2+$4+$6+$8+$10+$12+$14+$15) ev sk ti 'sum of all terms in P.R.','' u ($1/ky):($3+$5+$18+$11+$13+$14+$16) ev sk w l ti 'sum of all energies'

if ( plot==1 ) set o 
#unset log x
