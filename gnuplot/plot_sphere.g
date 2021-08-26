#to plot a yellow  sphere: 
set parametric
set isosamples 50,50
#set hidden

R = 1.   # radius of sphere
set urange [-pi/2:pi/2]
set vrange [0:2*pi]

system(" grep FOSSILB ../run/principal.dat > ../run/FossilINI.dat ")
system(" grep LOADB ../run/principal.dat > ../run/LoadINI.dat ")
system(" grep FOSSILR ../run/principal.dat > ../run/FossilROT.dat ")
system(" grep LOADR ../run/principal.dat > ../run/LoadROT.dat ")
system(" grep TOTALR ../run/principal.dat > ../run/TotalROT.dat ")
system(" grep TOTALB ../run/principal.dat > ../run/TotalINI.dat ")
system(" grep TOTALEND ../run/principal.dat > ../run/TotalEND.dat ")
system(" grep AXIS ../run/principal.dat > ../run/RotAxis.dat ")

file1 = '../run/FossilINI.dat' 
file2 = '../run/LoadINI.dat' 
file1b = '../run/FossilROT.dat' 
file2b = '../run/LoadROT.dat' 
file3b = '../run/TotalROT.dat'
file3 = '../run/TotalINI.dat'  
file4 = '../run/RotAxis.dat' 
file5 = '../run/TotalEND.dat' 

set xlabel 'x'
set ylabel 'y'
set zrange [0:1]
set object circle at 0.0,0.0,1.0 size 10
set key font "Verdana,22"

splot R*cos(u)*cos(v),R*cos(u)*sin(v),R*sin(u) w l lc rgb "yellow" ti 'sphere',\
 file1 u 1:2:3:4:5:6 w vectors lw 3 ti 'Fossil ini',\
 file2 u 1:2:3:4:5:6 w vectors lw 3 ti 'Load ini',\
 file1b u 1:2:3:4:5:6 w vectors lw 3 ti 'Fossil rot',\
 file2b u 1:2:3:4:5:6 w vectors lw 3 ti 'Load rot',\
 file3b u 1:2:3:4:5:6 w vectors lw 3 ti 'Total rot',\
 file4 u 1:2:3:4:5:6 w vectors lw 3 ti 'Rotation axis',\
 file5 u 1:2:3:4:5:6 w vectors lw 3 ti 'Total END',\
 0,0,1 w p lc rgb "black"
