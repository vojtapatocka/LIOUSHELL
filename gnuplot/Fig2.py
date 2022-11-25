import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pylab import figure, show, legend, xlabel
import os.path
import sys

yscl = 1.e6
MIAfile = "../run/MIA.dat"
colmax = np.loadtxt(MIAfile, usecols=(4))
lonmax = np.loadtxt(MIAfile, usecols=(5))
mommax = np.loadtxt(MIAfile, usecols=(1))
ite = np.loadtxt(MIAfile, usecols=(9))
Mlabel = np.genfromtxt(MIAfile, dtype='str', usecols=(8))
time = np.loadtxt(MIAfile, usecols=(0))
colmin = np.loadtxt(MIAfile, usecols=(6))
lonmin = np.loadtxt(MIAfile, usecols=(7))
mommin = np.loadtxt(MIAfile, usecols=(3))
mommed = np.loadtxt(MIAfile, usecols=(2))

omgfile = "../run/omega.dat"
colw = np.loadtxt(omgfile, usecols=(2))
lonw = np.loadtxt(omgfile, usecols=(3))
itw = np.loadtxt(omgfile, usecols=(8))
timew = np.loadtxt(omgfile, usecols=(0))
colL = np.loadtxt(omgfile, usecols=(5))
lonL = np.loadtxt(omgfile, usecols=(6))
momL = np.loadtxt(omgfile, usecols=(4))

def degsymb(x, y =None):
    return "{}$\degree$".format(str(int(x)))

i0 = 2
fosslit = False
if (Mlabel[3]=="FLOAD_+_F.B."):
    i0=4
    fosslit = True
twopanes = True
taxplot = False
iterations = False
MIAplot = False
fps = False

tidfile = "../run/tidax.dat"
if os.path.exists(tidfile):
    coltax = np.loadtxt(tidfile, usecols=(2))
    lontax = np.loadtxt(tidfile, usecols=(3))
    taxplot = True

stuff = sys.argv[1].split(',')
for item in stuff:
    exec(item)

plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'lines.linewidth': 4})
plt.rcParams["figure.figsize"] = [10, 9]
fig1 = figure()
ax2 = fig1.add_subplot(111, projection='polar')

# Longitude must be provided in radians, but colatitude stays in degrees
# EVOLUTIONS / ITERATIONS
maxit = max(ite)
if (maxit < 2. or not(iterations)):
    print("Showing evolutions")
    cmia = time[i0:]/yscl
    comg = timew[:]/yscl
else:
    print("Showing iterations")
    cmia = ite[i0:]/maxit
    comg = itw[:]

dlist = ["P60", "P60_hf1", "P60_hf01", "P60_hf001", "P30_hf1", "P15_hf1"]
clist = ["C2", "C8", "C5", "C6", "C1", "C0"]
tlist = ["6.0", "6.0", "6.0", "6.0", "3.0", "1.5"]
timel = np.linspace(0,14,57)
hflist = ["$\infty$", "1.0", "0.1", "0.01", "1.0", "1.0"]
ax2.scatter(lonmin[0]*np.pi/180., colmin[0], c='red', s=100, marker='s')
for dname in dlist:
    omgfile = "../../" + dname + "/run/omega.dat"
    colL = np.loadtxt(omgfile, usecols=(5))
    lonL = np.loadtxt(omgfile, usecols=(6))
    timew = np.loadtxt(omgfile, usecols=(0))
    MIAfile = "../../" + dname + "/run/MIA.dat"
    colLf = np.loadtxt(MIAfile, usecols=(4))
    lonLf = np.loadtxt(MIAfile, usecols=(5))
    lbl = tlist[dlist.index(dname)] + r"$\cdot 10^{17}$ kg"
    lstl = '-'
    if ('hf' in dname):
        lbl += ", $M_\mathrm{h}/M =$" + hflist[dlist.index(dname)]
        lstl='-'
    else: lbl += r", $M_\mathrm{h}/M \rightarrow$" + hflist[dlist.index(dname)]
    if ('60' not in dname): lstl=":"
    scolor = clist[dlist.index(dname)] # next(ax._get_lines.color_cycle) #['color']
    im2 = ax2.plot(lonL*np.pi/180., colL, linestyle=lstl, clip_on=False, lw=5, c=scolor, label=lbl)
    im2 = ax2.plot(lonLf[2]*np.pi/180., colLf[2], 'o', clip_on=False, c=scolor, markersize=10)
    idt = 0
    for i in range(0,len(timew)):
        if timew[i]/yscl >= timel[idt]:
            im2 = ax2.plot(lonL[i]*np.pi/180., colL[i], 'o', c='black', markersize=2)
            idt += 1
            if idt == len(timel): break
            print ("Marking: ",dname," time [Myr] ", timew[i]/yscl)

ax2.yaxis.set_major_formatter(ticker.FuncFormatter(degsymb))        
#ax2.set_title('d) bulge-frame view')
ax2.set_rmin(0.0)
ax2.set_rmax(1.1*max(colL))
ax2.set_thetamin(-90)
ax2.set_thetamax(0)
ax2.legend(fontsize = 15, loc='lower right', bbox_to_anchor=(1.1, -0.08), framealpha=1.0)
ax2.grid(True)
ax2.plot(np.linspace(0,-90,100)*np.pi/180, np.linspace(90,90,100), '--', linewidth=1, c='gray')
ax2.scatter(0, 0, marker='*', s=100, clip_on=False, c='black', label='rot. axis', zorder=10)
ax2.scatter(0, 90,marker='+', s=100, clip_on=False, c='black', label='tidal axis')
#plt.colorbar(im2,ax=ax2,format="%.0f", orientation="horizontal", label='time [My]',pad=0.5)

plt.savefig("../figs/Fig2.png", bbox_inches='tight')
plt.clf()

