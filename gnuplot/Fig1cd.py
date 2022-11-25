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
plt.rcParams["figure.figsize"] = [20, 9]
fig1 = figure()
ax = fig1.add_subplot(121, projection='polar')
ax2 = fig1.add_subplot(122, projection='polar')
ax3 = ax

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
if taxplot:
    lbl = (r"$o(t), 1.5{\cdot}10^{17}$ kg, $M_\mathrm{h}/M\rightarrow\infty$")
    #im3 = ax3.scatter(lontax*np.pi/180., coltax, marker='s', s=50, c=comg, cmap='viridis', label=lbl)

lbl = (r"$\omega(t), 1.5{\cdot}10^{17}$ kg, $M_\mathrm{h}/M\rightarrow\infty$")
#im1 = ax.scatter(lonw*np.pi/180., colw, marker='.', s=50, c=comg, cmap='viridis', label=lbl)
#plt.colorbar(im1,ax=ax,format="%.0f", orientation="horizontal", label='time [My]',pad=0.06)

dlist = ["P15", "P30", "P60"]
clist = ["C0", "C1", "C2"]
tlist = ["1.5", "3.0", "6.0"]
lwl = ["8", "5", "2"]
hflist = ["$\infty$", "$\infty$", "$\infty$"]
timel = np.linspace(0,14,57)
ax2.scatter(lonmin[0]*np.pi/180., colmin[0], c='red', s=100, marker='s')
maxcolw = max(colw)
for dname in dlist:
    omgfile = "../../" + dname + "/run/omega.dat"
    ind = dlist.index(dname)
    colL = np.loadtxt(omgfile, usecols=(5))
    lonL = np.loadtxt(omgfile, usecols=(6))
    timew = np.loadtxt(omgfile, usecols=(0))
    MIAfile = "../../" + dname + "/run/MIA.dat"
    colLf = np.loadtxt(MIAfile, usecols=(4))
    lonLf = np.loadtxt(MIAfile, usecols=(5))
    lbl = tlist[ind] + r"$\cdot 10^{17}$ kg"
    lstl = '-'
    if ('hf' in dname):
        lbl += ", $M_\mathrm{h}/M =$" + hflist[ind]
        lstl='dashed'
    #else: lbl += r", $M_\mathrm{h}/M \rightarrow$" + hflist[ind]
    scolor = clist[ind] # next(ax._get_lines.color_cycle) #['color']
    im2 = ax2.plot(lonL*np.pi/180., colL, linestyle=lstl, c=scolor, label=lbl)
    im2 = ax2.plot(lonLf[2]*np.pi/180., colLf[2], 'o', c=scolor, markersize=10)
    idt = 0
    for i in range(0,len(timew)):
        if timew[i]/yscl >= timel[idt]:
            im2 = ax2.plot(lonL[i]*np.pi/180., colL[i], 'o', c='black', markersize=2)
            idt += 1
            if idt == len(timel): break
            print ("Marking: ",dname," time [Myr] ", timew[i]/yscl)
    omgtax_comp = True #(dname == "P60_hf1")
    if omgtax_comp:
        tidfile = "../../" + dname + "/run/tidax.dat"
        colctax = np.loadtxt(tidfile, usecols=(2))
        lonctax = np.loadtxt(tidfile, usecols=(3))
        colcw = np.loadtxt(omgfile, usecols=(2))
        loncw = np.loadtxt(omgfile, usecols=(3))
        maxcolw = max(maxcolw, max(colcw))
        im3 = ax3.plot(lonctax*np.pi/180., colctax, '-',c=scolor, label="$o(t)$, "+lbl, lw=lwl[ind], zorder=-1)
        im1 = ax.plot(loncw*np.pi/180., colcw, ':',c=scolor, label="$\omega(t)$, "+lbl, lw=lwl[ind], zorder=-1) 

ax2.yaxis.set_major_formatter(ticker.FuncFormatter(degsymb))        
ax2.set_title('d) bulge-frame view')
ax2.set_rmin(0.0)
ax2.set_rmax(1.1*max(colL))
ax2.set_thetamin(-90)
ax2.set_thetamax(0)
ax2.legend(fontsize = 15, loc='lower right', bbox_to_anchor=(0.8, 0.7), framealpha=1.0)
ax2.grid(True)
ax2.plot(np.linspace(0,-90,100)*np.pi/180, np.linspace(90,90,100), '--', linewidth=1, c='gray')

# POINTS (negative load means that the largest moment, corresponding to the surface depression, is positive)
neg_load = (abs(mommax[0]) > abs(mommin[0]))
if neg_load:
    lbl = ("point-load position")
    ax.scatter(lonmax[0]*np.pi/180., colmax[0], s=100, marker='o', facecolors='none', edgecolors='red', label='n_'+lbl)
else:
    lbl = ("point load") # ($[1.5,3,6]{\cdot}10^{17}$ kg)")
    ax.scatter(lonmin[0]*np.pi/180., colmin[0], c='red', s=100, marker='s', label=lbl)
lbl = ("Eq. MIA")
#ax.scatter(lonmax[1]*np.pi/180., colmax[1], edgecolors='blue', facecolors='none', s=200, marker='o', label=lbl)
if fosslit:
    lbl = (r"load+f.b., MIA, $1.5{\cdot}10^{17}$ kg, $M_\mathrm{h}/M\rightarrow\infty$")
    #ax.scatter(lonmax[3]*np.pi/180., colmax[3], edgecolors='brown', facecolors='none', s=200, marker='*', label=lbl)
if taxplot:
    lbl = ("Eq. mIA (t=0)")
    ax3.legend()
    ax3.grid(True)
    if fosslit:
        lbl = (r"load+f.b., mIA, $1.5{\cdot}10^{17}$ kg, $M_\mathrm{h}/M\rightarrow\infty$")
        #ax3.scatter(lonmin[3]*np.pi/180., colmin[3], edgecolors='brown', facecolors='none', s=200, marker='P', label=lbl)

#ax.set_thetamin(-90)
#ax.set_thetamax(180)
ax.set_rlabel_position(180)
ax.plot(np.linspace(-90,180,200)*np.pi/180, np.linspace(90,90,200), '--', linewidth=1, c='gray')
ax.set_rmin(0.0)
ax.set_rmax(1.1*maxcolw)
if twopanes and taxplot:
    ax.set_rmax(1.1*max(max(colw),max(coltax)))
ax.grid(True)
#ax.set_rticks([5, 30, 60, 90])

plt.subplots_adjust(wspace=0.05)
ax.set_title("c) body-fixed view", pad=15.0)
ax.legend(fontsize=15, loc='lower left', bbox_to_anchor=(0.05, 0.05), frameon=True, handlelength=2, framealpha=1.0)
ax2.scatter(0, 0, marker='*', s=100, clip_on=False, c='black', label='rot. axis')
ax2.scatter(0, 90,marker='+', s=100, clip_on=False, c='black', label='tidal axis')
plt.savefig("../figs/Fig1cd.png", bbox_inches='tight')
plt.clf()
