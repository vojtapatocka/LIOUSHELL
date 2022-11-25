import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pylab import figure, show, legend, xlabel
import os.path
import sys

def degsymb(x, y =None):
    return "{}$\degree$".format(str(int(x)))

ky = 1.e6
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
if twopanes or not(taxplot):
    plt.rcParams["figure.figsize"] = [20, 9]
    fig1 = figure()
    ax = fig1.add_subplot(121, projection='polar')
    ax2 = fig1.add_subplot(122, projection='polar')
    ax3 = ax
else:
    plt.rcParams["figure.figsize"] = [22, 9]
    fig1 = figure()
    ax = fig1.add_subplot(131, projection='polar')
    ax2 = fig1.add_subplot(132, projection='polar')
    ax3 = fig1.add_subplot(133, projection='polar')

# Longitude must be provided in radians, but colatitude stays in degrees
# EVOLUTIONS / ITERATIONS
maxit = max(ite)
if (maxit < 2. or not(iterations)):
    print("Showing evolutions")
    cmia = time[i0:]/ky
    comg = timew[:]/ky
else:
    print("Showing iterations")
    cmia = ite[i0:]/maxit
    comg = itw[:]
if MIAplot:
    #print("Min %.10e max %.10e" % (min(cmia),max(cmia)))
    ax.scatter(lonmax[i0:]*np.pi/180., colmax[i0:], s=70, c=cmia, marker='x', cmap='cool', label='MIA (cool palette)')
if taxplot:
    if MIAplot:
        im4 = ax3.scatter(lonmin[i0:]*np.pi/180., colmin[i0:], s=70, c=cmia, marker='+', cmap='cool', label='mIA')
    if fps: lbl = ("$o$, f.p.: %.1f°, %.1f°" % (coltax[-1], lontax[-1]))
    else: lbl = (r"$o(t), 1.5{\cdot}10^{17}$ kg, $M_\mathrm{h}/M\rightarrow\infty$")
    im3 = ax3.scatter(lontax*np.pi/180., coltax, marker='s', s=50, c=comg, cmap='viridis', label=lbl)

if fps: lbl = ("$\omega$, f.p.: %.1f°, %.1f°" % (colw[-1],lonw[-1]))
else: lbl = (r"$\omega(t), 1.5{\cdot}10^{17}$ kg, $M_\mathrm{h}/M\rightarrow\infty$")
im1 = ax.scatter(lonw*np.pi/180., colw, marker='.', s=50, c=comg, cmap='viridis', label=lbl)
plt.colorbar(im1,ax=ax,format="%.0f", orientation="horizontal", label='time [My]',pad=0.06)

dlist = ["P300_hf001", "P300_hfP", "P200_hfP", "P100_hfP", "P300_hf1", "P300_hfinf"]
dlistfb = []
for dname in dlist: dlistfb.append(dname + '_fb')
# plotting also fb50 results (with dashed lines)
dlist += dlistfb 
clist = ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]
tlist = ["300", "300", "200", "100", "300", "300"]
hflist = ["$=0.01$", "$=0.1212$", "$=0.1212$", "$=0.1212$", "$=1.0$", r"$\rightarrow\infty$"]
#ax2.scatter(lonmin[0]*np.pi/180., colmin[0], c='red', s=100, marker='s')
for dname in dlist:
    omgfile = "../../" + str(dname) + "/run/omega.dat"
    if not(os.path.exists(omgfile)): continue
    idl = dlist.index(dname) % len(dlistfb)
    colL = np.loadtxt(omgfile, usecols=(5))
    lonL = np.loadtxt(omgfile, usecols=(6))
    timew = np.loadtxt(omgfile, usecols=(0))
    MIAfile = "../../" + dname + "/run/MIA.dat"
    colLf = np.loadtxt(MIAfile, usecols=(4))
    lonLf = np.loadtxt(MIAfile, usecols=(5))
    lstl = ':'
    lbl = ''
    if ('fb' in dname):
        lstl = '-'
        lbl = tlist[idl] + r" m, $M_\mathrm{h}/M$" + hflist[idl]
    scolor = clist[idl] # next(ax._get_lines.color_cycle) #['color']
    im2 = ax2.plot(lonL*np.pi/180., colL, linestyle=lstl, clip_on=False, c=scolor, label=lbl)
    print ("Eq. pos: %.3f %.3f" % (lonLf[2], colLf[2]))
    if(lstl=='-'): im2 = ax2.plot(lonLf[2]*np.pi/180., colLf[2], 'o', clip_on=False, c=scolor, markersize=10)
    
ax2.yaxis.set_major_formatter(ticker.FuncFormatter(degsymb))
ax2.set_rticks( np.linspace(0,  90, 4, endpoint=True) )
ax2.set_title('a) Positive load')
ax2.set_rmin(0.0)
ax2.set_rmax(1.1*max(colL))
ax2.set_thetamin(0)
ax2.set_thetamax(90)
#ax2.legend(fontsize = 15, loc='lower right', bbox_to_anchor=(1.08, 0.81), framealpha=1.0)
ax2.legend(fontsize = 15, loc='lower right', bbox_to_anchor=(0.65, 0.25), framealpha=1.0)
ax2.grid(True)
ax2.plot(np.linspace(0,-90,100)*np.pi/180, np.linspace(90,90,100), '--', linewidth=1, c='gray')
ax2.scatter(0, 0, marker='*', s=100, clip_on=False, c='black', label='rot. axis')
ax2.scatter(0, 90,marker='+', s=100, clip_on=False, c='black', label='tidal axis')
#plt.colorbar(im2,ax=ax2,format="%.0f", orientation="horizontal", label='time [My]',pad=0.5)

# POINTS (negative load means that the largest moment, corresponding to the surface depression, is positive)
neg_load = (abs(mommax[0]) > abs(mommin[0]))
if neg_load:
    if fps: lbl = ("load, %.1f°, %.1f°" % (colmax[0],lonmax[0]) )
    else: lbl = ("point-load position")
    ax.scatter(lonmax[0]*np.pi/180., colmax[0], s=100, marker='o', facecolors='none', edgecolors='red', label='n_'+lbl)
else:
    if fps: lbl = ("load, %.1f°, %.1f°" % (colmin[0],lonmin[0]) )
    else: lbl = ("point-load ($[1.5,3,6]{\cdot}10^{17}$ kg)")
    ax.scatter(lonmin[0]*np.pi/180., colmin[0], c='red', s=100, marker='s', label=lbl)
if fps: lbl = ("Eq. MIA, t=0: %.1f°, %.1f°" % (colmax[1], lonmax[1]))
else: lbl = ("Eq. MIA")
#ax.scatter(lonmax[1]*np.pi/180., colmax[1], edgecolors='blue', facecolors='none', s=200, marker='o', label=lbl)
if fosslit:
    if fps: lbl = ("fluid MIA: %.1f°, %.1f°" % (colmax[3], lonmax[3]))
    else: lbl = (r"load+f.b., MIA, $1.5{\cdot}10^{17}$ kg, $M_\mathrm{h}/M\rightarrow\infty$")
    ax.scatter(lonmax[3]*np.pi/180., colmax[3], edgecolors='brown', facecolors='none', s=200, marker='*', label=lbl)
if taxplot:
    if fps: lbl = ("Eq. mIA, t=0: %.1f°, %.1f°" % (colmin[1], lonmin[1]))
    else: lbl = ("Eq. mIA (t=0)")
    #ax3.scatter(lonmin[1]*np.pi/180., colmin[1], edgecolors='green', facecolors='none', s=200, marker='o', label=lbl)
    ax3.legend()
    ax3.set_rmax(1.1*max(coltax))
    ax3.plot(0., 90., 'o', c='black', clip_on=False)
    ax3.grid(True)
    if fosslit:
        if fps: lbl = ("load+F.B.: %.1f°, %.1f°" % (colmin[3], lonmin[3]))
        else: lbl = (r"load+f.b., mIA, $1.5{\cdot}10^{17}$ kg, $M_\mathrm{h}/M\rightarrow\infty$")
        ax3.scatter(lonmin[3]*np.pi/180., colmin[3], edgecolors='brown', facecolors='none', s=200, marker='P', label=lbl)

#ax.set_thetamin(-90)
#ax.set_thetamax(180)
ax.set_rlabel_position(180)
ax.plot(np.linspace(-90,180,200)*np.pi/180, np.linspace(90,90,200), '--', linewidth=1, c='gray')
ax.set_rmin(0.0)
ax.set_rmax(1.1*max(colw))
if twopanes and taxplot:
    ax.set_rmax(1.1*max(max(colw),max(coltax)))
ax.grid(True)
#ax.set_rticks([5, 30, 60, 90])

plt.subplots_adjust(wspace=0.05)
ax.set_title("c) body-fixed view", pad=15.0)
ax.legend(fontsize=15, loc='lower left', bbox_to_anchor=(-0.2, -0.07), frameon=True, handlelength=2, framealpha=1.0)
plt.savefig("../figs/Fig3a.png", bbox_inches='tight')
plt.clf()

plot_ball = False
if plot_ball:
    ball = (mommax[i0:] + mommed[i0:] + mommin[i0:]) / 3.
    plt.xscale('log')
    plt.xlabel('time [yr]')
    plt.ylabel('nonspherical I_{ii}, normalized')
    plt.plot(time[i0:],mommax[i0:] - ball, linewidth=2.0, label='Max')
    plt.plot(time[i0:],mommed[i0:] - ball, linewidth=2.0, label='med')
    plt.plot(time[i0:],mommin[i0:] - ball, linewidth=2.0, label='min')
    plt.savefig("../figs/MIA.png", bbox_inches='tight')
