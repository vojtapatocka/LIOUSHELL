import numpy as np
import matplotlib.pyplot as plt
from pylab import figure, show, legend, xlabel
import os.path
import sys

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

i0 = 3
fosslit = False
if (Mlabel[3]=="FLOAD_+_F.B."):
    i0=4
    fosslit = True
twopanes = False
taxplot = False
iterations = False
MIAplot = True

tidfile = "../run/tidax.dat"
if os.path.exists(tidfile):
    coltax = np.loadtxt(tidfile, usecols=(2))
    lontax = np.loadtxt(tidfile, usecols=(3))
    taxplot = True

stuff = sys.argv[1].split(',')
for item in stuff:
    exec(item)

plt.rcParams.update({'font.size': 16})
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
    im4 = ax3.scatter(lonmin[i0:]*np.pi/180., colmin[i0:], s=70, c=cmia, marker='+', cmap='cool', label='mIA (cool palette)')
    if fosslit:
        lbl = ("Load+F.B.: %.1f°, %.1f°" % (colmin[3], lonmin[3]))
        ax3.scatter(lonmin[3]*np.pi/180., colmin[3], edgecolors='brown', facecolors='none', s=200, marker='P', label=lbl)
    lbl = ("tidax, f.p.: %.1f°, %.1f°" % (coltax[-1], lontax[-1]))
    im3 = ax3.scatter(lontax*np.pi/180., coltax, marker='s', s=10, c=comg, cmap='viridis', label=lbl)
    
lbl = ("Load, f.p.: %.1f°, %.1f°" % (colL[-1],lonL[-1])) 
im2 = ax2.scatter(lonL*np.pi/180., colL, c=timew/ky, cmap='viridis', label=lbl)
ax2.set_title('bulge-frame view')
ax2.scatter(0, 0, s=100, marker='x', c='red', clip_on=False, label='rot. axis')
ax2.scatter(0, 90, s=100, marker='+', c='red', clip_on=False, label='tidal axis')
lbl = ("Load+F.B.: %.1f°, %.1f°" % (colmax[2], lonmax[2]))
ax2.scatter(lonmax[2]*np.pi/180., colmax[2], edgecolors='blue', facecolors='none', s=200, marker='o', label=lbl)
ax2.set_rmin(0.0)
ax2.set_rmax(1.1*max(colL))
thmin = 0.; thmax = 360.
#thmin = -90.; thmax = 0.
ax2.set_thetamin(thmin)
ax2.set_thetamax(thmax)
ax2.grid(True)
ax2.plot(np.linspace(thmin,thmax,200)*np.pi/180, np.linspace(90,90,200), '--', linewidth=1, c='gray')
ax2.legend()

lbl = ("$\omega$, f.p.: %.1f°, %.1f°" % (colw[-1],lonw[-1]))
im1 = ax.scatter(lonw*np.pi/180., colw, marker='.', c=comg, cmap='viridis', label=lbl)
plt.colorbar(im1, ax=ax3, format="%.0f", orientation="horizontal", label='time [My]', pad=0.01)

# POINTS (negative load means that the largest moment, corresponding to the surface depression, is positive)
neg_load = (abs(mommax[0]) > abs(mommin[0]))
if neg_load:
    lbl = ("load, %.1f°, %.1f°" % (colmax[0],lonmax[0]) )
    ax.scatter(lonmax[0]*np.pi/180., colmax[0], s=100, marker='o', facecolors='none', edgecolors='red', label='n_'+lbl)
else:
    lbl = ("load, %.1f°, %.1f°" % (colmin[0],lonmin[0]) )
    ax.scatter(lonmin[0]*np.pi/180., colmin[0], c='red', s=100, marker='s', label='p_'+lbl)
lbl = ("eq MIA ($t_\mathrm{gr}$): %.1f°, %.1f°" % (colmax[1], lonmax[1]))
ax.scatter(lonmax[1]*np.pi/180., colmax[1], edgecolors='blue', facecolors='none', s=200, marker='o', label=lbl)
if taxplot:
    lbl = ("eq mIA ($t_\mathrm{gr}$): %.1f°, %.1f°" % (colmin[1], lonmin[1]))
    ax3.scatter(lonmin[1]*np.pi/180., colmin[1], edgecolors='green', facecolors='none', s=200, marker='o', label=lbl)
    ax3.legend()
    ax3.plot(np.linspace(0,360,400)*np.pi/180, np.linspace(90,90,400), '--', linewidth=1, c='gray')
    ax3.set_rmax(1.1*max(coltax))
    ax3.grid(True)
if fosslit:
    lbl = ("Load+F.B.: %.1f°, %.1f°" % (colmax[3], lonmax[3]))
    ax.scatter(lonmax[3]*np.pi/180., colmax[3], edgecolors='brown', facecolors='none', s=200, marker='*', label=lbl)

ax.plot(np.linspace(0,360,400)*np.pi/180, np.linspace(90,90,400), '--', linewidth=1, c='gray')
ax.set_rmin(0.0)
ax.set_rmax(1.1*max(colw))
if twopanes and taxplot:
    ax.set_rmax(1.1*max(max(colw),max(coltax)))
#ax.set_rticks([5, 30, 60, 90])
#ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax.grid(True)

ax.set_title("body-fixed view")
ax.legend() #loc='upper left', bbox_to_anchor=(-0.1, 1.0), framealpha=1.0)
plt.savefig("../figs/polar.png", bbox_inches='tight')
plt.clf()

plt.rcParams["font.size"] = 18
plt.rcParams["figure.figsize"] = [12, 9]
fig1 = figure()
ball = (mommax[i0:] + mommed[i0:] + mommin[i0:]) / 3.
nrm = 1.e4
plt.xscale('log')
plt.xlabel('Time [yr]')
plt.xlim((1.1e5,4.e9))
plt.ylabel(r"Moments of Inertia [$\times 10^4$]")
plt.plot(time[i0:],(mommax[i0:] - ball) / ball*nrm, c='C0', linewidth=3.0, label='maximum')
plt.plot(time[i0:],(mommed[i0:] - ball) / ball*nrm, c='C1', linewidth=3.0, label='medium')
plt.plot(time[i0:],(mommin[i0:] - ball) / ball*nrm, c='C2', linewidth=3.0, label='minimum')
if len(compdir)>1:
    MIAfile = "../../" + compdir + "/run/MIA.dat"
    mommax2 = np.loadtxt(MIAfile, usecols=(1))
    time2 = np.loadtxt(MIAfile, usecols=(0))
    mommin2 = np.loadtxt(MIAfile, usecols=(3))
    mommed2 = np.loadtxt(MIAfile, usecols=(2))
    ball2 = (mommax2[i0:] + mommed2[i0:] + mommin2[i0:]) / 3.
    plt.plot(time2[i0:],(mommax2[i0:] - ball2) / ball2*nrm, '--', c='C0', linewidth=3.0)
    plt.plot(time2[i0:],(mommed2[i0:] - ball2) / ball2*nrm, '--', c='C1', linewidth=3.0)
    plt.plot(time2[i0:],(mommin2[i0:] - ball2) / ball2*nrm, '--', c='C2', linewidth=3.0)
plt.legend() #loc="center left")
plt.savefig("../figs/MIA.png", bbox_inches='tight')
