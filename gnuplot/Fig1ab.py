import numpy as np
import matplotlib.pyplot as plt
from pylab import figure, show, legend, xlabel
import os.path
import sys

stuff = sys.argv[1].split(',')
for item in stuff:
    exec(item)

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = [20, 9]
#plt.rcParams.update({'lines.linewidth': 5})
fig1 = figure()
ax = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)
dlist = ["P15", "P15_Hu", "P30", "P30_Hu", "P60", "P60_Hu"]
tlist = ["1.5", "1.5", "3.0", "3.0", "6.0", "6.0"]
clist = ["C0", "C3", "C1", "C4", "C2", "C5", "C6", "C7", "C8"]

for dname in dlist:
    My = 1.e6
    tpwfile = "../../" + dname + "/run/tpw_process.dat"
    time = np.loadtxt(tpwfile, usecols=(0)) / My
    omg = np.loadtxt(tpwfile, usecols=(1,2,3))
    tax = np.loadtxt(tpwfile, usecols=(6,7,8))

    rtd = 180./np.pi
    ampo = np.sqrt(np.power(omg[:,0],2) + np.power(omg[:,1],2) + np.power(omg[:,2],2))
    ampt = np.sqrt(np.power(tax[:,0],2) + np.power(tax[:,1],2) + np.power(tax[:,2],2))
    if ("Hu" in dname):
        lblc = ''; lbll = '';
        lstl = 'dotted'
        lwd = 7
    else:
        lblc = tlist[dlist.index(dname)] + r"$\cdot 10^{17}$ kg, colatitude"
        lbll = tlist[dlist.index(dname)] + r"$\cdot 10^{17}$ kg, longitude"
        lstl = '-'
        lwd = 5
    ind = dlist.index(dname)
    ax.plot(time[:], np.arccos(omg[:,2]/ampo[:])*rtd, linestyle=lstl, c=clist[ind], lw=lwd, label=lblc)
    ax.plot(time[:], np.arctan2(omg[:,1],omg[:,0])*rtd, linestyle=lstl, c=clist[ind+3], lw=lwd, label=lbll)
    ax2.plot(time[:], np.arccos(tax[:,2]/ampt[:])*rtd, linestyle=lstl, c=clist[ind], lw=lwd, label=lblc)
    ax2.plot(time[:], np.arctan2(tax[:,1],tax[:,0])*rtd, linestyle=lstl, c=clist[ind+3], lw=lwd, label=lbll)

ax.set_title("comparison of $\omega(t)$")
ax2.set_title("comparison of $o(t)$")
ax.set_ylabel('North pole [°]')
ax2.set_ylabel('Sub-host point [°]')
ax.set_xlabel('Time [My]')
ax2.set_xlabel('Time [My]')
ax.legend(fontsize=16) #loc='upper left', bbox_to_anchor=(-0.1, 1.0), framealpha=1.0)
#ax2.legend(fontsize=16)
plt.savefig("../figs/Fig1ab.png", bbox_inches='tight')
plt.clf()
