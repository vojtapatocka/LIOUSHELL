import numpy as np
import matplotlib.pyplot as plt
from pylab import figure, show, legend, xlabel
import os.path
import sys
from mpl_toolkits.mplot3d import Axes3D
import numpy.linalg as linalg

def generate_specific_rows(filePath, userows=[]):
    with open(filePath) as f:
        for i, line in enumerate(f):
            if i in userows:
                yield line

unit = 1.e29                
print ("LOADING DATA: You better check Jtensors.dat to see if the order of J constituents matches the one below")                
#Jfile = '../run/Jtensors.dat'
Jfile = '../Jtensors_N100fb50.dat'
gen = generate_specific_rows(Jfile, userows=[1,2,3])
hydb = np.loadtxt(gen, unpack='true') / unit
gen = generate_specific_rows(Jfile, userows=[5,6,7])
load = np.loadtxt(gen, unpack='true') / unit
gen = generate_specific_rows(Jfile, userows=[9,10,11])
eqMIA = np.loadtxt(gen, unpack='true') / unit
gen = generate_specific_rows(Jfile, userows=[13,14,15])
fload = np.loadtxt(gen, unpack='true') / unit
gen = generate_specific_rows(Jfile, userows=[17,18,19])
fossb = np.loadtxt(gen, unpack='true') / unit
#gen = generate_specific_rows(Jfile, userows=[21,22,23])
#final = np.loadtxt(gen, unpack='true') / unit
final = fload + fossb

# Removing the spherical parts of contributions with a non-zero mass
hydb -= np.trace(hydb)/3.*np.identity(3)

print("hydb", hydb)
print("fossb", fossb)

# shows an ellipsoid that represents the hydrostatic bulges (0 - ring, 1 - 3D wire, 2 - 2D wire)
def plot_ellipsoid(A,lbl,barva,typell):
    UnitA, s, rotation = linalg.svd(A)
    radii = s
    print("SVD radii, radii[2]/radii[1] ", radii, radii[2]/radii[1])
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    if typell>1:
        maxk = (s == np.max(s))
        radii[maxk] = 0.0
    if typell>0:
        v = np.linspace(0.0, np.pi, 100)
        x = radii[0] * np.outer(np.cos(u), np.sin(v))
        y = radii[1] * np.outer(np.sin(u), np.sin(v))
        z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
        for i in range(len(x)):
            for j in range(len(x)):
                [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotation)
        ax.plot_wireframe(x, y, z, rstride=10, cstride=10, color=barva, alpha=0.05)
    else:
        x = 0 * np.cos(u)
        y = radii[1] * np.cos(u)
        z = radii[2] * np.sin(u)
        for i in range(len(x)):
            [x[i],y[i],z[i]] = np.dot([x[i],y[i],z[i]], rotation)
        ax.plot(x, y, z,  color=barva, alpha=1.0, label=lbl)
        
def plot_eig(A,lbl,barva,thck):
    enum, evec = linalg.eig(A)
    maxk = (np.abs(enum) == np.max(np.abs(enum)))
    mink = np.logical_and(np.abs(enum) < 1.001*np.min(np.abs(enum)), np.abs(enum) > 0.999*np.min(np.abs(enum)))
    print(lbl," id_max, id_min, values: ", maxk, mink, enum)
    xco = enum[maxk]*evec[0,maxk]
    yco = enum[maxk]*evec[1,maxk]
    zco = enum[maxk]*evec[2,maxk]
    orig = np.zeros(len(xco))
    # ELONGATION: the largest eigenvalue is negative (elongation along the symmetry axis)
    if enum[maxk]<0.0:
        ax.quiver(orig, orig, orig, xco, yco, zco, color=barva, label=lbl)
        ax.quiver(orig, orig, orig, -xco, -yco, -zco, color=barva)
        # The inertia ellipsoid is not axially symmetrical
        if np.count_nonzero(mink)==1: 
            #medk = np.logical_and(maxk==False, mink==False)
            medk = (enum == np.max(enum))
            xco = enum[medk]*evec[0,medk]
            yco = enum[medk]*evec[1,medk]
            zco = enum[medk]*evec[2,medk]
            ax.quiver(xco, yco, zco, -xco, -yco, -zco, color=barva)
            ax.quiver(-xco, -yco, -zco, xco, yco, zco, color=barva)
    # FLATTENING: the largest eigenvalue is positive (flattening along the symmetry axis)
    else: plot_ellipsoid(A,lbl,barva,0)

fig = plt.figure(figsize=(8, 12))
#plt.rcParams["figure.figsize"] = [20, 16]
plt.rcParams.update({'font.size': 14})
ax = fig.add_subplot(111, projection='3d')

# MAIN CODE
plot_ellipsoid(hydb,'Hydrostatic bulge','blue',0)
plot_eig(fload,'Disc load','orange',2)
plot_eig(final,'Load + Fossil b.','green',4)

ax.scatter(0,0,hydb[2,2],marker='+',c='r',s=60,label='North-pole')
ax.scatter(hydb[0,0],0,0,marker='x',c='r',s=60,label='Sub-host point')
enum, evec = np.linalg.eig(final)
maxk = (enum == np.max(enum))
mink = (enum == np.min(enum))
medk = np.logical_and(maxk==False,mink==False)
fmaxi = np.zeros(3)
fmedi = np.zeros(3)
fmini = np.zeros(3)
omgini = np.array([0.,0.,1.])
taxini = np.array([1.,0.,0.])
# Somehow this works although maxk,mink,medk are fields and not identifiers
for i in range(0,3):
    fmaxi[i] = evec[i,maxk]
    fmini[i] = evec[i,mink]
    fmedi[i] = evec[i,medk]
# Flipping eigenvectors - choosing the one closer to the initial configuration
if (np.arccos(np.dot(fmaxi,omgini)) > 0.5*np.pi):
    fmaxi *= -1.
if (np.arccos(np.dot(fmini,taxini)) > 0.5*np.pi):
    fmini *= -1.
print("fmaxi", fmaxi)
print("fmedi", fmedi)
print("fmini", fmini)
print("final: maxk, mink, enum: ", maxk, mink, enum)

def Rotmat(axis, theta):
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def porange(avec, bvec, amp, npoints, fullsp):
    linco = np.zeros((3,npoints+1))
    for i in range(0,npoints+1):
        fc = i/npoints
        fcn = 1.0 - fc
        if fullsp==1:
            rotang = 2*np.pi * fc
            maxv = np.dot(Rotmat(bvec,rotang), avec) * amp
        else:
            maxc = avec[:]*fc + bvec[:]*fcn
            maxv = maxc[:]/np.power(np.dot(maxc,maxc),0.5) * amp
        linco[:,i] = maxv[:]
    return linco

# Plotting polar "TPW"
if np.count_nonzero(maxk)==2:
    omgspan = porange(fmaxi, fmini, hydb[2,2], 20, 1)
    ax.scatter(omgspan[0,:],omgspan[1,:],omgspan[2,:],marker='+',c='green',s=60)
else:
    omgdirp = porange(fmaxi, omgini, hydb[2,2], 40, 0)
    ax.plot(omgdirp[0,:],omgdirp[1,:],omgdirp[2,:],'-',c='black')
    ax.scatter(omgdirp[0,-1],omgdirp[1,-1],omgdirp[2,-1],marker='+',c='green',s=60)
    
# Plotting sub-host point "TPW"
if np.count_nonzero(mink)==2:
    taxspan = porange(fmini, fmaxi, hydb[0,0], 20, 1)
    ax.scatter(taxspan[0,:],taxspan[1,:],taxspan[2,:],marker='x',c='green',s=60)
else:
    taxdirp = porange(fmini, taxini, hydb[0,0], 40, 0)
    ax.plot(taxdirp[0,:],taxdirp[1,:],taxdirp[2,:],'--',c='black')
    ax.scatter(taxdirp[0,-1],taxdirp[1,-1],taxdirp[2,-1],marker='x',c='green',s=60)

ax.set_xlabel("\n$x\, [10^{29}$ kg m$^2]$", linespacing=2)
ax.set_ylabel("\n$y\, [10^{29}$ kg m$^2]$", linespacing=2)
ax.set_zlabel("\n$z$") # [10^{29}$ kg m$^2]$")
plt.legend(loc="upper left")
ax.text(0.7,0.8,4.5,'f) N 100 m, EL 50 km',fontsize=18)
ax.grid(False)
#ax.set_aspect('equal')
unirnge = np.max(hydb)
#plot_ellipsoid(unirnge*np.identity(3),1,'black')
ax.set_xlim(-unirnge/2,unirnge/2)
ax.set_ylim(-unirnge/2,unirnge/2)
ax.set_zlim(-unirnge/2,unirnge)
ax.view_init(elev=15., azim=230.)
plt.savefig('../figs/Ellipsoid.png',bbox_inches='tight')
#plt.show()

