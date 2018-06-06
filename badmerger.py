#  find the DF timescale for BHs I suspect are
# too massive and merge too quickly
import pynbody
import numpy as np

option= 1

if (option == 1):
    file = 'h148.cosmo50PLK.3072g3HbwK1BH.000640/h148.cosmo50PLK.3072g3HbwK1BH.000640'
    thisBH = 101863741
else:
    file = 'h148.cosmo50PLK.3072g3HbwK1BH.000974/h148.cosmo50PLK.3072g3HbwK1BH.000974'
    thisBH = 101863769

s = pynbody.load(file)
h = s.halos()
s.physical_units()
pynbody.analysis.halo.center(h[1],mode='ssc')
iords = s['iord']
iord = np.where(iords == thisBH)
r = pynbody.derived.r(s[iord])
print 'radius ',r
rfilt = pynbody.filt.Sphere(r,(0,0,0))
bulge = s[rfilt]
sigmax = np.std(bulge['vel'][0])
sigmay = np.std(bulge['vel'][1])
sigmaz = np.std(bulge['vel'][2])
sigma = np.sqrt(sigmax*sigmax+sigmay*sigmay+sigmaz*sigmaz)
sigma1D = sigma/np.sqrt(3.)
print 'sigma ',sigma1D
vbh = s['vel'][iord].in_units('cm s**-1')
vbh = ((vbh[0][0]*vbh[0][0])**2 + (vbh[0][1]*vbh[0][1])**2 + (vbh[0][2]*vbh[0][2])**2)**(0.5)
print 'velocity of BH ', vbh

mass = s['mass'][iord]
print 'actual mass ',mass
mass = mass/10.   #  change mass factor for DF approximations
print 'augmented mass ',mass
# radius is r
#bmax = r.in_units('cm')
bmax = s['eps'][iord].in_units('cm')
print 'bmax ',bmax
#bmin = 0.01
bmin = 6.67e-8*mass.in_units('g')/vbh**2.0  # in cm
print 'bmin ',bmin
#bmin = bmin.in_units(kpc) 

tdf = (19.0/(np.log(bmax / bmin)))*(r/5.)*(sigma1D/200.)*(1e8/mass)
print tdf[0],' Gyr  should be the dynamical time'


# what is the actual time to merger?
# readcol times.list.  time[9] - time[3] is the answer
import readcol
time = readcol.readcol('times.list')
time = time[:,0]
if (option == 1):
    actualtime = time[10] - time[3]
else:
    actualtime = time[9] - time[2]
print actualtime,' is the actual time'
