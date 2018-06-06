''' testing to see if dwarfs are reasonable '''
import pynbody
import matplotlib.pyplot as plt

file ='h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'
statfile = file+".amiga.stat"

import readcol
stats = readcol.readcol(statfile,twod=False,skipline=1)
haloid =  stats[0]

s = pynbody.load(file)
s.physical_units()
h = s.halos()
mags = []

"""
M,sigma,N = pynbody.analysis.hmf.halo_mass_function(s,log_M_min=5.0,log_M_max=13.0)

plt.plot(M,N)
plt.xlabel("Log Mass")
plt.ylabel("N of halos in bin")
plt.yscale("log")
plt.xscale("log")
#plt.xlim(5,13)
plt.show()
"""
i=0
limit = 1500.0

while len(h[haloid[i]].dark) > limit:
    print "halo ",haloid[i]
    if len(h[haloid[i]].s) == 0:
        i=i+1
        print "no stars"
        continue
    mag = pynbody.analysis.luminosity.halo_mag(h[haloid[i]].s)
    print mag
    mags.append(mag)
    i=i+1

print mags

### plotting ###
plt.hist(mags)
plt.xlabel("V-band Magnitude")
plt.show()
