# plot occupation fraction for a combined set of simulations
import pynbody
import matplotlib.pyplot as plt
import numpy as np
import readcol


### inputs ###
files = ['h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096','../h229/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096','../h242/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096','../h329/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096']
#files = ['h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096']

centrals = True
rlim = 1
cum=True

# setting up the mass bins
logmassbins = np.arange(15)/2.+6.5
MVirMin = logmassbins.min()
MVirMax = logmassbins.max()
bins = len(logmassbins)
dLogM = 0.5
ngalaxies = 1000
occfrac = []
MassBinFlag = np.zeros(bins)
allMVirs = []
BHFlag = []
kk=0

for file in files:
    stat = readcol.readcol(file+'.amiga.stat',skipline=1)
    Mvir = stat[:,5].astype('float')
    allMVirs.append(Mvir)
    print "MVir length = ",len(allMVirs)
    s = pynbody.load(file)
    h = s.halos()
    bhgrps = s.stars['amiga.grp'][(s.stars['tform']<0)]
#    print "calculating Min and Max halo masses..."
#    MvirMin = h[int(bhgrps.max())]['mass'].in_units('Msol').sum()
#    MvirMax = h[1]['mass'].in_units('Msol').sum()
#    print "Min: ", MvirMin, " Max: ", MvirMax
#    dLogM = (np.log10(MvirMax) - np.log10(MvirMin))/bins
    BHFlag.append(np.zeros(bhgrps.max()))
    ugrps = np.unique(bhgrps)
    print "there are ", len(ugrps), "halos with BHs"
    print "determining existence of central BHs..."
    for i in ugrps:
        if i <2: continue  # for now skip halo 1, for quickness
        print "halo ", i
        cen = pynbody.analysis.halo.shrink_sphere_center(h[i])
        h[i].stars['pos'] -= cen
        if centrals == True:
            if len(h[i].stars[((h[i].stars['tform']<0)&(h[i].stars['r'].in_units('kpc')<rlim))])>0: BHFlag[kk[i-1]] = 1
        else:
            BHFlag[kk[i-1]] = 1
        h[i].stars['pos'] += cen
        #  which mass bin is it in?
        loghalomass = np.log10(Mvir[i])
        print loghalomass," halo mass"
        found = False
        j=0
        while (found == False or j > 1000):
            #print j
           # print  MVirMin+dLogM*np.float(j), MVirMin+dLogM*np.float(j+1)
            if (loghalomass > MVirMin+dLogM*np.float(j) and loghalomass <= MVirMin+dLogM*np.float(j+1)):
                #print j," now j"
                #print loghalomass, " mass between ",MVirMin+dLogM*np.float(j), MVirMin+dLogM*np.float(j+1)
                MassBinFlag[j] += 1
                found = True
            j += 1

    print MassBinFlag
    kk += 1  # file counter

#print len(MassBinFlag),len(allMVirs),bins

masses = [j for i in allMVirs for j in i]
BHFlag = [j for i in BHFlag for j in i]

occNum, Mbins = np.histogram(np.log10(masses),bins=bins,weights = BHFlag,range=(6.5,13))
HNum, Mbins2 = np.histogram(np.log10(masses),bins=bins,range=(6.5,13))
if cum==True:
    occFrac = np.cumsum(occNum)/np.cumsum(HNum).astype('float')
else:
    occFrac = occNum/HNum.astype('float')
    
	

#Mbins,occfrac = bhanalysis.plotOccFrac(file)
# Mbins has one extra element compared to occfrac
# take the average as the data point for plotting
Mbinsmean = (Mbins - np.roll(Mbins,1))/2.0 + Mbins
print Mbinsmean

Mbins = Mbinsmean[1:]
print len(Mbins),len(occFrac)

plt.plot(Mbins,occFrac)
plt.xlabel("Mass (Msun)")
plt.ylabel("MBH Occupation Fraction")
plt.show()

