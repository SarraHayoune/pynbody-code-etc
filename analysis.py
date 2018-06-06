import pynbody
import numpy as np
import pandas as pd
import os
pd.set_option('display.max_column', None)
pd.set_option('display.max_rows', None)
import readcol
import gc

""" 
Let's figure out how to make the bhhalo.sav file in python
What do I need to do?
amiga grp files already exist. I need to:

- find the BHs in each snapshot
- write a mark file (Michael does this?  Michael does everything)
- find the halo centers
- find BH radial distance
- for each BH at each timestep record the halo and BH data

"""

files = readcol.readcol('files.list')
files = files[:,0]
nfiles = len(files)

def findBH(snap):
    BHfilter = pynbody.filt.LowPass('tform',0.0)
    BH = snap.stars[BHfilter]
    return BH

def findBHhalos(snap):
    BH = findBH(snap)
    BHhalos = BH['amiga.grp']
    return BHhalos

def getz(snap):
    return snap.properties['z']

def gettime(snap):
    return pynbody.analysis.cosmology.age(snap)

def writeBHmark(snap,filename):
    if not os.path.exists(filename+'.bhmark'):
        f = open(filename+".bhmark","w")
        f.write(str(len(snap))+' '+str(len(snap.gas))+' '+str(len(snap.star))+'\n')
        bhind = np.where(s.stars['tform']<0)
        bhindreal = bhind+len(snap.dark)+len(snap.gas)+1
        for ii in range(len(bhindreal)):
            f.write(str(bhindreal[ii])+'\n')
        f.close()
    else: print "bhmark file exists"

    
# initialize dataframe
columns = ['mass','haloid','BHpos','BHvel','redshift','time','bhiord','halodist','Mvir','Mstar','Mgas']
bhinfo = pd.DataFrame(columns=columns)

# create the dataframe by going thru all the files
for file in files:
    # load in the file and get the halo info
    s = pynbody.load(file)
    h = s.halos()
    s.physical_units()
    writeBHmark(s,file)
    # find the BHs
    BH = findBH(s)
    BHhalos = findBHhalos(s)
    nBH = len(BHhalos)
    # want to go through each halo one by one but not have to repeat-load or repeat-center any
    sortedhaloinds = np.argsort(BHhalos)
    print sortedhaloinds
    print BHhalos[sortedhaloinds]
    halo = 0  # initialize what halo we are on
 
    for i in sortedhaloinds:
        # which halo are we on?  need to center 
        currenthalo = BHhalos[i]
        print 'current halo: ',currenthalo
        if currenthalo != halo:  # need to center on new halo
            print "new halo calcs"
            halo = currenthalo
            pynbody.analysis.halo.center(h[currenthalo],mode='ssc')
            starmass = h[currenthalo].s['mass'].sum()
            gasmass = h[currenthalo].g['mass'].sum()
            virialmass = starmass+gasmass+h[currenthalo].d['mass'].sum()
            
            data = [[BH['mass'][i],BHhalos[i],BH['pos'][i].in_units('kpc'),BH['vel'][i],getz(s),gettime(s),BH['iord'][i],BH['r'][i],virialmass,starmass,gasmass]]
            info = pd.DataFrame(data,columns=columns)
            bhinfo = bhinfo.append(info)

    del(BH)
    del(h)
    del(s)
    gc.collect()
bhinfo.to_csv('bhinfo.csv')
    
