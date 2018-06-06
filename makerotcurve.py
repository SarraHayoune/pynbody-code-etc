import pynbody
import numpy as np

files = np.loadtxt("files.list",unpack=True,dtype='str')
nfiles = len(files)
#file = 'h148.cosmo50PLK.3072g3HbwK1BH.002048/h148.cosmo50PLK.3072g3HbwK1BH.002048'

for file in files:
    s=pynbody.load(file)
    h=s.halos()
    BHindex = pynbody.filt.LowPass('tform',0.0)
    BHs = s.star[BHindex]
    halos = BHs['amiga.grp']
    print halos

    s.physical_units()
    #import matplotlib.pyplot as plt
    halos = np.unique(halos)
    print halos

    for i in range(len(halos)):
        pynbody.analysis.angmom.faceon(h[halos[i]])
        p = pynbody.analysis.profile.Profile(h[halos[i]].s)
        rotfile = file+'.'+str(halos[i])+".rot"
        print "making "+rotfile
        f=open(rotfile,"w")
        print len(p['rbins'])," length of profile"
        for j in range(len(p['rbins'])):
            thingtowrite = p['rbins'][j],p['v_circ'][j]
            #print thingtowrite
            thingtowrite = str(thingtowrite)[1:-1]+'\n'
            f.write(thingtowrite)
        f.close()
        
        #plt.plot(p['rbins'],p['v_circ'])
        #plt.xlabel('radius')
        #plt.ylabel('v_circ')
        #plt.xscale('log')
        #plt.title('halo '+str(halos[i]))
        #plt.show()

