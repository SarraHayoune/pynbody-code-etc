''' calculating the velocity dispersion of the central region '''
# convert this to specifically calculating the DF-related sigma
# only need this for z=0
import pynbody
import matplotlib.pyplot as plt
import glob
import numpy as np
import os

filename = 'h148.cosmo50PLK.3072g3HbwK1BH.004096'

file = glob.glob(filename+'/'+filename+'.[0-9]*.std')
#file = glob.glob('h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096.[0-9]*.std')

# get a list of directories
#import subprocess
#p = subprocess.Popen(["ls h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096.iord |cut -d '/' -f1"],stdout=subprocess.PIPE,shell=True)
#(dirs,err)=p.communicate()
dirs = filename

def findBH(snap):
    BHfilter = pynbody.filt.LowPass('tform',0.0)
    BH = snap.stars[BHfilter]
    return BH

tinyhalolimit = 100  # star particles 
currentdir =  ''

f = open(dirs+"/dfsigmas.dat","w")

for i in range(len(file)):
    s = pynbody.load(file[i])
    bh = findBH(s)
    nbh = len(bh)
    print "number of BHs = ",nbh
    if len(s.stars) < tinyhalolimit:  # tiny halos don't get kinematics
        print "skipping tiny halo "
        continue
    s.physical_units()
    #print " numbers of particles ", len(s.g),len(s.d),len(s.s)
    # get the halo id from the file name
    splitfile = file[i].split(".")
    haloid = splitfile[7]
    #print haloid," haloid"

    pynbody.analysis.halo.center(s,mode='ssc')
    r = pynbody.derived.r(bh)
    # loop through each BH 
    j=0  # need to keep track of BH index too
    for radius in r:
        print radius
        if radius < 0.1:
            j +=1
            continue
        #    radius = '1.5 kpc'  
        rfilt = pynbody.filt.Sphere(radius,(0,0,0))
        #print max(pynbody.derived.r(s.stars[rfilt]))
        bulge = s[rfilt]
        sigmax = np.std(bulge['vel'][0])
        sigmay = np.std(bulge['vel'][1])
        sigmaz = np.std(bulge['vel'][2])
        sigma = np.sqrt(sigmax*sigmax+sigmay*sigmay+sigmaz*sigmaz)
        sigma1D = sigma/np.sqrt(3.)
        print sigma1D," km/s for halo ",haloid," radius ",radius
    
        print "BH mass: ", bh[j]['mass']
        BHmassstring = str(bh[j]['mass'])
        BHmassstring = BHmassstring[1:-1]  # strip off the brackets
        print BHmassstring
        strradius = str(radius)
#        BHiord = bh['iord'][j]
#        BHiord = str(BHiord[1:-1])

        data = str(sigma1D)+"    "+haloid+"  "+strradius+"    "+BHmassstring+"\n"
        f.write(data)
        j+=1
        
f.close()
