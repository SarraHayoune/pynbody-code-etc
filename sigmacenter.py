''' calculating the velocity dispersion of the central region '''

import pynbody
import matplotlib.pyplot as plt
import glob
import numpy as np
import os

file = glob.glob('h148.cosmo50PLK.3072g3HbwK1BH.00????/h148.cosmo50PLK.3072g3HbwK1BH.00????.[0-9]*.std')
#file = glob.glob('h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096.[0-9]*.std')

# get a list of directories
import subprocess
p = subprocess.Popen(["ls h148.cosmo50PLK.3072g3HbwK1BH.00????/h148.cosmo50PLK.3072g3HbwK1BH.00????.iord |cut -d '/' -f1"],stdout=subprocess.PIPE,shell=True)
(dirs,err)=p.communicate()
#print dirs


tinyhalolimit = 100  # star particles 
currentdir =  ''

for i in range(len(file)):
    # first, check if we're in a new timestep. if so, 
    # close prior file and make a new one
    thisdir = file[i].split("/")
    if thisdir[0] != currentdir:
        if i != 0:
            f.close()
        currentdir = thisdir[0]
        print "new dir ",thisdir[0]
        f = open(thisdir[0]+"/sigmas.dat","w")

 #   if (i == 0):
 #       continue  # skip halo 1
    s = pynbody.load(file[i])
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
    #pynbody.plot.sph.image(s.g,qty='rho',units='g cm^-2',width=100,cmap='Greys')
    #pynbody.plot.stars.render(s,width='10 kpc')
    #radius = pynbody.derived.r(s.stars)
    #print min(radius),max(radius)
    radius = '1.5 kpc'  
    rfilt = pynbody.filt.Sphere(radius,(0,0,0))
    #print max(pynbody.derived.r(s.stars[rfilt]))
    bulge = s.stars[rfilt]
    sigmax = np.std(bulge['vel'][0])
    sigmay = np.std(bulge['vel'][1])
    sigmaz = np.std(bulge['vel'][2])
    sigma = np.sqrt(sigmax*sigmax+sigmay*sigmay+sigmaz*sigmaz)
    sigma1D = sigma/np.sqrt(3.)
    print sigma1D," km/s for halo ",haloid
    # pynbody can maybe do this?
    #v_disp = pynbody.derived.v_disp(s.stars[rfilt])
    #print v_disp
    # this does not work.  oh well moving on
    
    ### ok, find the BHs  ###
    BHfilter = pynbody.filt.LowPass('tform',0.0)
    BH = s.stars[BHfilter]
    print len(BH)," BHs in this halo"
    for j in range(len(BH)):
        print "BH mass: ", BH[j]['mass']
        BHmassstring = str(BH[j]['mass'])
        BHmassstring = BHmassstring[1:-1]  # strip off the brackets
        ## is the BH in the center?  ##
        if pynbody.derived.r(BH[j]) < 0.2:  # is the Bh within 1 kpc
            central = "yes"
        else:
            central = "no"
        data = str(sigma1D)+"   "+BHmassstring+"   "+haloid+"  "+central+"\n"
        f.write(data)
    
f.close()
