import matplotlib.pyplot as plt
import matplotlib
import readcol
import subprocess
import numpy as np

### function to make 2nd legend, from stack exchange ###
def create_proxy(label):
    line = matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='black',
                mec='none', marker=r'$\mathregular{{{}}}$'.format(label))
    return line

##########################################
#### read in the McConnell and Ma data ###
##########################################
mmsigma,mmsigmaerr,mmbh,mmexp,mmerr1,mmerr2 = readcol.readcol('mcconnellma.dat',twod=False)
mmbh = mmbh*10.0**mmexp
mmerr1 = mmerr1*10**mmexp
mmerr2 = mmerr2*10**mmexp
alpha = 8.39  # 8.39: late type  8.32: total
beta = 5.05  # 5.05: late type  5.64:  total
alphaerr = 0.06   # 0.06:  late type  0.05:  total
betaerr = 1.16    # 1.16: late type  0.32:  total

### McConnell and Ma line ###
sigmaline = np.arange(500)+1.
massline = alpha + beta*np.log10(sigmaline/200.)
massline = 10**massline


files = ['h148.cosmo50PLK.3072g3HbwK1BH.004096/sigmas.dat',
'../h229/h229.cosmo50PLK.3072gst5HbwK1BH.004096/sigmas.dat',
         '../h242/h242.cosmo50PLK.3072gst5HbwK1BH.004096/sigmas.dat',
         '../cptmarvel/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/sigmas.dat',
         '/data/elektra/elektra.cosmo25cmb.4096g5HbwK1BH.004096/sigmas.dat','/data/rogue/rogue.cosmo25cmb.4096g5HbwK1BH.004096/sigmas.dat','../storm/storm.cosmo25cmb.4096g5HbwK1BH.004096/sigmas.dat']
#files = ['h148.cosmo50PLK.3072g3HbwK1BH.004096/sigmas.dat']

satsymbol ='o'
fieldsymbol = '*'
labels = [satsymbol,fieldsymbol]
proxies = [create_proxy(item) for item in labels]
words = ["satellites","field"]

for file in files:
    print file
    sigma,mass,haloid,central = readcol.readcol(file,twod=False)
    sigma = np.array(sigma)
    print max(sigma)
    mass = np.array(mass)
    yay = np.where(sigma < 100.)
    sigma = sigma[yay]
    mass = mass[yay]
    central = central[yay]
    isitcentral=[]
    for i in range(len(central)):
        if central[i] == 'yes':
            isitcentral.append(i)

    if file == files[3] or file == files[4] or file == files[5]:
        symbol = fieldsymbol
    else:
        symbol = satsymbol



    plt.plot(sigma,mass,'b'+symbol,label='non-central BHs' if file == files[0] else '')
    plt.plot(sigma[isitcentral],mass[isitcentral],"r"+symbol,label='central BHs' if file == files[0] else '')

plt.plot(mmsigma,mmbh,"go",label='McConnell & Ma')
plt.plot(sigmaline,massline,'b--')
plt.yscale("log")
plt.xscale("log")
plt.xlim(4,500)
plt.ylim(1e4,1e11)
plt.xlabel("$\sigma$ (km/s)",fontsize=20)
plt.ylabel("BH Mass ($M_{\odot}$)",fontsize=20)
#plt.title("all sims z=0 dwarfs")
plt.legend(loc='upper left',numpoints = 1)
#plt.legend(proxies, words, numpoints=1,loc='lower right')
#plt.legend([satsymbol,fieldsymbol],["satellites","field"],loc="lower right",scatterpoints=1)
plt.show()
#plt.savefig("msigma.png")
