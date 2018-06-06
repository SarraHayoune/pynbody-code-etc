import pynbody
import glob

filelist = glob.glob('h148.cosmo50PLK.3072g3HbwK1BH.00????/h148.cosmo50PLK.3072g3HbwK1BH.00????.[0-9]*.std')


#file = 'h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096.3.std'

for file in filelist:
    s=pynbody.load(file)
    if len(s.star) == 0 or len(s.gas) == 0 or len(s.dark) == 0:
        continue
    s.physical_units()
    pynbody.analysis.decomp(s,aligned=True,angmom_size='1 kpc')
    pynbody.array.SimArray.write(s.star['decomp'],overwrite=True)
