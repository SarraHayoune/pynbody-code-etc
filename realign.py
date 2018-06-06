import pynbody
import numpy
lines = [line.strip() for line in open('file.list')]

for i in lines:
    infile = i
    outfile = i + '.updated'
    s=pynbody.load(infile)
    h=s.halos()
    s.physical_units()
    pynbody.analysis.angmom.faceon(h[1])
    s.write(fmt=pynbody.tipsy.TipsySnap, filename=outfile)
