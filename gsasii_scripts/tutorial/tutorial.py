"""This is a test code for scripting powder diffraction analysis (phase fraction, lattice strain) in GSASII for batch processing.

This script follows the tutorial provided by Jackson O'Donnell at:
subversion.xray.aps.anl.gov/pyGSAS/Tutorials/PythonScript/Script.html

This code was authored by:
Nathan S Johnson
Ph.D Student
Dept. of Mechanical Engineering
Colorado School of Mines"""

#These are the minimum packages you need just to run GSASII from the cmd line
#The second line (sys.path.insert) adds the GSASII scriptable directory to your current path
import os, sys
sys.path.insert(0,os.path.expanduser("~/g2conda/GSASII/"))
import GSASIIscriptable as G2sc

#Prints refinement results
def HistStats(gpx):
    print(u"*** profile Rwp, "+os.path.split(gpx.filename)[1])
    for hist in gpx.histograms():
        print("/t{:20s}:{:.2f}".format(hist.name,hist.get_wR()))
    print("")
    gpx.save()

#To create a new project use the following command
#NOTE: The tutorial says to use newgpx = 'filename.gpx' to create a new file
#This approach WILL NOT work -- instead use filename = 'filename.gpx'
gpx = G2sc.G2Project(filename='PBSO4.gpx')

#I am not entirely sure on the nature of G2Project variables (gpx, in this case) -- whether they are a true python variable, a class, etc

#Define a data directory and a work directory
#Data directory stores the data (profound)
datadir = "~/git/gsasii_scripts/tutorial/datadir/"
workdir = "~/git/gsasii_scripts/tutorial/datadir/"

#Let's add a histogram
#The author of the tutorial uses os.path.join to look for files in the data and work directories -- not what I would have used, but keeping it for now for the sake of following directions
 
hist1=gpx.add_powder_histogram(os.path.join(datadir,"PBSO4.XRA"),os.path.join(datadir,"INST_XRY.PRM"))

#Now let's add a phase from a .cif file to our project
phase0 = gpx.add_phase(os.path.join(datadir,"PbSO4-Wyckoff.cif"),phasename="PbSo4",histograms=[hist1])

#Here is a hack found by the author of this tutorial for setting the number of refinement cycles
gpx.data['Controls']['data']['max cyc'] = 3

#Now we're going to start defining refinement steps
# The first refinement step will be for the background
# This is going to refine the background (assuming a Chebyschev fn?) with 3 terms in the fit
refdict0 = {"set":{"Background":{"no. coeffs":3,"refine":True}},"output":'step4.gpx',"call":HistStats,"callargs":[gpx]}

#now we'll refine the unit cell 
refdict1 = {"set":{"Cell":True},"output":'step5.gpx','call':HistStats} #NOTE: this will set the refinement TRUE for all phases currently in .gpx

#Next is a refinement of the hydrostatic strain
refdict2 = {"set":{"HStrain":True},"histograms":[hist1],"phases":[phase0],"output":'step6.gpx',"call":HistStats}

#Next is a refinement of the microstrain
refdict3 = {"set":{"Mustrain":{"type":"isotropic","refine":True},"Size":{"type":"isotropic","refine":True}},"histograms":[hist1],"output":'step7.gpx',"call":HistStats}

dictList = [refdict1,refdict2,refdict3]
gpx.do_refinements(dictList)

HistStats(gpx)

x = gpx.histogram(0).getdata('x')
y = gpx.histogram(0).getdata('ycalc')
import matplotlib.pyplot as plt
plt.plot(x,y)
plt.show()
