import sys
import numpy as np
import math
import pandas
import glob
from matplotlib import pyplot as plt

# Why is it so difficult to comment your code, Nathan?
# This code plots the intensity of diffraction peaks as a function of time
# It operates on integrated diffraction spectra (I vs 2theta, for example)
# and it plots subsequent spectra such that evolution in peak intensity or 
# position can be seen through time

def plot_intensities():

    hertz = input('Detector Hertz?')
    hertz = float(hertz)
    wavey = input('Wavelength (nm)?')
    
    x = []
    intensity =[]

    colnames= ['x','intensity','garbage','garbage','garbage']
    for file_path in sorted(glob.glob('*.csv')):
        with open(file_path,'r') as powder:
            data = pandas.read_csv(powder,names=colnames)
	    store_x = data.x
	    store_int = data.intensity
	    x.append(store_x[13:len(store_x)])
	    intensity.append(store_int[13:len(store_int)])

    #Turn the arays of 2-theta (a.k.a. x) and intensity into numpy arrays
    # The arrays are arranged as follows: each row is a 2-theta position for a given histogram, each column is a new histogram
    intensity = np.c_[intensity]
    intensity = np.transpose(intensity)
    intensity = intensity.astype(np.float)
    x = np.c_[x]
    x = x.astype(np.float)

    d = x[1,:]
    d = (d/180)*math.pi
    lower = (wavey)/(2*math.sin(d[1]))
    upper = (wavey)/(2*math.sin(d[len(d)-1]))	 
		
    fig, ax = plt.subplots()
    im = ax.imshow(intensity,cmap='gist_stern',extent = [0,1/hertz*len(d),upper,lower],aspect = 15)
    plt.xlabel('Time (s)')
    plt.ylabel('d Spacing (Ang.)')
    plt.show()

if __name__ == '__main__':
    plot_intensities()
