import sys
import os
import numpy as np
from matplotlib import pyplot as plt
import diffcomp
import read_ge

def image(ceria,dark):
    ims = read_ge.read_ge(ceria)
    darks = read_ge.read_ge(dark)
    
    #Only going to take the first image
    im  = ims[0,:,:]
    dark = darks[0,:,:]
    im = im - dark
    
    sum = np.sum(im,axis=0)
    plt.plot(sum)
    plt.show()




if __name__ == '__main__':
	image(sys.argv[1],sys.argv[2])
