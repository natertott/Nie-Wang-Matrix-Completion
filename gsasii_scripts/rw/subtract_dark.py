##########################
# Nathan S Johnson
# Python script for subtracting dark files from diffraction images
# 11/12/2018
########################

import sys
import numpy as np
import scipy.io
from matplotlib import pyplot as plt

def read_ge(filename):
    """
    Reads the binary output from GE detectors.

    :param filename: The name of the file that contains the dataself.
    :type filename: str
    :returns: Frames contained in the GE file.
    :rtype: numpy.ndarray
    """   
    with open(filename) as ifs:
        ifs.seek(8192)
        data = np.fromfile(ifs, np.uint16)
    return np.reshape(data, (-1, 2048, 2048))

def subtract(diff,dark):
    diffall = read_ge(diff)
    
    # The read_ge function outputs all the stacked images into stacked arrays
    # We only really want the first array
    diffimg = diffall[0,:,:]
    np.shape(diffimg)   

    darkall = read_ge(dark)
    darkimg = darkall[0,:,:]
    thresh = diffimg - darkimg
    scipy.io.savemat('tresh.mat', mdict={'thresh':thresh})
    #plt.imshow(thresh)
    #plt.show() 

    #plt.plot(thresh)
    #plt.show()	
    


if __name__ == '__main__':
	subtract(sys.argv[1],sys.argv[2])
