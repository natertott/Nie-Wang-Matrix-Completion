import sys
import numpy as np
from glob import glob
from matplotlib import pyplot as plt
from PIL import Image

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


def showorsave(ge_path):
    """
    Will either show or save (or both) .ge files obtained from the 1ID Beamline     at  the Advanced Photon Source at Argonne National Laboratory.

    :param outfile: the filename of the image you want to save
    :param sors: whether or not the user wants to show, save, or both
    
    Code will either display each frame of a .ge file using plt.imshow()
    or it will save the images using plt.savefig
    or both.
    """

    outfile = input('Output Filename?')
    sors = input('Show(1), Save(2), or Both(3)?')
    i = 0
    if sors ==1:
	for data in read_ge(ge_path):
    	    print i
    	    plt.imshow(data)
    	    plt.show()
    	    i+=1
    if sors ==2:
        for data in read_ge(ge_path):
            print i
	    cmap = plt.cm.bone
	    image = cmap(data)
	    plt.imsave(outfile + str(i) + '.png',image)
            i+=1
    if sors ==3:
        for data in read_ge(ge_path):
            print i
            plt.imshow(data)
            cmap = plt.cm.bone
	    image = cmap(data)
	    plt.savefig(outfile + str(i) + '.png',image)
            i+=1
    else:
	print('Not a valid option')


    
# combine the data from a list of files.
def join_ge(filelist):
    return np.concatenate([d for d in ge_iterator(filelist)])

if __name__ == '__main__':
   showorsave(sys.argv[1]) 
