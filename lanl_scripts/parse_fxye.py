import sys
import csv
import glob

# This file takes .fxye files exported from GSASII and re-writes them 
# so that they can be read into GSAS

# Currently not working properly. Something is wrong with the spacing of 
# the output file. 

def parser():

    inputs = glob.glob('*.fxye')
    index = 0

    with open('experiment.fxye','wb') as outfile:
        for f in inputs:
            with open(f,'rb') as infile:
                data = infile.readlines()
		print data
                data[0] = 'test experiment histogram #' +  str(index) + '\n'
                index +=1
                data[1] = 'Instrument parameter file:1ID.prm \n'
                data[2] = 'BANK ' + str(index) + ' 2250 2250 CONS 225.00 0.74 0 0 FXYE \n'
                outfile.writelines(data)
		
if __name__ == '__main__':
    parser()
