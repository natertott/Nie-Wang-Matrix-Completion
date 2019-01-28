import sys
import csv
import numpy as np
from matplotlib import pyplot as plt
def temp(csv_path):
    hertz = input('Detector Hertz?')
    hertz = float(hertz)

    raw = get_data(csv_path)
    size = raw.shape
    count =-1
    for item in raw[0,:]:
	count += 1
	if item == "1::a":
	    latt = raw[1:size[0]-1,count]
  
    for ii in np.arange(len(latt)):
	if latt[ii] == 'None':
	    latt[ii] = 0
	    ranger = ii    

    strain = np.zeros(len(latt))
    latt = latt.astype(np.float)
    strain = 100*(3.25-latt)/(3.25)

    print strain
 
    roots = np.zeros((len(strain),3))
    for ii in np.arange(len(strain)):
        roots[ii,:] =  np.roots([-2.978*(10**-10), 1.031*(10**-6), 9.472*(10**-4), -0.358- strain[ii]])

#    print roots
    temperature = roots[:,2]
    temperature = temperature - 273
    temperature = abs(temperature)

    #for ii in np.arange(ranger):
    #	temperature[ii] = 0

    plt.plot(np.arange(len(temperature))/hertz,temperature)
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (C)')
    plt.show()
     
def get_data(csv_path):
    results = []
    with open(csv_path) as csvfile:
        reader = csv.reader(csvfile) # change contents to floats
        for row in reader: # each row is a list
            results.append(row)

    return np.array(results)

if __name__ == '__main__':
    temp(sys.argv[1])
