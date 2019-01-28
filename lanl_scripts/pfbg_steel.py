#making this script executable
#!/usr/bin/env python


import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

# This script is for pulling phase fraction calculations out of GSASII .lst output files
# Since GSASII gives somewhat nonsensical phase fraction values, on occassion, this script will normalize all values to add up to 1
# It is currently written to only handle 304L scans that have ferrite, austenite and liquid phase fractions in them

def pfvt(csv_path):

    hertz = input('Detector Hertz?')
    hertz = float(hertz)
    raw = get_data(csv_path)
    size = raw.shape
    ranger = 0

    count = -1
    for item in raw[0,:]:
        count += 1
        if item == "0:*:Scale":
            ferrite = raw[2:size[0]-1,count]
        elif item == "1:*:Scale":
            austenite = raw[2:size[0]-1,count]
        elif item == ":*:BkPkint;0":
            liq = raw[2:size[0]-1,count]
    for ii in np.arange(len(liq)):
	if liq[ii] == 'None':
		liq[ii] = 0
    for ii in np.arange(len(ferrite)):
        if ferrite[ii] == 'None':
            ferrite[ii] = 0
    for ii in np.arange(len(austenite)):
        if austenite[ii] == 'None':
            ferrite[ii] = 0

    

    ferrite = ferrite.astype(np.float)
    austenite = austenite.astype(np.float)
    liq = liq.astype(np.float)
    for ii in np.arange(len(ferrite)):
        if ferrite[ii] < 0:
                ferrite[ii] = 0

    for ii in np.arange(len(austenite)):
       if austenite[ii] < 0:
               austenite[ii] = 0

    for ii in np.arange(len(liq)):
        if liq[ii] < 0:
                liq[ii] = 0


    fudge = ferrite[size[0]-4]+austenite[size[0]-4]
    liq = (liq/max(liq))*fudge

    ferrite_frac = np.zeros(len(ferrite))
    austenite_frac = np.zeros(len(austenite))
    liq_frac = np.zeros(len(liq))

    for i in np.arange(ranger,len(ferrite)):
        ferrite_frac[i] = ferrite[i]/(ferrite[i] + austenite[i] + liq[i])
        austenite_frac[i] = austenite[i]/(ferrite[i] + austenite[i] + liq[i])
        liq_frac[i] = liq[i]/(ferrite[i] + austenite[i] + liq[i])

    sum_mat = ferrite_frac + austenite_frac + liq_frac

    time = np.arange(0,len(ferrite_frac))
    time = time/hertz


    plt.plot(time,ferrite_frac, 'r-',time,austenite_frac,'g*',time,liq_frac,'b*',time,sum_mat)
    plt.xlabel('Time (sec)')
    plt.ylabel('Phase fraction')
    plt.title(csv_path)
    plt.legend(['ferrite_frac','austenite_frac','liq_frac'])
    plt.show()

#    if hertz == 10:
#        for ii in np.arange(len(time)):
#            if time[ii] == 2.6:
#                count = ii

#        plt.plot(time[0:count],alpha_frac[0:count],'r-',time[0:count],beta_frac[0:count],'b.',time[0:count],liq_frac[0:count],'g*',time[0:count],sum_mat[0:count])
#        plt.xlabel('Time (sec)')
#        plt.ylabel('Phase fraction')
#        plt.title(csv_path)
#        plt.legend(['alpha_frac','beta_frac','liq_frac'])
#        plt.show()

def get_data(csv_path):
    results = []
    with open(csv_path) as csvfile:
        reader = csv.reader(csvfile) # change contents to floats
        for row in reader: # each row is a list
            results.append(row)

    return np.array(results)

if __name__ == "__main__":
    pfvt(sys.argv[1])
