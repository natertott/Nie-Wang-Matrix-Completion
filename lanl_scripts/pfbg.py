#making this script executable
#!/usr/bin/env python


import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

# This script is for pulling phase fraction calculations out of GSASII .lst output files
# Since GSASII gives somewhat nonsensical phase fraction values, on occassion, this script will normalize all values to add up to 1
# It is currently written to only handle Ti-6Al-4V scans that have alpha, beta, and liquid phase fractions in them

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
            alpha = raw[2:size[0]-1,count]
        elif item == "1:*:Scale":
            beta = raw[2:size[0]-1,count]
        elif item == ":*:BkPkint;0":
            liq = raw[2:size[0]-1,count]
    for ii in np.arange(len(liq)):
	if liq[ii] == 'None':
		liq[ii] = 0
    for ii in np.arange(len(alpha)):
        if alpha[ii] == 'None':
            alpha[ii] = 0
    for ii in np.arange(len(beta)):
        if beta[ii] == 'None':
            beta[ii] = 0

    

    alpha = alpha.astype(np.float)
    beta = beta.astype(np.float)
    liq = liq.astype(np.float)
    for ii in np.arange(len(alpha)):
        if alpha[ii] < 0:
                alpha[ii] = 0

    for ii in np.arange(len(beta)):
       if beta[ii] < 0:
               beta[ii] = 0

    for ii in np.arange(len(liq)):
        if liq[ii] < 0:
                liq[ii] = 0


    fudge = alpha[size[0]-4]+beta[size[0]-4]
    liq = (liq/max(liq))*fudge

    alpha_frac = np.zeros(len(alpha))
    beta_frac = np.zeros(len(beta))
    liq_frac = np.zeros(len(liq))

    for i in np.arange(ranger,len(alpha)):
        alpha_frac[i] = alpha[i]/(alpha[i] + beta[i] + liq[i])
        beta_frac[i] = beta[i]/(alpha[i] + beta[i] + liq[i])
        liq_frac[i] = liq[i]/(alpha[i] + beta[i] + liq[i])

    sum_mat = alpha_frac + liq_frac + beta_frac

    time = np.arange(0,len(alpha_frac))
    time = time/hertz


    plt.plot(time,alpha_frac, 'r-',time,beta_frac,'g*',time,liq_frac,'b*',time,sum_mat)
    plt.xlabel('Time (sec)')
    plt.ylabel('Phase fraction')
    plt.title(csv_path)
    plt.legend(['alpha_frac','beta_frac','liq_frac'])
    plt.show()

    if hertz == 10:
        for ii in np.arange(len(time)):
            if time[ii] == 2.6:
                count = ii

        plt.plot(time[0:count],alpha_frac[0:count],'r-',time[0:count],beta_frac[0:count],'b.',time[0:count],liq_frac[0:count],'g*',time[0:count],sum_mat[0:count])
        plt.xlabel('Time (sec)')
        plt.ylabel('Phase fraction')
        plt.title(csv_path)
        plt.legend(['alpha_frac','beta_frac','liq_frac'])
        plt.show()

def get_data(csv_path):
    results = []
    with open(csv_path) as csvfile:
        reader = csv.reader(csvfile) # change contents to floats
        for row in reader: # each row is a list
            results.append(row)

    return np.array(results)

if __name__ == "__main__":
    pfvt(sys.argv[1])
