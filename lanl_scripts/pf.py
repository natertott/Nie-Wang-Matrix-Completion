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
    raw = get_data(csv_path)
    size = raw.shape

    count = -1
    for item in raw[0,:]:
        count += 1
        if item == "0:*:Scale":
            alpha_frac = raw[1:size[0]-1,count]
        elif item == "1:*:Scale":
            beta_frac = raw[1:size[0]-1,count]

    alpha_frac = alpha_frac.astype(np.float)
    beta_frac = beta_frac.astype(np.float)

    max_alpha = alpha_frac[size[0]-3]
    max_beta = beta_frac[size[0]-3]

    for ii in np.arange(len(alpha_frac)):
        if alpha_frac[ii] < 0:
            alpha_frac[ii] = 0

    for ii in np.arange(len(beta_frac)):
        if beta_frac[ii] < 0:
            beta_frac[ii] = 0

    alpha_frac = alpha_frac/(max_alpha + max_beta)
    beta_frac = beta_frac/(max_alpha + max_beta)


    liq_frac = np.ones(len(alpha_frac)) - alpha_frac - beta_frac

    for i in np.arange(len(alpha_frac)):
        if alpha_frac[i] > 1:
            alpha_frac[i] = 1

    for ii in np.arange(len(beta_frac)):
        if beta_frac[ii] > 1:
            beta_frac[ii] = 1

    for ii in np.arange(len(liq_frac)):
        if liq_frac[ii] > 1:
            liq_frac[ii] == 1

    for ii in np.arange(len(liq_frac)):
        if liq_frac[ii] < 0:
            liq_frac[ii] = 0

    sum_mat = alpha_frac + liq_frac + beta_frac

    time = np.arange(0,len(alpha_frac))
    time = time/10.


    plt.plot(time,alpha_frac, 'r-',time,beta_frac,'g*',time,liq_frac,'b*',time,sum_mat)
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
