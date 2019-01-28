#making this script executable
#!/usr/bin/env python


import sys
import numpy as np
import matplotlib.pyplot as plt
import re

# This script is for pulling phase fraction calculations out of GSASII .lst output files
# Since GSASII gives somewhat nonsensical phase fraction values, on occassion, this script will normalize all values to add up to 1
# It is currently written to only handle Ti-6Al-4V scans that have alpha, beta, and liquid phase fractions in them


def pfvt(lst_path):
    raw = get_ti_fracs(lst_path)
    split = len(raw)/2

    alpha_frac = np.array(raw[0:split])
    beta_frac = np.array(raw[split:len(raw)])

    for num in np.arange(len(alpha_frac)):
        if alpha_frac[num] < 0:
            alpha_frac[num] = 0

    #hold = np.arange(len(alpha_frac))
    #plt.plot(hold,alpha_frac,hold,beta_frac)
    #plt.show()

    for num in np.arange(len(alpha_frac)):
        if alpha_frac[num] < 0:
            alpha_frac[num] = 0

    for num in np.arange(len(beta_frac)):
        if beta_frac[num] < 0:
            beta_frac[num] = 0

    #beta_frac = beta_frac/max(beta_frac)

    max_alpha = alpha_frac[len(alpha_frac)-1]
    max_beta = beta_frac[len(beta_frac)-1]



    alpha_frac = alpha_frac/(max_alpha + max_beta)
    beta_frac = beta_frac/(max_alpha + max_beta)

    for num in np.arange(len(alpha_frac)):
        if alpha_frac[num] > 1:
            alpha_frac[num] = 1

    for num in np.arange(len(beta_frac)):
        if beta_frac[num] > 1:
            beta_frac[num] = 1

    liq_frac = np.ones(len(alpha_frac)) - alpha_frac - beta_frac

    #return alpha_frac, beta_frac, liq_frac

    # for num in np.arange(len(liq_frac)):
    #     if liq_frac[num] < 0:
    #         liq_frac[num] = 0

    sum_mat = alpha_frac + liq_frac + beta_frac

    det = input("Detector Hertz?")
    if det == "10":
        time = np.arange(0,len(alpha_frac))
        time = time/20.
    elif det == "200":
        time = np.arange(0,len(alpha_frac))
        time = time/200.


    plt.plot(time,alpha_frac, 'r-',time,beta_frac,'g-',time,liq_frac,'b.',time,sum_mat)
    plt.xlabel('Time (sec)')
    plt.ylabel('Phase fraction')
    plt.title(lst_path)
    plt.legend(['alpha_frac','beta_frac','liq_frac'])
    plt.show()

    return alpha_frac, beta_frac, liq_frac

def get_ti_fracs(lst_path):
    f = open(lst_path,"r").read()

    where = re.compile('Phase\sfraction\s\s:[\s\-0-9.0-9]*\s')
    matches = where.findall(f)
    matches = ''.join(matches)

    pf = re.compile('[\-|\d][\d|\.][\.0-9]*\d')
    frac = pf.findall(matches)

    raw = np.asarray(frac,dtype='float')
    return raw

if __name__ == "__main__":
    pfvt(sys.argv[1])
