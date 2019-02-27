import sys,os
import numpy as np
from matplotlib import pyplot as plt
import csv

def is_empty(any_structure):
    if any_structure:
        return 1
    else:
        return 0

# Process all data from a GSASII .csv sequential refinement results file
def get_data(csv_path):
    results = []
    rows = 0
    with open(csv_path) as csvfile:
        reader = csv.reader(csvfile) # change contents to floats
        for row in reader: # each row is a list
            rows = rows+1
            results.append(row)

    return np.asarray(results)

# Process all raw data returned from get_data()
def alldata(csv_path,hertz):
    hertz = float(hertz)

    raw = get_data(csv_path)
    size = raw.shape

    # This is the dictionary we will pass to all subsequent functions
    # All data pulled from GSASII will be stored here

    data = {'hertz':hertz}

    #This is going to store the raw data
    data_dict = {}

    count = -1
    for item in raw[0,:]:
        count += 1
        if item == "0:*:Scale":
            data_dict["alpha"] = raw[2:size[0]-1,count]
        elif item == "1:*:Scale":
            data_dict["beta"] = raw[2:size[0]-1,count]
        elif item == ":*:BkPkint;0":
            data_dict["liq"] = raw[2:size[0]-1,count]
        elif item == "0::a":
            data_dict["alpha_a"] = raw[2:size[0]-1,count]
        elif item == "0::c":
            data_dict["alpha_c"] = raw[2:size[0]-1,count]
        elif item == "1::a":
            data_dict["beta_a"] = raw[2:size[0]-1,count]
        elif item == "0::Vol":
            data_dict["alpha_vol"] = raw[2:size[0]-1,count]
        elif item == "1::Vol":
            data_dict["beta_vol"] = raw[2:size[0]-1,count]

    for k, v in data_dict.items():
        for ii in np.arange(0,len(v)):
            if v[ii] == 'None':
                data_dict[k][ii] = 0.0
        data_dict[k] = data_dict[k].astype(float)

    index = np.asarray(np.arange(len(data_dict["alpha"])))
    data_dict["index"] = index

    #Store the raw values
    data['raw'] = data_dict


    """At this point, every data list is equally long
    Some of the values in those lists are 0 (previously 'None')
    We need to keep track of events that happen for subsets of the data and
    when then happen

    Many of entries in the lattice paramater and volume lists are 0
    We need to pull out the entries which are not 0 and note their associated index
    This will inform of us when these lists start"""

    # I am going to store the lattice parameter, volume values in this dict
    lattice = {}

    for k,v in data_dict.items():
        # We only want the lattice p's and volume
        if k == 'index' or k == 'liq' or k == 'alpha' or k == 'beta':
            continue
        else:
            v = v.astype(float)
            store = np.nonzero(v) #Find the index of all values not equal to
            lattice[k] = {k:data_dict[k][store].astype(float),k + str('_index'):data_dict['index'][store]} #store the values

    #Update the data dict with our selected lattice parameters
    data['lattice'] = lattice
    return data

########################################
def phasefrac(data):
    # This is a function to calculate phase fractions (if you want to call it that)
    # It takes as input the full data dict returned from alldata()

    #Hacks. Hacks everywhere.
    #First let's pull the values we need out of the data dictionary

    pf = data['raw']
    alpha = pf['alpha']
    beta = pf['beta']

    #check if the liquid phase is present
    hack = 0
    if "liq" in pf.keys():
        liq = pf['liq']
        hack = 1

    #For calculating time from indices
    hertz = data['hertz']

    #Store the length of the arrays
    size = alpha.size

    # Here I am calculating the phase fraction using Adrian's method
    for ii in np.arange(len(alpha)):
        if alpha[ii] < 0:
                alpha[ii] = 0

    for ii in np.arange(len(beta)):
       if beta[ii] < 0:
               beta[ii] = 0

    if hack == 1:
        for ii in np.arange(len(liq)):
            if liq[ii] < 0:
                liq[ii] = 0
        #This is a 'fudge factor' which we will use to scale the background peak
        # intensity fit of the liquid amorphous scatter.

        # The maximum intensity of the liquid phase will be scaled by the sum of the peak
        # intensity for alpha + beta in the final frame
        fudge = alpha[size-1]+beta[size-1]
        liq = (liq/max(liq))*fudge

    alpha_frac = np.zeros(len(alpha))
    beta_frac = np.zeros(len(beta))
    liq_frac = np.zeros(len(liq))

    # I'm going to try and get rid of some of the noise that is present in the pre-welder region
    for ii in np.arange(len(alpha)):
        if alpha[ii] + beta[ii] <0.00001:
            alpha[ii] = 0.000001
            beta[ii] = 0.000001
            liq[ii] = 100

    for i in np.arange(len(alpha)):
        alpha_frac[i] = alpha[i]/(alpha[i] + beta[i] + liq[i])
        beta_frac[i] = beta[i]/(alpha[i] + beta[i] + liq[i])
        liq_frac[i] = liq[i]/(alpha[i] + beta[i] + liq[i])

    sum_mat = alpha_frac + liq_frac + beta_frac
    time = np.arange(0,len(alpha_frac))
    time = time/hertz
    pf_dict = {'alpha_frac':alpha_frac,'beta_frac':beta_frac,'liq_frac':liq_frac,'time':time}
    return pf_dict

######################################################
def strain(data):
    # This function calculates the strain on each phase's lattice
    lattice = data['lattice']
    strain = {}

    for k,v in lattice.items():
        # Each key in lattice is tied to another dict
        # These dicts should always have two keys themselves:
        #   - the lattice parameters
        #   - the associated index
        # They should always be in that order too

        # I'll pull out the names of each key in v
        # Then use those names to pull the phase and index, respectively
        names = list(v.keys())
        phase = v[names[0]]
        index = v[names[1]]

        #Define a datum to measure strain
        datum = phase[len(phase)-1]
        if datum != 0:
            eps = (phase - datum)/datum

        # Store the strain values under a key of the same name as the key
        # in the lattice dict
        strain[names[0]] = {k+ '_strain':eps, k + '_index':index}

    return strain
