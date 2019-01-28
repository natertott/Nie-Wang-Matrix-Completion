import sys
import numpy as np
from matplotlib import pyplot as plt
import csv

def alldata(csv_path):
    hertz = input('Detector Hertz?')
    hertz = float(hertz)

    title = input('Title?')
    raw = get_data(csv_path)
    size = raw.shape

    count = -1
    for item in raw[0,:]:
        count += 1
        if item == "0:*:Scale":
            alpha = raw[2:size[0]-1,count]
        elif item == "1:*:Scale":
            beta = raw[2:size[0]-1,count]
        elif item == ":*:BkPkint;0":
            liq = raw[2:size[0]-1,count]
        elif item == "0::a":
	        alpha_a = raw[2:size[0]-1,count]
        elif item == "1::a":
	        beta_a = raw[2:size[0]-1,count]
        elif item == "0::Vol":
	        alpha_vol = raw[2:size[0]-1,count]
        elif item == "1::Vol":
	        beta_vol = raw[2:size[0]-1,count]

    for ii in np.arange(len(alpha)):
        if alpha[ii] == 'None':
            alpha[ii] = 0.0
    for ii in np.arange(len(beta)):
        if beta[ii] =='None':
            beta[ii] = 0
    for ii in np.arange(len(liq)):
        if liq[ii] == 'None':
            liq[ii] = 0
    for ii in np.arange(len(alpha_a)):
        if alpha_a[ii] == 'None':
            alpha_a[ii] = float(0)
    for ii in np.arange(len(beta_a)):
        if beta_a[ii] == 'None':
            beta_a[ii] = 0
    for ii in np.arange(len(alpha_vol)):
        if alpha_vol[ii] == 'None':
            alpha_vol[ii] = 0
    for ii in np.arange(len(beta_vol)):
        if beta_vol[ii] == 'None':
            beta_vol[ii] = 0

    alpha = alpha.astype(np.float)
    beta = beta.astype(np.float)
    liq = liq.astype(np.float)
    alpha_a = alpha_a.astype(np.float)
    beta_a = beta_a.astype(np.float)
    alpha_vol = alpha_vol.astype(np.float)
    beta_vol = beta_vol.astype(np.float)
    index = np.arange(len(alpha))
    index = index/hertz


# Here I am calculating the phase fraction using Adrian's method
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


#Plotting the phase fraction
    plt.plot(time,alpha_frac, 'r-',time,beta_frac,'g*',time,liq_frac,'b*',time,sum_mat)
    plt.xlabel('Time (sec)')
    plt.ylabel('Phase fraction')
    plt.title(title + 'Relative Phase Fractions')
    plt.legend(['alpha_frac','beta_frac','liq_frac'])
    plt.show()



#Putting the volume fractions into their own arrays
#If I don't fit the volume, it is going to be 0 and throw off the scaling of the graphs
#So I'm only going to select the data that is non-zero

    plot_alpha_vol = []
    alpha_vol_index = []
    for ii in np.arange(len(alpha_vol)):
        if alpha_vol[ii] != 0:
            plot_alpha_vol.append(alpha_vol[ii])
            alpha_vol_index.append(index[ii])

    plot_beta_vol = []
    beta_vol_index = []
    for ii in np.arange(len(beta_vol)):
        if beta_vol[ii] != 0:
            plot_beta_vol.append(beta_vol[ii])
            beta_vol_index.append(index[ii])

# Plotting the volume fraction as a function of time
    plt.plot(alpha_vol_index,plot_alpha_vol,'r*',beta_vol_index,plot_beta_vol,'b.')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Unit Cell Volume (Angstroms^3)')
    plt.legend(['alpha','beta'])
    plt.title(title + ' Unit Cell Volume Change')
    plt.show()

# Putting the lattice parameters into their own (properly formatted) arrays
    plot_alpha_a = []
    alpha_a_index = []
    for ii in np.arange(len(alpha_a)):
        if alpha_a[ii] != 0:
            plot_alpha_a.append(alpha_a[ii])
            alpha_a_index.append(index[ii])

    plot_beta_a = []
    beta_a_index = []
    for ii in np.arange(len(beta_a)):
        if beta_a[ii] != 0:
            plot_beta_a.append(beta_a[ii])
            beta_a_index.append(index[ii])

    # Later we'll need to do some computation with these, so let's convert them to numpy arrays of floats
    plot_alpha_a = np.asarray(plot_alpha_a)
    plot_beta_a = np.asarray(plot_beta_a)

# Plotting the lattice parameters
    plt.plot(alpha_a_index[5:len(alpha_a_index)],plot_alpha_a[5:len(alpha_a_index)],'r*',beta_a_index[5:len(alpha_a_index)],plot_beta_a[5:len(alpha_a_index)],'b.')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Lattice parameter a')
    plt.legend(['alpha a','beta'])
    plt.title(title + ' Lattice Parameter Change')
    plt.show()

#
#Let's calculate some strains!
#The values in beta_a and alpha_a are lattice parameters, adjusted for hydrostatic strains
    # the unstrained alpha a parameter is ~2.95 based on the cif file I am using
    #the unstrained alpha c parameter is ~4.686 based on the cif file I am using
    # the unstrained beta a parameter is ~3.25 based on the cif file I am using

    alpha_a_strain = 100* (plot_alpha_a - 2.860)/(2.860)
    beta_strain = 100* (plot_beta_a -3.652)/(3.652)

    plt.plot(alpha_a_index[5:len(alpha_a_index)],alpha_a_strain[5:len(alpha_a_index)],'r*',beta_a_index[5:len(alpha_a_index)],beta_strain[5:len(alpha_a_index)],'b.')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Eng Strain')
    plt.legend(['alpha a','beta'])
    plt.title(title + ' Lattice Parameter Strains')
    plt.show()


    # This method is from 'Thermophysical Properties of Matter Volume 12: Thermal Expansion, Metallic Elements and Alloys'
    # pg. 1278 has polynomial fits for linear expansion of 304L Stainless Steel with a Composition of (16-19)Cr,(6-12)Ni, remainder Fe and X, as a function of temperature

    # We're gonna compute the temp as a function of both ferrite (alpha) and austenite (beta) strains
    roots_alpha = np.zeros((len(alpha_a_strain),3))
    roots_beta = np.zeros((len(beta_strain),3))
    for ii in np.arange(len(alpha_a_strain)):
	roots_alpha[ii,:] = np.roots([-2.978*(10**-10), 1.181*(10**-6), 9.472*(10**-4), -0.358-alpha_a_strain[ii]])

    for ii in np.arange(len(beta_strain)):
	roots_beta[ii,:] = np.roots([-2.978*(10**-10), 1.181*(10**-6), 9.472*(10**-4), -0.358 - beta_strain[ii]])



    temperature_alpha = roots_alpha[:,2]
    temperature_alpha = temperature_alpha - 273 #the polynomial fit is for Kelvin, this converts to Celsius
    temperature_alpha = abs(temperature_alpha)
    #temperature_beta = roots_beta[:,2]
    #temperature_beta = temperature_beta - 273
    #temperature_beta = abs(temperature_beta)

    #plt.plot(alpha_a_index,temperature_alpha,'r.',beta_a_index,temperature_beta,'b.')
    #plt.xlabel('Time')
    #plt.ylabel('Temperature (C)')
    #plt.legend(['Ferrite Computed Temp','Austenite Comptued Temp'])
    #plt.show()

    liquidus = np.ones(len(temperature_alpha))*1454 #this is the liquidus temperature, roughly
    ferrite_transus = np.ones(len(temperature_alpha))*1424 #ferrite phase transus, roughly
    austenite_solidus = np.ones(len(temperature_alpha))*1300 #austenite phase solidus, roughly
    # I'm going to plot both temperature and phase fractions at the same time, on diffrent axes
    # I don't think the temperature visualization makes much sense w/o also looking at the temperature

    fig, ax1 = plt.subplots()
    ax1.plot(time, alpha_frac, 'r-', time, beta_frac, 'g*', time, liq_frac, 'b*')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Phase fraction')
    ax1.legend(['Ferrite','Austenite','Liquid'])

    ax2 = ax1.twinx()
    ax2.plot(alpha_a_index,temperature_alpha,'m.',alpha_a_index,liquidus,'m',alpha_a_index,ferrite_transus,'m',alpha_a_index,austenite_solidus,'m')
    ax2.set_ylabel('Temperature (deg C)')
    ax2.tick_params('y',colors='m')
    fig.tight_layout()
    plt.title(title)
    plt.show()
    
    # Going to plot strain versus temperature
    plt.plot(temperature_alpha, alpha_a_strain,'r*')
    plt.xlabel('Temperature (C)')
    plt.ylabel('Strain')
    plt.title(title + 'Temperature Dependent Strain ')    
    plt.show() 

def get_data(csv_path):
    results = []
    with open(csv_path) as csvfile:
        reader = csv.reader(csvfile) # change contents to floats
        for row in reader: # each row is a list
            results.append(row)

    return np.array(results)

if __name__ == '__main__':
    alldata(sys.argv[1])
