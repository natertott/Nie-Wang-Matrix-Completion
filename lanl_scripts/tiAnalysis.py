import sys
import numpy as np
from matplotlib import pyplot as plt
import csv

def alldata(csv_path,hertz,name):
    hertz = float(hertz)

    title = str(name)
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
        elif item == "0::c":
            alpha_c = raw[2:size[0]-1,count]
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
    for ii in np.arange(len(alpha_c)):
        if alpha_c[ii] == 'None':
            alpha_c[ii] = float(0)
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
    alpha_c = alpha_c.astype(np.float)
    beta_a = beta_a.astype(np.float)
    alpha_vol = alpha_vol.astype(np.float)
    beta_vol = beta_vol.astype(np.float)
    index = np.arange(len(alpha))
    index = index/hertz

    data_dict = {'alpha':alpha,'beta':beta,'liq':liq,'alpha_a':alpha_a,'alpha_c':alpha_c,'beta_a':beta_a,'alpha_vol':alpha_vol,'beta_vol':beta_vol,'index':index,'hertz':hertz}
    return data_dict

########################################
def phasefrac(data_dict):
    #Hacks. Hacks everywhere.
    #First let's pull the values we need out of the data dictionary

    alpha = data_dict['alpha']
    beta = data_dict['beta']
    liq = data_dict['liq']
    hertz = data_dict['hertz']

    #Store the length of the arrays
    size = alpha.size

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
    fudge = alpha[size-4]+beta[size-4]
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

###################################################
def lat_param(data_dict):
    # First we'll read in the data that we need from the output of alldata()
    # Then we'll play this game of getting it into numpy arrays
    alpha_a = data_dict['alpha_a']
    alpha_c = data_dict['alpha_c']
    beta_a = data_dict['beta_a']
    index = data_dict['index']

    # Putting the lattice parameters into their own (properly formatted) arrays
    plot_alpha_a = []
    alpha_a_index = []
    for ii in np.arange(len(alpha_a)):
        if alpha_a[ii] != 0:
            plot_alpha_a.append(alpha_a[ii])
            alpha_a_index.append(index[ii])

    plot_alpha_c = []
    alpha_c_index = []
    for ii in np.arange(len(alpha_c)):
        if alpha_c[ii] != 0:
            plot_alpha_c.append(alpha_c[ii])
            alpha_c_index.append(index[ii])

    plot_beta_a = []
    beta_a_index = []
    for ii in np.arange(len(beta_a)):
        if beta_a[ii] != 0:
            plot_beta_a.append(beta_a[ii])
            beta_a_index.append(index[ii])

    # Later we'll need to do some computation with these, so let's convert them to numpy arrays of floats
    plot_alpha_a = np.asarray(plot_alpha_a)
    plot_alpha_c = np.asarray(plot_alpha_c)
    plot_beta_a = np.asarray(plot_beta_a)
    alpha_a_index = np.asarray(alpha_a_index)
    alpha_c_index = np.asarray(alpha_c_index)
    beta_a_index = np.asarray(beta_a_index)

    plot_param_dict = {'plotalpha':{'a':plot_alpha_a,'aindex':alpha_a_index,'c':plot_alpha_c,'cindex':alpha_c_index},'plotbeta':{'a':plot_beta_a,'aindex':beta_a_index}}
    return plot_param_dict

######################################################
def strain(plot_param_dict):
    plot_alpha_a = plot_param_dict['plotalpha']['a']
    alpha_a_index = plot_param_dict['plotalpha']['aindex']
    plot_alpha_c = plot_param_dict['plotalpha']['c']
    alpha_c_index = plot_param_dict['plotalpha']['cindex']
    plot_beta_a = plot_param_dict['plotbeta']['a']
    beta_a_index = plot_param_dict['plotbeta']['aindex']

    #Let's calculate some strains!
    #The values in beta_a and alpha_a are lattice parameters, adjusted for hydrostatic strains
    aa_datum = plot_alpha_a[len(plot_alpha_a)-1]
    ac_datum = plot_alpha_c[len(plot_alpha_c)-1]
    ba_datum = plot_beta_a[len(plot_beta_a)-1]

    alpha_a_strain = (plot_alpha_a - aa_datum)/(aa_datum)
    alpha_c_strain = (plot_alpha_c- ac_datum)/(ac_datum)
    beta_strain =(plot_beta_a - ba_datum)/(ba_datum)
    strain_dict = {'alpha':{'a':alpha_a_strain,'c':alpha_c_strain},'beta':beta_strain}
    return strain_dict

####################################################
def tiTemp(strain_dict,plot_param_dict):
    #Reading in the necessary data to perform the computation
    beta_a = plot_param_dict['plotbeta']['a']
    beta_strain = strain_dict['beta']
    index = plot_param_dict['plotbeta']['aindex']

    #### Calculating temperature
    #Elmer and Palmer report in Mat. Sci & Eng. A, 391, (2005), 104-113 that the coefficient of linear thermal expansion for the a parameter in alpha Ti is alpha_v = 9.7*1-^-6
    coef1 = 5.8 * (10**-6)
    coef2 = 5.8 * (10**-6)
    print(beta_strain[10])
    print(beta_a[10])
    beta_temp = np.zeros(len(beta_strain))
    for nn in np.arange(len(beta_strain)):
        if beta_a[nn] > 3.26:
            beta_temp[nn] = beta_strain[nn]/coef2
        else:
            beta_temp[nn] = beta_strain[nn]/coef1

    plt.plot(beta_temp,'r*')
    plt.show()

    # # This method is from 'Thermophysical Properties of Matter Volume 12: Thermal Expansion, Metallic Elements and Alloys'
    # # pg. 1278 has polynomial fits for linear expansion of Ti-6Al-4V alloys as a function of temperature
    # roots = np.zeros((len(beta_strain),3))
    # for ii in np.arange(len(beta_strain)):
    #     roots[ii,:] = np.roots([-1.994*(10**-10), 5.807*(10**-7), 5.992*(10**-4), -0.22 - beta_strain[ii]])
    #
    # temperature = roots[:,2]
    # temperature = temperature - 273 #the polynomial fit is for Kelvin, this converts to Celsius
    # temperature = abs(temperature)
    # for ii in np.arange(len(temperature)):
    #     if temperature[ii] > 1600:
    #         temperature[ii] = 1600
    #
    # transus = np.ones(len(temperature))*980 #this is the beta transus temperature
    # martensite_s = np.ones(len(temperature))*780 #martensite start temperature
    # martensite_f = np.ones(len(temperature))*650
    # # I'm going to plot both temperature and phase fractions at the same time, on diffrent axes
    # # I don't think the temperature visualization makes much sense w/o also looking at the temperature
    #
    # fig, ax1 = plt.subplots()
    # ax1.plot(time, alpha_frac, 'r-', time, beta_frac, 'g*', time, liq_frac, 'b*')
    # ax1.set_xlabel('Time (s)')
    # ax1.set_ylabel('Phase fraction')
    # ax1.legend(['Alpha','Beta','Liquid'])
    #
    # ax2 = ax1.twinx()
    # ax2.plot(beta_a_index,temperature,'m.',beta_a_index,transus,'m',beta_a_index,martensite_s,'m',beta_a_index,martensite_f,'m')
    # ax2.set_ylabel('Temperature (deg C)')
    # ax2.tick_params('y',colors='m')
    # fig.tight_layout()
    # plt.title(title)
    # plt.show()

###########################################
def pf_plot(pf_dict,title):
    #Plotting the phase fraction
    #Takes as input the dictionary returned by phasefrac()
    alpha_frac = pf_dict['alpha_frac']
    beta_frac = pf_dict['beta_frac']
    liq_frac = pf_dict['liq_frac']
    time = pf_dict['time']

    #Now for plotting
    plt.plot(time,alpha_frac, 'r-',time,beta_frac,'g*',time,liq_frac,'b*')
    plt.xlabel('Time (sec)')
    plt.ylabel('Phase fraction')
    plt.title(title + 'Relative Phase Fractions')
    plt.legend(['Alpha Phase','Beta Phase','Liquid'],loc=2)
    plt.show()

###################################################
def lat_param_plot(plot_param_dict,title):
    #Read in all the data first
    plot_alpha_a = plot_param_dict['plotalpha']['a']
    alpha_a_index = plot_param_dict['plotalpha']['aindex']
    plot_alpha_c = plot_param_dict['plotalpha']['c']
    alpha_c_index = plot_param_dict['plotalpha']['cindex']
    plot_beta_a = plot_param_dict['plotbeta']['a']
    beta_a_index = plot_param_dict['plotbeta']['aindex']

    # Plotting the lattice parameters
    plt.plot(alpha_a_index,plot_alpha_a,'r*',alpha_c_index,plot_alpha_c,'c+',beta_a_index,plot_beta_a,'b.')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Lattice parameter a')
    plt.legend(['alpha a','alpha c','beta'])
    plt.title(title + ' Lattice Parameters')
    plt.show()

    #Why won't this plot???
    plt.plot(alpha_a_index,plot_alpha_c/plot_alpha_a,'r*')
    plt.xlabel('Time (seconds)')
    plt.ylabel('c//a')
    plt.legend('alpha')
    plt.title(title + ' Lattice Parameter Ratio for alpha Phase')
    plt.show()

##################################################
def vf_plot(data_dict,title):
    #Pulling out the data values that we need
    alpha_vol = data_dict['alpha_vol']
    beta_vol = data_dict['beta_vol']
    index = data_dict['index']

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

def get_data(csv_path):
    results = []
    with open(csv_path) as csvfile:
        reader = csv.reader(csvfile) # change contents to floats
        for row in reader: # each row is a list
            results.append(row)

    return np.array(results)
