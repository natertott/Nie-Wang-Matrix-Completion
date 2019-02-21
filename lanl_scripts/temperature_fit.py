import sys,os
sys.path.insert(0,os.path.expanduser("~/git/lanl_scripts/"))
import tiAnalysis as ti
import numpy as np
from matplotlib import pyplot as plt

def main(csv_path,title,hertz):
    # First we'll do the basic data processing to get the values we need
    data = ti.alldata(csv_path,hertz)
    lattice = data['lattice']
    strain = ti.strain(data)

    # Going to need the raw beta lattice parameter Later
    beta_a = lattice['beta_a']['beta_a']

    # Pull out the lattice strain values for beta
    beta = strain['beta_a']
    beta_strain = beta['beta_a_strain']
    beta_index = beta['beta_a_index']

    #### Calculating temperature
    #Elmer and Palmer report in Mat. Sci & Eng. A, 391, (2005), 104-113 that the coefficient of linear thermal expansion for the a parameter in alpha Ti is alpha_v = 9.7*1-^-6
    coef1 = 5.8 * (10**-5)
    coef2 = 9.2 * (10**-6)
    beta_temp = np.zeros(len(beta_strain))
    count = 0
    for nn in beta_a:
        if nn > 3.26:
            beta_temp[count] = beta_strain[count]/coef1
            count = count+1
        else:
            beta_temp[count] = beta_strain[count]/coef2
            count = count+1

    for nn in np.arange(len(beta_temp)):
        if beta_temp[nn] < 0:
            beta_temp[nn] = 0

    #Now we'll try to do the same thing using the alpha lattice parameter
    coef = 9.7 * (10**-6)
    alpha = strain['alpha_a']
    alpha_strain = alpha['alpha_a_strain']
    alpha_temp = np.zeros(len(alpha_strain))
    alpha_index = alpha['alpha_a_index']
    alpha_temp = alpha_strain/coef

    plt.plot(alpha_index,alpha_temp,'r.',beta_index,beta_temp,'b.')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.legend(['alpha','beta'])
    plt.savefig('temperature.png')

    ### Calculating macroscopic strain fit
    # I am going to regress the macroscopic strain polynomial on a vector of the beta and alpha strains at different times

    # # This method is from 'Thermophysical Properties of Matter Volume 12: Thermal Expansion, Metallic Elements and Alloys'
    # # pg. 1278 has polynomial fits for linear expansion of Ti-6Al-4V alloys as a function of temperature
    #roots = np.zeros((len(strain),3))
    #for ii in np.arange(len(strains)):
    #     roots[ii,:] = np.roots([-1.994*(10**-10), 5.807*(10**-7), 5.992*(10**-4), -0.22 - beta_strain[ii]])

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3])
