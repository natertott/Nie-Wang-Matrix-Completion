import sys,os
sys.path.insert(0,os.path.expanduser("~/git/lanl_scripts/"))
import tiAnalysis as ti
import numpy as np
from matplotlib import pyplot as plt

def main(csv_path,title,hertz):
    hertz = float(hertz)

    # First we'll do the basic data processing to get the values we need
    data = ti.alldata(csv_path,hertz)
    lattice = data['lattice']
    strain = ti.strain(data)

    raw = data['raw']
    alpha_a = raw['alpha_a']
    beta_a = raw['beta_a']

    #plt.plot(np.true_divide(np.arange(len(alpha_a),hertz)))

    # Might as well get the relative phase fractions
    phasefrac = ti.phasefrac(data)
    alpha = phasefrac['alpha_frac']
    beta = phasefrac['beta_frac']
    liq = phasefrac['liq_frac']
    time = phasefrac['time']

    plt.figure()
    plt.plot(time,liq,'k.',time,alpha,'r.',time,beta,'b.')
    plt.xlabel('Time (s)')
    plt.ylabel('Fraction')
    plt.legend(['Liquid','Alpha','Beta'])
    plt.savefig(str(title) + '_pf.png')

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


    #######
    # This one is going to be kind of weird
    # Let's take the first n strain values of beta and the last m strain values from alpha
    # Then we'll find the roots of a polynomial fit to that
    # Then plot the temperature

    bs = beta_strain[13:103]
    bi = beta_index[13:103]
    alphs = alpha_strain[len(alpha_strain)-100:len(alpha_strain)]
    alphsi = alpha_index[len(alpha_index)-100:len(alpha_index)]
    macrostrain = np.zeros([len(bs) + len(alphs)])
    macroindex = np.zeros([len(bi) + len(alphsi)])
    macrostrain[0:len(bs)] = bs
    macrostrain[len(bs):len(macrostrain)] = alphs
    macroindex[0:len(bi)] = bi
    macroindex[len(bi):len(macroindex)] = alphsi


    # # This method is from 'Thermophysical Properties of Matter Volume 12: Thermal Expansion, Metallic Elements and Alloys'
    # # pg. 1278 has polynomial fits for linear expansion of Ti-6Al-4V alloys as a function of temperature
    roots = np.zeros((len(macrostrain),3))
    for ii in np.arange(len(macrostrain)):
         roots[ii,:] = np.roots([-1.994*(10**-10), 5.807*(10**-7), 5.992*(10**-4), -0.22 - macrostrain[ii]])

    # Now we'll try the same thing for just beta strain, just alpha strain
    beta_roots = np.zeros((len(beta_strain),3))
    for ii in np.arange(len(beta_strain)):
        beta_roots[ii,:] = np.roots([-1.994*(10**-10), 5.807*(10**-7), 5.992*(10**-4), -0.22 - beta_strain[ii]])

    alpha_roots = np.zeros((len(alpha_strain),3))
    for ii in np.arange(len(alpha_strain)):
        alpha_roots[ii,:] = np.roots([-1.994*(10**-10), 5.807*(10**-7), 5.992*(10**-4), -0.22 - alpha_strain[ii]])

    ### Let's plot up the polynomial we're fitting
    temps = np.arange(1,2000)
    strains = -1.994*(10**-10)*(temps**3) + 5.807*(10**-7)*(temps**2) + 5.992*(10**-4)*temps -0.22

    ### Okay, we've done a lot of different types of analyses, let's put them together

    plt.figure()
    plt.plot(temps-273,strains,'b')
    plt.xlabel('Temperature (C)')
    plt.ylabel('Macroscopic Strain (%)')
    plt.legend(['-(1.994e-10)*(T^3) + (5.807e-7)*(T^2) + (5.992e-4)*T - 0.22'],loc=1)
    plt.savefig('Ti64_eps_T_polynomial.tif')

    plt.figure()
    plt.plot(np.true_divide(beta_index,hertz),beta_strain,'b.',np.true_divide(beta_index[13:103],hertz),beta_strain[13:103],'g.',np.true_divide(alpha_index,hertz),alpha_strain,'r.',np.true_divide(alpha_index[len(alpha_index)-100:len(alpha_index)],hertz),alpha_strain[len(alpha_strain)-100:len(alpha_strain)],'g.')
    plt.xlabel('Time (s)')
    plt.ylabel('Strain')
    #plt.ylim([0,0.04])
    plt.legend(['alpha','alpha fit vals','beta','beta fit vals'])
    plt.title('Measured Lattice Strains (a)')
    plt.savefig(str(title) + '_strains.png')

    plt.figure()
    plt.plot(np.true_divide(beta_index,hertz),beta_temp,'b.',np.true_divide(alpha_index,hertz),alpha_temp,'r.')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (C)')
    plt.legend(['Beta','Alpha'])
    plt.title('Temperature Measurements from Coefficients of Thermal Expansion')
    plt.savefig(str(title) + '_coef_temp.png')

    plt.figure()
    plt.plot(np.true_divide(macroindex,hertz),abs(roots[:,2])-273)
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (C)')
    plt.title('Polynomial Fit of Temperature from Separate Phase Strains')
    plt.savefig(str(title) + '_poly_temp.png')

    plt.figure()
    plt.plot(np.true_divide(beta_index,hertz),abs(beta_roots[:,2])-273,'b.',np.true_divide(alpha_index,hertz),abs(alpha_roots[:,2])-273,'r.')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (C)')
    plt.legend(['Beta','Alpha'])
    plt.savefig(str(title) + '_indiv_poly_fit.png')


if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3])
