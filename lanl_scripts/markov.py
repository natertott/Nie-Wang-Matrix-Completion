import sys
import numpy as np
from matplotlib import pyplot as plt
import pf_regex

def plot_svd(lst_path):
    mat = pf_regex.pfvt(lst_path)
    #zeroth returned argument is alpha_frac, first is beta_frac, second is liq_frac
    fracs = np.column_stack((mat[0],mat[1],mat[2]))
    index = np.arange(0,len(fracs))
    smat = np.zeros([len(fracs),3])

    for i in index:

        store = np.array([fracs[i,0],fracs[i,1],fracs[i,2]])

        mat =np.zeros([9,6])
        mat[0,0] = fracs[i,0]
        mat[1,0] = fracs[i,1]
        mat[2,0] = fracs[i,2]
        mat[3,1] = fracs[i,0]
        mat[4,1] = fracs[i,1]
        mat[5,1] = fracs[i,2]
        mat[6,2] = fracs[i,0]
        mat[7,2] = fracs[i,1]
        mat[8,2] = fracs[i,2]
        mat[0,3] = 1
        mat[3,3] = 1
        mat[6,3] = 1
        mat[1,4] = 1
        mat[4,4] = 1
        mat[7,4] = 1
        mat[2,5] = 1
        mat[5,5] = 1
        mat[8,5] = 1

        store = np.linalg.lstsq(mat.T, np.array([store[0],store[1],store[2],1,1,1]))
        solut = np.asarray(store[0])
        P = np.array([[solut[0],solut[1],solut[2]],[solut[3],solut[4],solut[5]],[solut[6],solut[7],solut[8]]])
        u,s,v, = np.linalg.svd(P)
        smat[i,0] = s[0]
        smat[i,1] = s[1]
        smat[i,2] = s[2]


    time = index/200.
    for ii in index:
        plt.plot(smat[ii,:])

    plt.xlabel('Time')
    plt.ylabel('SV Magnitude')
    plt.title(sys.argv[1])
    plt.show()



if __name__ == '__main__':
    plot_svd(sys.argv[1])
