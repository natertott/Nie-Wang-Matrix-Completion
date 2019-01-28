import sys
import numpy as np
import pf_regex

def write_text(lst_path,write_path):
    mat = pf_regex.pfvt(lst_path)
    fracs = np.column_stack((mat[0],mat[1],mat[2]))
    alpha = fracs[:,0]
    beta = fracs[:,1]
    liq = fracs[:,2]

    det = input('Detector Type')
    pos = input('Position')
    speed = input('Speed')

    f = open(write_path+".txt","w+")
    f.write("Detector Type: " + det + "\n")
    f.write("Beam Position: " + pos + "\n")
    f.write("Wire feed speed: " + speed + "\n")
    f.write("Alpha phase fraction (relative) \n")
    np.savetxt(f,alpha)
    f.write("Beta phase fraction (relative) \n")
    np.savetxt(f,beta)
    f.write("Liquid phase fraction (relative) \n")
    np.savetxt(f,liq)

    f.close()


if __name__ == "__main__":
    write_text(sys.argv[1],sys.argv[2])
