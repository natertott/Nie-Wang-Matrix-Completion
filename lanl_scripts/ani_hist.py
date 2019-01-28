import sys
import numpy as np
import math
import pandas
import glob
from matplotlib import pyplot as plt
from matplotlib import animation as animation

# This is a code to animate temperature data obtained from diffraction analysis

hertz = input('Detector Hertz?')
hertz = float(hertz)
wavey = input('Wavelength (nm)?')
fig, ax = plt.subplots()
ln, = ax.plot([],[],animated=True)

x = []
intensity =[]

colnames= ['x','intensity','garbage','garbage','garbage']
for file_path in sorted(glob.glob('*.csv')):
    with open(file_path,'r') as powder:
        data = pandas.read_csv(powder,names=colnames)
        store_x = data.x
        store_int = data.intensity
        x.append(store_x[13:len(store_x)])
        intensity.append(store_int[13:len(store_int)])


x = np.asarray(x,dtype='float')
intensity = np.asarray(intensity,dtype='float')

d = np.zeros(len(x[0,:]))
for ii in np.arange(len(x[0,:])):
    d[ii] = wavey/(2*math.sin(x[0,ii]))

ln, = ax.plot(x[0,:],intensity[0,:])
ax.set_ylim(0,max(intensity[len(intensity)-1,:]))

print d

def update_histogram(frame):
    ln.set_xdata(x[frame,:])
    ln.set_ydata(intensity[frame,:])
    return ln,

plt.xlabel('2 theta')
plt.ylabel('Intensity')
anim = animation.FuncAnimation(fig,update_histogram,frames=90,interval=100, blit = False)
anim.save('solidification.mp4',writer='ffmpeg')

plt.show()
