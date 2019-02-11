####################################
# Nathan S Johnson
# DIFFraction matrix COMPletion (diffcomp)
# 11/7/2018
###################################

import sys
import numpy as np
from matplotlib import pyplot as plt
import NieWangMatComp as matcomp
import scipy as sc
import math
from itertools import product

def read_ge(filename):
    """
    Reads the binary output from GE detectors.

    :param filename: The name of the file that contains the dataself.
    :type filename: str
    :returns: Frames contained in the GE file.
    :rtype: numpy.ndarray
    """
    with open(filename) as ifs:
        ifs.seek(8192)
        data = np.fromfile(ifs, np.uint16)
    return np.reshape(data, (-1, 2048, 2048))

############
def subtract(diff,dark):
    diffall = read_ge(diff)

    # The read_ge function outputs all the stacked images into stacked arrays
    # We only really want the first array
    diffimg = diffall[0,:,:]
    np.shape(diffimg)

    darkall = read_ge(dark)
    darkimg = darkall[0,:,:]
    thresh = diffimg - darkimg
    plt.imsave('Original Image',thresh)
    return thresh
    #plt.imshow(thresh)
    #plt.show()

    #plt.plot(thresh)
    #plt.show()


#############
## This is a newton-raphson minimization function that will minimize the entires of M in the main algorithm
# We will be taking observed values from the matrix M and returning a matrix Emin, which will be used to update the
# final update matrix X

def newraps0(M,lam,p): #takes a matrix, the lam parameter from the main algorithm, and the power p as args
    x,y = M.shape
    Emin=np.zeros_like(M) # minimized matrix to be output
    error = 0
    nu = ((lam*p*(1-p))//(2))**(1//(2-p)) #this is our initial 'guess' at what the minimized value will be
    for i,j in product(np.arange(x),np.arange(y)):
        b = 2*nu - 2 * M[i,j] + lam*p*nu**(p-1)
        if b < 0:
            sigma = 2*nu
            iter = range(10)
            for cc in iter:
                f = 2*sigma - 2*M[i,j] + lam*p*sigma**(p-1)
                #quadratic function to be minimized, based on the decoupled version of the overall problem
                df = 2 + lam*p*(p-1)*sigma**(p-2) #derivative of f
                sigma = sigma - f/df #update sigma
            error = sigma**2 - 2*sigma*M[i,j] + lam*sigma**p
            #error estimation of our 'guess' at the quadratic's minimum, should be a negative value
            if error > 0:
                sigma = 0
        else:
            sigma = 0
        Emin[i,j] = sigma
    return Emin

################
# This is a newton-raphson function that takes the singular values of the M matrix (from the main algorithm) and
# uses them to estimate the singular values of the main algorithms output, X

# Note that the singular values are the diagonal entries of S in USV^T = svd(M), the reduced SVD of M

def newraps(s,lam,p): #s is a 1D array of the singular values of M, lam is a parameter from the main algorithm, p is the power
    smin = np.zeros_like(s) #output array
    store = np.zeros(3)
    store[0] = 0
    error = np.zeros(3)
    error[0] = 0
    count = 0
    for alpha in s:
        #first we estimate the minimum of the quadratic from the RHS
        nu = (lam*p*(1-p)/2)**(1/(2-p)) #this is our initial 'guess'
        b = 2*nu - 2*alpha + lam*p*alpha**(p-1)
        if b<0:
            sigma = 2*nu

            it = range(10)
            for i in it:
                f = 2*sigma -2*alpha + lam*p*sigma**(p-1) #quadratic to be minimized over
                df = 2 + lam*p*(p-1)*sigma**(p-2) #derivative of f
                sigma = sigma - f/df #update sigma
            store[1] = sigma
            error[1] = sigma**2 - 2*sigma*alpha + lam*abs(sigma)**p #error estimation of our guess
        else:
            store[1]=1
            error[1]=float('inf')

        #next we estimate the minimum of the quadratic function from the LHS
        nu = - nu #our initial guess
        b = 2*nu - 2*alpha + lam*p*abs(alpha)**(p-1)
        if b>0:
            sigma = 2*nu

            it = range(10)
            for j in it:
                f = 2*sigma -2*alpha - lam*p*abs(sigma)**(p-1)
                df = 2 + lam*p*(p-1)*abs(sigma)**(p-2)
                sigma = sigma - f/df
            store[2] = sigma
            error[2] = sigma**2 - 2*sigma*alpha + lam*abs(sigma)**p
        else:
            store[2]=1
            error[2]=float('inf')
        ide = np.where(error==error.min())

        #the following code is a work-around I implemented because I was having trouble dealing with different kinds of data types
        #essentially, whenever the minimization function cannot find a minimum, the np.where function above returns a tuple whose first entry is an array with the values [0 1 2]
        #this means that every value in store had the same error (inf)
        #however, you cannot assign an array as an index in store, so I built this workaround to detect when the minimization function fails, and automatically set that smin value to be 0

        try:
            smin[count] = store[ide]
        except ValueError:
            smin[count] = store[0]
        except IndexError:
            smin[count] = store[0]
        count += 1
    return smin
pass



######################################
#the following is the main body of the matrix completion algorithm

## D = original matrix to be completed
## D0 = matrix with only observed values from D
## X = matrix which will serve as the completed version of D
## Xo = you get the idea

def matcomp(D,r,p):
    [a,b] = D.shape #get the shape of the data matrix being analyzed
    x = range(a)
    y = range(b)
    rho = 1.2 # parameter for augmented lagrangian method
    mu = 0.1 # parameter for augmented lagrangian method
    inmu = 1/mu #obvi

    X = np.zeros((a,b))
    Lo = np.zeros((a,b)) # Lo = Eo - Xo + Do
    Sigma = np.zeros((a,b)) # Sigma = X - Z
    Eo = np.zeros((a,b)) # Minimized version of Xo - Do - (1/mu)Lo
    Z = np.zeros((a,b)) # Minimized SVD of csv

    ## Define a matrix A that satisfies the following condition: if csv(i,j) is observed, then A(i,j) is 1; else, A(i,j)
    ## is 0
    A = np.zeros((a,b))
    for i in x:
        for j in y:
            if D[i,j] != 0:
                A[i,j] = 1

    Xo = D*A # This is where we store the updated versions of observed entires from D
    Do = D*A
    plt.imsave('Do.png',Do)
    Nope = np.zeros((a,b)) + 1
    Nope = Nope - A
    Xn = D*Nope # This is where we store the currently-unobserved entries that we will soon be finding

    counter = np.arange(63)
   
    for count in counter: #this is our convergence condition
	print(count)
        #Hold X, Eo constant
        #Optimize on Z
        M = X + inmu*Sigma # temporary matrix we are using to find the SVD of the expression shown
        U, S, V = np.linalg.svd(M,full_matrices = False) #the reduced SVD of M
        lam = 2*r*inmu #parameter for Newton-Raphson method
        Stemp = newraps(S,lam,p) #Find the optimized values for the singular values of M
        Z = np.dot(np.dot(U,np.diag(Stemp)),V) #Update Z

        #Now we hold X, Z constant
        # Optimize on Eo
        M = Xo - Do - inmu*Lo       # temporary matrix to be operated upon
        lam = 2*r*inmu      #same as above
        Eo = newraps0(M,lam,p)       #Updated Eo

        ## Finally, we hold Eo, Z constant
        ##Optimize on X
        N = Z - inmu*Sigma        #Temporary
        No = N*A
        Mo = Eo + Do + inmu*Lo      #Temporary
        Xn = N*Nope     #storing the new values of the previously unobserved entires
        X0 = (Mo + No)/2
        X = Xo + Xn #Updating our output matrix

        Sigma = Sigma + mu*(X - Z)
        Lo = Lo + mu*(Eo - Xo + Do)
        mu = rho*mu
    return X

#######################
def complete(diff,dark,fname):
    #ge_path is the path to the .ge file
    #first we need to read the file in
    orig = read_ge(diff)
    original = orig[0,:,:]

    #background subtraction
    sub = subtract(diff,dark)
    
    # We need to subsample the original image to produce a spare matrix
    # First generate a random set of coordinates to either remove or keep
    # I'm gonna try first to just take a few rows out of the image
    
    these =[np.random.randint(0,2047,400000), np.random.randint(0,2047,400000)] 
    input = np.zeros([2048,2048])
    input[these] = sub[these]
    input_norml = np.divide(input,np.amax(input))
    plt.imsave(str(fname)+'Input Image.png',input_norml,cmap='Spectral')

   

    # Alright let's try to do some matrix completion
    r = 1.1 #this is the first input parameter for matrix completion, should always be greater than 1
    rho = 1.5 #the other parameter for matrix completion

    comp = matcomp(input,r,rho)

    for i in np.arange(0,2047):
	for j in np.arange(0,2047):
	    if abs(comp[i,j]) > 50000:
		comp[i,j] = 0
    comp[2047,2047] = 0

    plt.imsave(str(fname) +'Completed Matrix',comp,cmap='Spectral')
    difference = comp - original 
    plt.imsave(str(fname)+'Diffrence of Images',difference,cmap='Spectral')

if __name__ == '__main__':
    complete(sys.argv[1], sys.argv[2], sys.argv[3])
