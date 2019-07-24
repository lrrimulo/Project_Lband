#!/usr/bin/env python
"""


"""

from __future__ import print_function
import numpy as np
import emcee
import matplotlib.pyplot as plt

### First, define the probability distribution 
### that you would like to sample.

def lnprob(x, mu, sigg):
    ### One dimensional gaussian distribution.
    
    diff = x-mu
    return -0.5*(diff/sigg)**2.

def lnprobsine(x):
    ### Returns A+sin^2(3*pi*x)
    
    A=0.00 ### Try A<0 to see emcee giving an error message!
    if 0. <= x <= 1.:
        return np.log(A+np.sin(x*3.*np.pi)**2)
    else:
        return -np.inf


def lnMBdistribution(x,aa):
    ### Returns the Maxwell-Boltzmann's distribution
    if x >= 0.:
        return np.log((2./np.pi)**0.5/aa**3.*x**2.*\
            np.exp(-0.5*x**2./aa**2.))
    else:
        return -np.inf


####################################


nwalkers = 2500 ### We'll sample with 2500 walkers.
Nchain = 1000 ### A chain with a 1000 elements looks good enough!
ndim = 1    ### Number of dimensions of the distribution. 
            ### (In this simple example, it is one.)

### Below, uncomment one distribution with its parameters:
disttype = "gaussian"; means = 0.5; sigg = 0.1
#disttype = "sine"
#disttype = "MB"; aa=0.3

### Choose an initial set of positions for the walkers.
### (It is a list of ndim-dimensional arrays.)
p0 = [0.1+0.001*np.random.rand(ndim) for i in xrange(nwalkers)]




### Initialize the sampler with the chosen specs.
if disttype == "gaussian":
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, \
        args=[means, sigg], a=2, threads=1)
if disttype == "sine":
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprobsine, \
        a=2, threads=1)
if disttype == "MB":
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnMBdistribution, \
        args=[aa], a=2, threads=1)


### Running MCMC:
print("SAMPLING...")
pos, prob, state = sampler.run_mcmc(p0, Nchain)
print("SAMPLING DONE.")
print("")


if 1==2:
    ### 'sampler.chain' contains the Markov chains generated. 
    print("'sampler.chain' is an array with the following shape: ",\
        sampler.chain.shape)
    print(len(sampler.chain)," is the number of walkers/chains")
    print(len(sampler.chain[0])," is the size of the chains")
    print(len(sampler.chain[0][0])," is the number of dimensions.")
    print("")
    print("'sampler.lnprobability' is an array with the foll. shape: ",\
        sampler.lnprobability.shape)
    print(len(sampler.lnprobability)," is the number of walkers/chains")
    print(len(sampler.lnprobability[0])," is the size of the chains")    
    


hlims=[-0.2,1.2] ### These are the limits for the histograms below.
if 1==2:
    ### This is a collection of histograms following the spread of the 
    ### walkers. Each plot contains the position of all walkers at a 
    ### specific position in the chain. It is interesting to see that
    ### the walkers will converge to the desired distribution, 
    ### independent on their initial positions! 
    for i in range(Nchain):
        plt.figure()
        plt.hist(np.array([sampler.chain[j][i][0] \
            for j in xrange(0,len(sampler.chain))]), \
            100, color="k", histtype="step", range=hlims)
        plt.title("Chain position: {0:d}".format(i+1))
        plt.show()


if 1==1:
    plt.figure()
    arraytohist=[]
    for i in xrange(200,Nchain):    ### Don't start the counting from 0!
                                    ### Exclude the burnin phase!
        for j in xrange(0,len(sampler.chain)):
            arraytohist.append(sampler.chain[j][i][0])
    arraytohist=np.array(arraytohist)
    plt.hist(arraytohist, \
        100, color="k", histtype="step", range=hlims, density=True)
    plt.title("Distribution")

    ### I decided to study a little more the case of the 
    ### MB distribution:
    if disttype == "MB":
        xd=np.array([hlims[0]+float(i)*(hlims[1]-hlims[0])/300. \
                    for i in xrange(0,301)])
        yd=np.array([np.exp(lnMBdistribution(xd[i],aa)) \
            for i in xrange(0,len(xd))])
        plt.plot(xd,yd,color="blue")
        
        print("RESULTS FOR MAXWELL-BOLTZMANN'S DISTRIBUTION:")
        print("a = ",aa)
        print("")
        
        vmp=sampler.flatchain[np.argmax(sampler.flatlnprobability)][0]
        print("most prob speed (correct) = ",np.sqrt(2.*aa**2.))
        print("most prob speed (obtained) = ",vmp)
        print("")        
        vmean=np.mean(arraytohist)
        print("mean speed (correct) = ",2.*np.sqrt(2./np.pi*aa**2.))
        print("mean speed (obtained) = ",vmean)
        print("")
        vstd=np.std(arraytohist)
        vrms=np.sqrt(vstd**2.+vmean**2.)
        print("rms speed (correct) = ",np.sqrt(3.*aa**2.))
        print("rms speed (obtained) = ",vrms)
        print("")
        

    plt.show()



