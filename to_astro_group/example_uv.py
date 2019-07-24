#!/usr/bin/env python
"""

"""

from __future__ import print_function
import numpy as np
import emcee
import matplotlib.pyplot as plt
import corner

### First, define the probability distribution 
### that you would like to sample.


### Posterior for measuring x alone
def lnposteriorx1st(uvec, x, sigmax):
    lnpr=lnprior1(uvec)
    if np.isinf(-lnpr):
        return lnpr
    else:
        return lnpr+lnlikelihoodx(uvec, x, sigmax)

def lnlikelihoodx(uvec, x, sigmax):
    xmod=(uvec[0]**2.+uvec[1]**2.)**0.5
    
    return -0.5*((xmod-x)/sigmax)**2.



### Posterior for measuring y alone
def lnposteriory1st(uvec, y, sigmay):
    lnpr=lnprior1(uvec)
    if np.isinf(-lnpr):
        return lnpr
    else:
        return lnpr+lnlikelihoody(uvec, y, sigmay)

def lnlikelihoody(uvec, y, sigmay):
    ymod=2./np.pi*np.arctan(uvec[1]/uvec[0])
    
    return -0.5*((ymod-y)/sigmay)**2.



### Posterior for measuring x and y and assuming statistical independence.
def lnposteriorxy1st(uvec, x, sigmax, y, sigmay):
    lnpr=lnprior1(uvec)
    if np.isinf(-lnpr):
        return lnpr
    else:
        return lnpr+lnlikelihoodxy(uvec, x, sigmax, y, sigmay)

def lnlikelihoodxy(uvec, x, sigmax, y, sigmay):
    xmod=(uvec[0]**2.+uvec[1]**2.)**0.5
    ymod=2./np.pi*np.arctan(uvec[1]/uvec[0])
    
    return -0.5*((xmod-x)/sigmax)**2.-0.5*((ymod-y)/sigmay)**2.




### The two priors I invented
def lnprior1(uvec): 
    if (0. <= uvec[0] <= 1./(2.**0.5)) and \
            (0. <= uvec[1] <= 1./(2.**0.5)):
        return 0.
    else:
        return -np.inf
        
def lnprior2(uvec):
    if (0. <= uvec[0] <= 1./(2.**0.5)) and \
            (0. <= uvec[1] <= 1./(2.**0.5)):
        return np.log( 1./((uvec[0]-1./(2.**0.5))**2.+(uvec[1]-1./(2.**0.5))**2.) )
    else:
        return -np.inf    

#########################




nwalkers = 2500 ### We'll sample with 2500 walkers.
nchain=1000
nburnin=100
ndim = 2


### Prior knowledge: u and v can only be in this range:
uvlims= [
        [0.,1./(2.**0.5)],
        [0.,1./(2.**0.5)]
        ]

### These are the assumed observations:
### (Choose 0 <= x <= 1 and 0 <= y <= 1.) 
obs_x=0.5 
obs_sigmax=0.05
obs_y=0.8
obs_sigmay=0.05

### 
constrain_x=True
constrain_y=False
constrain_xy=False







##########################

### Choose an initial set of positions for the walkers.
### (It is a list of ndim-dimensional lists.)
p0 = [  [np.random.uniform(uvlims[0][0],uvlims[0][1]),
        np.random.uniform(uvlims[1][0],uvlims[1][1])] 
            for i in xrange(nwalkers)]



### Initialize the sampler with the chosen specs.
samplerx = emcee.EnsembleSampler(nwalkers, ndim, lnposteriorx1st, \
    args=[obs_x, obs_sigmax], a=2, threads=1)
samplery = emcee.EnsembleSampler(nwalkers, ndim, lnposteriory1st, \
    args=[obs_y, obs_sigmay], a=2, threads=1)
samplerxy = emcee.EnsembleSampler(nwalkers, ndim, lnposteriorxy1st, \
    args=[obs_x, obs_sigmax, obs_y, obs_sigmay], a=2, threads=1)


if constrain_x == True:
    print("SAMPLING for observable x...")
    pos, prob, state = samplerx.run_mcmc(p0, nchain)
    print("SAMPLING DONE.")

if constrain_y == True:
    print("SAMPLING for observable y...")
    pos, prob, state = samplery.run_mcmc(p0, nchain)
    print("SAMPLING DONE.")

if constrain_xy == True:
    print("SAMPLING for observables x and y...")
    pos, prob, state = samplerxy.run_mcmc(p0, nchain)
    print("SAMPLING DONE.")



### Uncheck this to see the evolution of the probabilities 
### for 'measuring x alone'
if 1==2 and constrain_x == True:
    lnprobabilityx=[]
    lnprobabilityx.append(samplerx.lnprobability)

    poslnprobx=np.arange(1,len(lnprobabilityx[0][0])+1)
    fig=plt.figure(figsize=(6,6))
    ax=plt.subplot(1,1,1)
        
    for i in xrange(0,len(lnprobabilityx[0])):
        plt.plot(poslnprobx,[np.arcsinh(lnprobabilityx[0][i][j]) \
            for j in range(0,len(lnprobabilityx[0][0]))], \
            color='black', linewidth=0.05, linestyle='-')
    plt.ylabel("$\\arcsinh(\\ln(\\mathrm{probx}))$")
    plt.xlabel("position in the chain")
    plt.show()






### Here, I do the corner plots

rangetriang=[(uvlims[0][0]-0.1*(uvlims[0][1]-uvlims[0][0]),\
                uvlims[0][1]+0.1*(uvlims[0][1]-uvlims[0][0])),
                (uvlims[1][0]-0.1*(uvlims[1][1]-uvlims[1][0]),\
                uvlims[1][1]+0.1*(uvlims[1][1]-uvlims[1][0]))]

if constrain_x == True:
    samplesx = samplerx.chain[:, nburnin:, :].reshape((-1, ndim))
    fig = corner.corner(samplesx, labels=["$u$", "$v$"], color="red", \
        range=rangetriang, bins=60)
    fig.savefig("uv_triangx.png")

if constrain_y == True:
    samplesy = samplery.chain[:, nburnin:, :].reshape((-1, ndim))
    fig = corner.corner(samplesy, labels=["$u$", "$v$"], color="blue", \
        range=rangetriang, bins=60)
    fig.savefig("uv_triangy.png")

if constrain_xy == True:
    samplesxy = samplerxy.chain[:, nburnin:, :].reshape((-1, ndim))
    fig = corner.corner(samplesxy, labels=["$u$", "$v$"], color="purple", \
        range=rangetriang, bins=60)
    fig.savefig("uv_triangxy.png")



