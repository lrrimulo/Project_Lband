#!/usr/bin/env python
"""

"""

from __future__ import print_function
import numpy as np
import emcee
import matplotlib.pyplot as plt
import corner
import pyhdust.lrr as lrr

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



### Posterior for measuring x and y and assuming statistical 
### independence.
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
        return np.log( 1./((uvec[0]-1./(2.**0.5))**2.\
                    +(uvec[1]-1./(2.**0.5))**2.) )
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
obs_y=0.2
obs_sigmay=0.05

### 
constrain_x=False
constrain_y=True
constrain_xy=False

###
constrain_int_x=False
constrain_int_y=False
constrain_int_xy=False




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
    fig = corner.corner(samplesxy, labels=["$u$", "$v$"], 
        color="purple", range=rangetriang, bins=60)
    fig.savefig("uv_triangxy.png")



















###################################################


def x_uv(u,v):
    
    return (u**2.+v**2.)**0.5
    
def y_uv(u,v):
    
    return 2./np.pi*np.arctan(v/u)







### Posterior for measuring x alone
def lnposteriorx2nd(uvec, x, sigmax, axis, x_values):
    lnpr=lnprior1(uvec)
    if np.isinf(-lnpr):
        return lnpr
    else:
        return lnpr+lnlikelihoodx_int(uvec, x, sigmax, axis, x_values)

def lnlikelihoodx_int(uvec, x, sigmax, axis, x_values):
    xmod=lrr.interpLinND([uvec[0],uvec[1]],axis,x_values,allow_extrapolation="no")
    if np.isnan(xmod):
        return -np.inf
    else:
        return -0.5*((xmod-x)/sigmax)**2.



### Posterior for measuring y alone
def lnposteriory2nd(uvec, y, sigmay, axis, y_values):
    lnpr=lnprior1(uvec)
    if np.isinf(-lnpr):
        return lnpr
    else:
        return lnpr+lnlikelihoody_int(uvec, y, sigmay, axis, y_values)

def lnlikelihoody_int(uvec, y, sigmay, axis, y_values):
    ymod=lrr.interpLinND([uvec[0],uvec[1]],axis,y_values,allow_extrapolation="no")

    if np.isnan(ymod):
        return -np.inf
    else:
        return -0.5*((ymod-y)/sigmay)**2.



### Posterior for measuring x and y and assuming statistical 
### independence.
def lnposteriorxy2nd(uvec, x, sigmax, y, sigmay, axis, x_values, \
                            y_values):
    lnpr=lnprior1(uvec)
    if np.isinf(-lnpr):
        return lnpr
    else:
        return lnpr+lnlikelihoodxy_int(uvec, x, sigmax, y, sigmay, \
                            axis, x_values, y_values)

def lnlikelihoodxy_int(uvec, x, sigmax, y, sigmay, axis, x_values, \
                            y_values):
    xmod=lrr.interpLinND([uvec[0],uvec[1]],axis,x_values,allow_extrapolation="no")
    ymod=lrr.interpLinND([uvec[0],uvec[1]],axis,y_values,allow_extrapolation="no")

    if np.isnan(xmod) or np.isnan(ymod):
        return -np.inf
    else:
        return -0.5*((xmod-x)/sigmax)**2.-0.5*((ymod-y)/sigmay)**2.










umin=0.01
umax=1./(2.**0.5)
Nu=12
u_axis=[umin+(umax-umin)*float(i)/float(Nu-1) for i in range(0,Nu)]

vmin=0.0
vmax=1./(2.**0.5)
Nv=9
v_axis=[vmin+(vmax-vmin)*float(i)/float(Nv-1) for i in range(0,Nv)]


def is_in_region(u,v,typ=None):

    if typ == None:
        return "no"
    
    if typ == "circular":
        ucent=0.4
        vcent=0.4
        radius=0.2
        rn=np.sqrt((u-ucent)**2.+(v-vcent)**2.)
        if rn <= radius:
            return "yes"
        else:
            return "no"
    
    if typ == "above_line":
        v0=0.2
        A=1.0
        vline=v0+A*u
        if v > vline:
            return "yes"
        else:
            return "no"    
    
    



### Atributing values to every point in the grid.
x_values=[]
y_values=[]
for i in xrange(0,len(u_axis)):
    for j in xrange(0,len(v_axis)):
        isinregion=is_in_region(u_axis[i],v_axis[j],typ=None)
        if isinregion == "no":
            auxi_x=x_uv(u_axis[i],v_axis[j])
            auxi_y=y_uv(u_axis[i],v_axis[j])
        else:
            auxi_x=np.nan
            auxi_y=np.nan            
        x_values.append(auxi_x)
        y_values.append(auxi_y)

### The list 'axis', which enters the routine:
axis=[u_axis,v_axis]

### Several examples to be printed on screen:
if 1==2:
    print("")
    point=np.array([0.2,0.3]) 
    print("EXAMPLE: Point in 2D space inside the grid.")
    print("Point: ",point)
    interpx=lrr.interpLinND(point,axis,x_values)
    print("The interpolated value for x: ",interpx)
    f_example=x_uv(point[0],point[1])
    print("The real value for x: ",f_example)
    interpy=lrr.interpLinND(point,axis,y_values)
    print("The interpolated value for y: ",interpy)
    f_example=y_uv(point[0],point[1])
    print("The real value for y: ",f_example)
    import sys; sys.exit()


### Choose an initial set of positions for the walkers.
### (It is a list of ndim-dimensional lists.)
p0 = [  [np.random.uniform(uvlims[0][0],uvlims[0][1]),
        np.random.uniform(uvlims[1][0],uvlims[1][1])] 
            for i in xrange(nwalkers)]



### Initialize the sampler with the chosen specs.
samplerx_int = emcee.EnsembleSampler(nwalkers, ndim, lnposteriorx2nd, \
    args=[obs_x, obs_sigmax, axis, x_values], a=2, threads=1)
samplery_int = emcee.EnsembleSampler(nwalkers, ndim, lnposteriory2nd, \
    args=[obs_y, obs_sigmay, axis, y_values], a=2, threads=1)
samplerxy_int = emcee.EnsembleSampler(nwalkers, ndim, lnposteriorxy2nd, \
    args=[obs_x, obs_sigmax, obs_y, obs_sigmay, axis, x_values, \
            y_values], a=2, threads=1)


if constrain_int_x == True:
    print("SAMPLING for observable x...")
    pos, prob, state = samplerx_int.run_mcmc(p0, nchain)
    print("SAMPLING DONE.")

if constrain_int_y == True:
    print("SAMPLING for observable y...")
    pos, prob, state = samplery_int.run_mcmc(p0, nchain)
    print("SAMPLING DONE.")

if constrain_int_xy == True:
    print("SAMPLING for observables x and y...")
    pos, prob, state = samplerxy_int.run_mcmc(p0, nchain)
    print("SAMPLING DONE.")









### Here, I do the corner plots

rangetriang=[(uvlims[0][0]-0.1*(uvlims[0][1]-uvlims[0][0]),\
                uvlims[0][1]+0.1*(uvlims[0][1]-uvlims[0][0])),
                (uvlims[1][0]-0.1*(uvlims[1][1]-uvlims[1][0]),\
                uvlims[1][1]+0.1*(uvlims[1][1]-uvlims[1][0]))]

if constrain_int_x == True:
    samplesx_int = samplerx_int.chain[:, nburnin:, :].reshape((-1, ndim))
    fig = corner.corner(samplesx_int, labels=["$u$", "$v$"], color="red", \
        range=rangetriang, bins=60)
    fig.savefig("uv_triangx_int.png")

if constrain_int_y == True:
    samplesy_int = samplery_int.chain[:, nburnin:, :].reshape((-1, ndim))
    fig = corner.corner(samplesy_int, labels=["$u$", "$v$"], color="blue", \
        range=rangetriang, bins=60)
    fig.savefig("uv_triangy_int.png")

if constrain_int_xy == True:
    samplesxy_int = samplerxy_int.chain[:, nburnin:, :].reshape((-1, ndim))
    fig = corner.corner(samplesxy_int, labels=["$u$", "$v$"], 
        color="purple", range=rangetriang, bins=60)
    fig.savefig("uv_triangxy_int.png")


