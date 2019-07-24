#!/usr/bin/env python
"""

"""

from __future__ import print_function
import numpy as np
import emcee
import matplotlib.pyplot as plt
import random
import corner





### Probability distributions for PART 1:
def ln_likelihood1c(M, mu_c,sigma_c,rho_c,r_c,sigmar_c,rhor_c,\
        xlim_c,ylim_c,masslim):
    mu=np.array([M[0],M[1]])
    r=np.array([M[2],M[3]])
    mass=M[4]
    
    ln_mu=ln_likelihood1_mu(mu, mu_c, sigma_c, rho_c)
    ln_posc=ln_likelihood1_posc(r, mass, r_c, sigmar_c, rhor_c, \
        xlim_c, ylim_c, masslim)
    ln_mass=ln_likelihood1_massc(mass, masslim)
    
    return ln_mu+ln_posc+ln_mass


def ln_likelihood1f(M, mu_f,sigma_f,rho_f,xlim_c,ylim_c,masslim):
    mu=np.array([M[0],M[1]])
    r=np.array([M[2],M[3]])
    mass=M[4]
    
    ln_mu=ln_likelihood1_mu(mu, mu_f, sigma_f, rho_f)
    ln_posf=ln_likelihood1_posf(r, xlim_c, ylim_c)
    ln_mass=ln_likelihood1_massf(mass, masslim)
    
    return ln_mu+ln_posf+ln_mass


def ln_likelihood1_mu(mu, mupar, sigmapar, rhopar):
    diff = mu-mupar

    return -(1./2./(1.-rhopar**2))*(\
        (diff[0]/sigmapar[0])**2+(diff[1]/sigmapar[1])**2-\
        2.*rhopar*(diff[0]/sigmapar[0])*(diff[1]/sigmapar[1]))

def ln_likelihood1_posc(r, mass, rpar, sigmarpar, rhorpar, \
        xlim_c, ylim_c, masslim):
    ### The spacial distribution of the stars from the cluster is 
    ### also assumed to be a bivariate normal distribution.

    if (not (xlim_c[0] <= r[0] <= xlim_c[1])) or \
            (not (ylim_c[0] <= r[1] <= ylim_c[1])) or \
            (not (masslim[0] <= mass <= masslim[1])):
        return -np.inf
    
    diff = r-rpar

    return -(1./2./(1.-rhorpar**2))*(\
        (diff[0]/sigmarpar[0])**2+(diff[1]/sigmarpar[1])**2-\
        2.*rhorpar*(diff[0]/sigmarpar[0])*(diff[1]/sigmarpar[1]))

def ln_likelihood1_posf(r, xlim_c, ylim_c):
    ### Any position in the field is equally probable, 
    ### for a star from the field.
    
    if (not (xlim_c[0] <= r[0] <= xlim_c[1])) or \
            (not (ylim_c[0] <= r[1] <= ylim_c[1])):
        return -np.inf
    else:
        return 0.   
        

def ln_likelihood1_massc(mass, masslim):
    ### The mass distribution of the cluster follows a well known IMF.
    
    if not (masslim[0] <= mass <= masslim[1]):
        return -np.inf
    
    return -2.3*np.log(mass)

def ln_likelihood1_massf(mass, masslim):
    ### The mass distribution of the field is more steeper than the 
    ### IMF of the cluster
    
    if not (masslim[0] <= mass <= masslim[1]):
        return -np.inf
    
    return (-2.3-3.0)*np.log(mass)






##################
### Probability distributions for PART 2:


def ln_posterior_cluster(M, mulist, Mlim):

    lnprior=ln_prior(M, Mlim)
    if np.isinf(-lnprior):
        return -np.inf
    else:
    
        Mcluster=[M[0],M[1],M[2],M[3],M[4],M[5]]
        Mfield=[M[6],M[7],M[8],M[9],M[10],1.-M[5]]
    

    
        ln_post=0.
        for i in range(0,len(mulist)):
            ln_post+=np.log(    likelihood_cluster(Mcluster, \
                                    mulist[i][0], mulist[i][1])+\
                                likelihood_cluster(Mfield, \
                                    mulist[i][0], mulist[i][1])\
                                )
        ln_post+=lnprior
                    
                       
    return ln_post




def likelihood_cluster(Mreduced, muxpar, muypar):
    mux = Mreduced[0]
    muy = Mreduced[1]
    sigmax = Mreduced[2]
    sigmay = Mreduced[3]
    rho = Mreduced[4]
    n = Mreduced[5]
    #
    diffx = mux-muxpar
    diffy = muy-muypar

    return n/(2.*np.pi*sigmax*sigmay*(1.-rho**2)**0.5)*\
        np.exp(-(1./2./(1.-rho**2))*(\
        (diffx/sigmax)**2+(diffy/sigmay)**2-\
        2.*rho*(diffx/sigmax)*(diffy/sigmay)))

def ln_prior(M,Mlim):
    
    if not (Mlim[0][0] <= M[0] <= Mlim[0][1]) or \
            not (Mlim[1][0] <= M[1] <= Mlim[1][1]) or \
            not (Mlim[2][0] <= M[2] <= Mlim[2][1]) or \
            not (Mlim[3][0] <= M[3] <= Mlim[3][1]) or \
            not (Mlim[4][0] <= M[4] <= Mlim[4][1]) or \
            not (Mlim[5][0] <= M[5] <= Mlim[5][1]) or \
            not (Mlim[6][0] <= M[6] <= Mlim[6][1]) or \
            not (Mlim[7][0] <= M[7] <= Mlim[7][1]) or \
            not (Mlim[8][0] <= M[8] <= Mlim[8][1]) or \
            not (Mlim[9][0] <= M[9] <= Mlim[9][1]) or \
            not (Mlim[10][0] <= M[10] <= Mlim[10][1]):
        return -np.inf
    else:
        return 0.
    
    
    













### PART 1: GENERATING THE FAKES CLUSTER AND FIELD
part1 = False
nburnin=1000
nsampling=1


### N cluster
N_c=291
### mu cluster
mu_c=np.array([0.310,0.722])
sigma_c=np.array([0.077,0.077])
rho_c=0.

### N field
N_f=584
### mu field
mu_f=np.array([0.516,0.464])
sigma_f=np.array([0.155,0.155])
rho_f=-0.414

### pos cluster
r_c=np.array([0.000,0.000])
sigmar_c=np.array([0.3,0.3])
rhor_c=-0.4
xlim_c=np.array([-1.0,1.0])
ylim_c=np.array([-1.0,1.0])

### mass
masslim=np.array([0.5,100.0])






### PART 2: OBTAINING CLUSTER AND FIELD PARAMETERS FROM MEASUREMENTS
part2 = True
nwalkers2=650
nburnin2=100
nsampling2=100
outputtitle="example_steven_analysis.txt"


Mlim =  [
        [-1.0,1.0],         ### mux_cluster
        [-1.0,1.0],         ### muy_cluster
        [0.0,0.2],          ### sigmax_cluster
        [0.0,0.2],          ### sigmay_cluster
        [-0.0001,0.0001],   ### rho_cluster
        [0.1,1.0],          ### n_cluster

        [-1.0,1.0],         ### mux_field
        [-1.0,1.0],         ### muy_field
        [0.0,0.4],          ### sigmax_field
        [0.0,0.4],          ### sigmay_field
        [-1.0,1.0],         ### rho_field
        [0.1,1.0]           ### n_field
        ]


### PART 3: MEMBERSHIP PROBABILITY
part3 = False








if part1 == True:
    
    doingnow = ["cluster","field"]

    for ido in range(0,len(doingnow)):
    
        ### mux, muy
        ### rx, ry
        ### mass
        ndim = 5 
        if doingnow[ido] == "cluster":
            nwalkers = N_c*2
        else: 
            nwalkers = N_f*2
    
        ### Choose an initial set of positions for the walkers.
        ### (It is a list of ndim-dimensional lists.)
        p0 =    [
                [random.uniform(-1.,1.),                ### mux
                random.uniform(-1.,1.),                 ### muy
                random.uniform(xlim_c[0],xlim_c[1]),    ### rx
                random.uniform(ylim_c[0],ylim_c[1]),    ### ry
                random.uniform(0.5,100.)]               ### mass
                for i in xrange(nwalkers)]


        ### Initialize the sampler with the chosen specs.
        if doingnow[ido] == "cluster":
            argumentos=[mu_c,sigma_c,rho_c,r_c,sigmar_c,rhor_c,\
                xlim_c,ylim_c,masslim]
            sampler = emcee.EnsembleSampler(nwalkers, ndim, \
                ln_likelihood1c, args=argumentos, a=2, threads=1)
        else:
            argumentos=[mu_f,sigma_f,rho_f,xlim_c,ylim_c,masslim]
            sampler = emcee.EnsembleSampler(nwalkers, ndim, \
                ln_likelihood1f, args=argumentos, a=2, threads=1)
        
        lnprobability=[]
        ### Run 'nburnin' steps as a burn-in.
        print("PART 1: BURNIN for "+doingnow[ido]+"...")
        pos, prob, state = sampler.run_mcmc(p0, nburnin)
        lnprobability.append(sampler.lnprobability)
        ### Reset the chain to remove the burn-in samples.
        sampler.reset()
        ### Starting from the final position in the burn-in chain, 
        print("PART 1: SAMPLING for "+doingnow[ido]+"...")
        sampler.run_mcmc(pos, nsampling, rstate0=state)
        print("PART 1: SAMPLING for "+doingnow[ido]+" DONE.")
        lnprobability.append(sampler.lnprobability)
        

        if 1==2:
            ### Check, from the probabilities, if it has converged.
            poslnprob=np.arange(1,len(lnprobability[0][0])+1)
            fig=plt.figure(figsize=(6,6))
            ax=plt.subplot(1,1,1)
        
            for i in xrange(0,len(lnprobability[0])):
                plt.plot(poslnprob,[np.arcsinh(lnprobability[0][i][j]) \
                    for j in range(0,len(lnprobability[0][0]))], \
                    color='black', linewidth=0.05, linestyle='-')
            plt.show()
        

        ### I will use the last iteration of the Markov chain as the 
        ### values of proper motions, positions and masses of the fakes
        ### cluster and field:
        f0 = open("example_steven_"+doingnow[ido]+".txt", "w")
        for j in xrange(0,len(sampler.chain)/2):
            f0.write(doingnow[ido]+str(j+1)+" "+\
                str(sampler.chain[j][-1][0])+" "+\
                str(sampler.chain[j][-1][1])+" "+\
                str(sampler.chain[j][-1][2])+" "+\
                str(sampler.chain[j][-1][3])+" "+\
                str(sampler.chain[j][-1][4])+"\n")
        f0.close()

        sampler.reset()






### This second part could be used with real data of a cluster+field. 
### Very little modifications are needed.
if part2 == True:
    
    doingnow = ["cluster","field"]

    ### Reading files from cluster and field and shuffling them:
    files=["example_steven_cluster.txt","example_steven_field.txt"]

    data=[]
    for ifiles in range(0,len(files)):
        f1=open(files[ifiles],"r")
        linhas=f1.readlines()
        f1.close()

        for i in range(0,len(linhas)):
            data.append([linhas[i].split()[0],linhas[i].split()[1],\
                linhas[i].split()[2],linhas[i].split()[3],\
                linhas[i].split()[4],linhas[i].split()[5]])
            
    np.random.shuffle(data)
    
    
    ### Uncheck this if you want to see the proper motions and a map of 
    ### the cluster and field
    if 1==1:
        fig=plt.figure(figsize=(6,6))
        ax=plt.subplot(1,1,1)
        for i in xrange(0,len(data)):
            if "cluster" in data[i][0]:
                plt.scatter([float(data[i][1])],[float(data[i][2])], \
                        color='red', marker='+')
            else:
                plt.scatter([float(data[i][1])],[float(data[i][2])], \
                        color='black', marker='+')
        plt.xlabel("$\mu_x$")
        plt.ylabel("$\mu_y$")
        plt.savefig("propermotions.png")
        plt.close(fig)


        fig=plt.figure(figsize=(6,6))
        ax=plt.subplot(1,1,1)
        for i in xrange(0,len(data)):
            if "cluster" in data[i][0]:
                plt.scatter([float(data[i][3])],[float(data[i][4])], \
                    color='red', marker='o', \
                    s=2.+8.*np.log10(float(data[i][5])))
            else:
                plt.scatter([float(data[i][3])],[float(data[i][4])], \
                    color='black', marker='o', \
                    s=2.+8.*np.log10(float(data[i][5])))
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.savefig("clusterimage.png")
        plt.close(fig)
        
        import sys; sys.exit()




    ### For every star of the cluster+field, read their pairs of proper
    ### motion. (The positions and masses are not used in the fitting 
    ### of the probability distribution law.)
    mulist=[]
    for i in range(0,len(data)):
        mulist.append([float(data[i][1]),float(data[i][2])])

    ndim = 11
    nwalkers = nwalkers2
    
    ### Choose an initial set of positions for the walkers.
    p0 =    [
            [random.uniform(Mlim[0][0],Mlim[0][1]),     ### mux_cluster
            random.uniform(Mlim[1][0],Mlim[1][1]),      ### muy_cluster
            random.uniform(Mlim[2][0],Mlim[2][1]),      ### sigmax_cluster
            random.uniform(Mlim[3][0],Mlim[3][1]),      ### sigmay_cluster
            random.uniform(Mlim[4][0],Mlim[4][1]),      ### rho_cluster
            random.uniform(Mlim[5][0],Mlim[5][1]),      ### n_cluster
            random.uniform(Mlim[6][0],Mlim[6][1]),      ### mux_field
            random.uniform(Mlim[7][0],Mlim[7][1]),      ### muy_field
            random.uniform(Mlim[8][0],Mlim[8][1]),      ### sigmax_field
            random.uniform(Mlim[9][0],Mlim[9][1]),      ### sigmay_field
            random.uniform(Mlim[10][0],Mlim[10][1])]    ### rho_field
            for i in xrange(nwalkers)]


    lnprobability=[]
    ### Initialize the sampler with the chosen specs.
    argumentos=[mulist, Mlim]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, \
        ln_posterior_cluster, args=argumentos, a=2, threads=1)
    ### Run steps as a burn-in.
    print("PART 2: BURNIN (this may take some time!)...")
    pos, prob, state = sampler.run_mcmc(p0, nburnin2)
    lnprobability.append(sampler.lnprobability)
    ### Reset the chain to remove the burn-in samples.
    sampler.reset()
    ### Starting from the final position in the burn-in chain, 
    print("PART 2: SAMPLING (this may take some time!)...")
    pos, prob, state = sampler.run_mcmc(pos, nsampling2, rstate0=state)
    print("PART 2: SAMPLING DONE.")
    lnprobability.append(sampler.lnprobability)




    ### Writing the walkers in an external file.
    ### Each line contains one component of all of the walkers.
    ### The last line contains the probabilities.
    f0 = open(outputtitle, "w")
    [f0.write(str(sampler.flatchain[j,0])+" ") \
        for j in xrange(0,len(sampler.flatchain))];f0.write("\n")
    [f0.write(str(sampler.flatchain[j,1])+" ") \
        for j in xrange(0,len(sampler.flatchain))];f0.write("\n")
    [f0.write(str(sampler.flatchain[j,2])+" ") \
        for j in xrange(0,len(sampler.flatchain))];f0.write("\n")
    [f0.write(str(sampler.flatchain[j,3])+" ") \
        for j in xrange(0,len(sampler.flatchain))];f0.write("\n")
    [f0.write(str(sampler.flatchain[j,4])+" ") \
        for j in xrange(0,len(sampler.flatchain))];f0.write("\n")
    [f0.write(str(sampler.flatchain[j,5])+" ") \
        for j in xrange(0,len(sampler.flatchain))];f0.write("\n")
    [f0.write(str(sampler.flatchain[j,6])+" ") \
        for j in xrange(0,len(sampler.flatchain))];f0.write("\n")
    [f0.write(str(sampler.flatchain[j,7])+" ") \
        for j in xrange(0,len(sampler.flatchain))];f0.write("\n")
    [f0.write(str(sampler.flatchain[j,8])+" ") \
        for j in xrange(0,len(sampler.flatchain))];f0.write("\n")
    [f0.write(str(sampler.flatchain[j,9])+" ") \
        for j in xrange(0,len(sampler.flatchain))];f0.write("\n")
    [f0.write(str(sampler.flatchain[j,10])+" ") \
        for j in xrange(0,len(sampler.flatchain))];f0.write("\n")
    [f0.write(str(sampler.flatlnprobability[j])+" ") \
        for j in xrange(0,len(sampler.flatchain))];f0.write("\n")
    f0.close()


    if 1==1:
        poslnprob=np.arange(1,len(lnprobability[0][0])+1)
        fig=plt.figure(figsize=(6,6))
        ax=plt.subplot(1,1,1)
        
        for i in xrange(0,len(lnprobability[0])):
            plt.plot(poslnprob,[np.arcsinh(lnprobability[0][i][j]) for j in range(0,len(lnprobability[0][0]))], color='black', linewidth=0.05, linestyle='-')
        fig.savefig("probabilities_part2.png")


        samples = sampler.chain[:, :, :].reshape((-1, ndim))

        fig = corner.corner(samples,\
            labels=["$\\mu_{xc}$","$\\mu_{yc}$",\
            "$\\sigma_{xc}$","$\\sigma_{yc}$","$\\rho_{c}$","$n_{c}$",\
            "$\\mu_{xf}$","$\\mu_{yf}$","$\\sigma_{xf}$","$\\sigma_{yf}$",\
            "$\\rho_{f}$"],\
            truths=[mu_c[0],mu_c[1], sigma_c[0], sigma_c[1], rho_c,\
                float(N_c)/float(N_c+N_f), mu_f[0],mu_f[1], sigma_f[0],\
                sigma_f[1], rho_f, float(N_f)/float(N_c+N_f)])
        fig.savefig("triangle_part2.png")


















if part3 == True:
    
    outputfile3="example_steven_analysis.txt"



    f0 = open(outputfile3, "r")
    linhas=f0.readlines()
    f0.close()
    
    particles=[]
    for i in range(0,len(linhas)):
        linha=linhas[i].split()
        particles.append([float(linha[j]) for j in range(0,len(linha))])
    
    
    
    
    print("MAXIMUM LIKELIHOOD")
    argmaxprob=np.argmax(particles[11])
    mlparticles=np.array([particles[0][argmaxprob],\
                            particles[1][argmaxprob],\
                            particles[2][argmaxprob],\
                            particles[3][argmaxprob],\
                            particles[4][argmaxprob],\
                            particles[5][argmaxprob],\
                            particles[6][argmaxprob],\
                            particles[7][argmaxprob],\
                            particles[8][argmaxprob],\
                            particles[9][argmaxprob],\
                            particles[10][argmaxprob]])
    print("mux_cluster    = "+str(mlparticles[0]))
    print("muy_cluster    = "+str(mlparticles[1]))
    print("sigmax_cluster = "+str(mlparticles[2]))
    print("sigmay_cluster = "+str(mlparticles[3]))
    print("rho_cluster    = "+str(mlparticles[4]))
    print("n_cluster      = "+str(mlparticles[5]))
    print("")
    print("mux_field      = "+str(mlparticles[6]))
    print("muy_field      = "+str(mlparticles[7]))
    print("sigmax_field   = "+str(mlparticles[8]))
    print("sigmay_field   = "+str(mlparticles[9]))
    print("rho_field      = "+str(mlparticles[10]))

    print("")
    print("")
    
    print("MEDIAN")
    print("mux_cluster    = "+str(np.percentile(particles[0],50))+" +- ("+str(np.percentile(particles[0],84)-np.percentile(particles[0],50))+" , "+str(np.percentile(particles[0],50)-np.percentile(particles[0],16))+")")
    print("muy_cluster    = "+str(np.percentile(particles[1],50))+" +- ("+str(np.percentile(particles[1],84)-np.percentile(particles[1],50))+" , "+str(np.percentile(particles[1],50)-np.percentile(particles[1],16))+")")
    print("sigmax_cluster = "+str(np.percentile(particles[2],50))+" +- ("+str(np.percentile(particles[2],84)-np.percentile(particles[2],50))+" , "+str(np.percentile(particles[2],50)-np.percentile(particles[2],16))+")")
    print("sigmay_cluster = "+str(np.percentile(particles[3],50))+" +- ("+str(np.percentile(particles[3],84)-np.percentile(particles[3],50))+" , "+str(np.percentile(particles[3],50)-np.percentile(particles[3],16))+")")
    print("rho_cluster    = "+str(np.percentile(particles[4],50))+" +- ("+str(np.percentile(particles[4],84)-np.percentile(particles[4],50))+" , "+str(np.percentile(particles[4],50)-np.percentile(particles[4],16))+")")
    print("n_cluster      = "+str(np.percentile(particles[5],50))+" +- ("+str(np.percentile(particles[5],84)-np.percentile(particles[5],50))+" , "+str(np.percentile(particles[5],50)-np.percentile(particles[5],16))+")")
    print("")
    print("mux_field      = "+str(np.percentile(particles[6],50))+" +- ("+str(np.percentile(particles[6],84)-np.percentile(particles[6],50))+" , "+str(np.percentile(particles[6],50)-np.percentile(particles[6],16))+")")
    print("muy_field      = "+str(np.percentile(particles[7],50))+" +- ("+str(np.percentile(particles[7],84)-np.percentile(particles[7],50))+" , "+str(np.percentile(particles[7],50)-np.percentile(particles[7],16))+")")
    print("sigmax_field   = "+str(np.percentile(particles[8],50))+" +- ("+str(np.percentile(particles[8],84)-np.percentile(particles[8],50))+" , "+str(np.percentile(particles[8],50)-np.percentile(particles[8],16))+")")
    print("sigmay_field   = "+str(np.percentile(particles[9],50))+" +- ("+str(np.percentile(particles[9],84)-np.percentile(particles[9],50))+" , "+str(np.percentile(particles[9],50)-np.percentile(particles[9],16))+")")
    print("rho_field      = "+str(np.percentile(particles[10],50))+" +- ("+str(np.percentile(particles[10],84)-np.percentile(particles[10],50))+" , "+str(np.percentile(particles[10],50)-np.percentile(particles[10],16))+")")


    if 1==2:
        samples = np.zeros((len(particles[0]),11))
        for i in range(0,11):
            for j in range(0,len(particles[0])):
                samples[j,i]=particles[i][j]

        fig = corner.corner(samples,\
            labels=["$\\mu_{xc}$","$\\mu_{yc}$",\
            "$\\sigma_{xc}$","$\\sigma_{yc}$","$\\rho_{c}$","$n_{c}$",\
            "$\\mu_{xf}$","$\\mu_{yf}$","$\\sigma_{xf}$","$\\sigma_{yf}$",\
            "$\\rho_{f}$"])
        fig.savefig("triangle.png")



### Reading files from cluster and field and shuffling them:
    files3=["example_steven_cluster.txt","example_steven_field.txt"]

    data3=[]
    for ifiles in range(0,len(files3)):
        f1=open(files3[ifiles],"r")
        linhas=f1.readlines()
        f1.close()

        for i in range(0,len(linhas)):
            data3.append([linhas[i].split()[0],linhas[i].split()[1],\
                linhas[i].split()[2],linhas[i].split()[3],\
                linhas[i].split()[4],linhas[i].split()[5]])
            
    np.random.shuffle(data3)


    MP=[]
    outputMP="example_steven_MP.txt"
    f2=open(outputMP,"w")
    Mreducedc=np.array([mlparticles[0],mlparticles[1],mlparticles[2],\
            mlparticles[3],mlparticles[4],mlparticles[5]])
    Mreducedf=np.array([mlparticles[6],mlparticles[7],mlparticles[8],\
            mlparticles[9],mlparticles[10],1.-mlparticles[5]])
    for i in range(0,len(data3)):
        [f2.write(str(data3[i][j])+" ") for j in range(0,len(data3[i]))]
        muxnow=float(data3[i][1])
        muynow=float(data3[i][2])
        psic=likelihood_cluster(Mreducedc, muxnow, muynow)
        psif=likelihood_cluster(Mreducedf, muxnow, muynow)
        f2.write(" "+str(psic/(psic+psif)))
        f2.write("\n")
    f2.close()

    # TODO: histograms of probability for cluster and field

    f2=open(outputMP,"r")
    linhas=f2.readlines()
    f2.close()
    
    percent_cluster=[]
    percent_field=[]
    for i in range(0,len(linhas)):
        linha=linhas[i].split()
        if "cluster" in linha[0]:
            percent_cluster.append(float(linha[6]))
        else:
            percent_field.append(float(linha[6]))
    percent_cluster=np.array(percent_cluster)
    percent_field=np.array(percent_field)
    
    
    
    plt.hist(percent_cluster, bins=20, color="blue", alpha=0.5) 
    plt.hist(percent_field, bins=20, color="red", alpha=0.5)
    plt.title("Membership probability of stars")
    plt.show()
    
    
    



