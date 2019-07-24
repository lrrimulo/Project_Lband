"""

"""

import glob as glob
import numpy as np
import pyhdust as hdt
import pyhdust.phc as phc
import pyhdust.lrr as lrr
import pyhdust.lrr.roche_singlestar as rss
import pyhdust.spectools as spt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt






###########################
import read_everything

files_fullsed_new, files_source_new, files_temps_new, fullsed_contents, \
        fullsed_path, source_path, temps_path, dist_std = \
        read_everything.read_everything()

outputfile = "observables_BeAtlas.txt"





#############################
### Obtaining S/N ratios for some intervals

def powerlaw(x,C,alpha):
    """
    Assuming a power-law for the fitting of a portion of the SED.
    """
    return C*x**alpha

SNratios=[]
lbda = [3.41,3.93,0.50]
lbdb = [3.47,4.00,0.55]

for ifile in xrange(0,len(files_fullsed_new)):
    print("Obtaining S/N ratios for file")
    print(files_fullsed_new[ifile][0])

    contents=fullsed_contents[ifile][1]

    auxiSNratios=[]
    for incs in range(0,len(contents[1])):
        xlp=contents[2][incs,:,0]   ### lambda [microns]
        ylp=contents[2][incs,:,1]   ### HDUST's flux [microns^-1]

        if len(lbda) == len(lbdb) and len(lbda) > 0:
            auxi2SNratios=[]
            for ilbda in range(0,len(lbda)):
                idx=[]
                for ii in range(0,len(xlp)):
                    if lbda[ilbda] <= xlp[ii] <= lbdb[ilbda]:
                        idx.append(ii)
                idx=np.array(idx)
                
                xlp_ajust = np.array(xlp[idx])
                ylp_ajust = np.array(ylp[idx])
                
        
                if len(idx) > 1:
                    popt, pcov = curve_fit(powerlaw, xlp_ajust, ylp_ajust)
                    ylp_ajusted = powerlaw(xlp_ajust,popt[0],popt[1])
                    diff2 = (ylp_ajust-ylp_ajusted)**2.
                    noise = np.sqrt(diff2.mean())
                    signal = powerlaw(0.5*(lbda[ilbda]+lbdb[ilbda]),\
                                popt[0],popt[1])
                    auxi2SNratios.append([lbda[ilbda],lbdb[ilbda],signal/noise])
                else:
                    auxi2SNratios.append([lbda[ilbda],lbdb[ilbda],np.nan])
                    
        auxiSNratios.append(auxi2SNratios)
    SNratios.append(auxiSNratios)





#############################
### Obtaining line observables

### Central lambda [microns], width [km/s]
### WARNING!: If you add or remove lines from the lists below, 
### you will have to change the writing to an external file 
### at the end of this program.
lbc_Bralpha=[4.052248,1000.]
lbc_Pfgamma=[3.740536,1000.]
lbc_Brgamma=[2.1661178,1000.]
lbc_Humphreyoth=    [
                    [4.020843,1000.,"Humphrey14"],
                    [3.907525,1000.,"Humphrey15"],
                    [3.819428,1000.,"Humphrey16"],
                    #[3.749370,1000.,"Humphrey17"],  ### This line merges with Pfgamma
                    [3.692611,1000.,"Humphrey18"],
                    [3.645901,1000.,"Humphrey19"],
                    [3.606946,1000.,"Humphrey20"],
                    [3.574082,1000.,"Humphrey21"],
                    [3.546079,1000.,"Humphrey22"],
                    [3.522003,1000.,"Humphrey23"],
                    [3.501142,1000.,"Humphrey24"],
                    [3.482938,1000.,"Humphrey25"]   ### The last Humphrey 
                                                    ### transition calculated 
                                                    ### by HDUST
                    ]
lbc_Bracketoth= [
                [0.656461,1000.,"Halpha"],
                [0.486271,1000.,"Hbeta"],
                [0.434169,1000.,"Hgamma"]
                #[0.410289,1000.,"Hdelta"]
                ]


### Limits of the BL and RL bands, as defined by 
### Mennickent et al. 2009PASP..121..125M
lamb1_BL=3.41 ; lamb2_BL=3.47
lamb1_RL=3.93 ; lamb2_RL=4.00
    
Bralpha=[]
Pfgamma=[]
Brgamma=[]
Humphreys=[]
Brackets=[]
BL=[]
RL=[]


def gaussian_fit(x,sigma2,A):
    return 0.+A/np.sqrt(2.*np.pi*sigma2)*np.exp(-0.5*(x-0.)**2./sigma2)
    
def double_gaussian_fit(x,meanV,sigma2V,V,meanR,sigma2R,R):
    return 0.+\
        abs(V)/np.sqrt(2.*np.pi*sigma2V)*np.exp(-0.5*(x-meanV)**2./sigma2V)+\
        abs(R)/np.sqrt(2.*np.pi*sigma2R)*np.exp(-0.5*(x-meanR)**2./sigma2R)

def obtaining_flux_ew_PS(contents,lbc,hwidth):
    """
    A procedure for calculating the flux, EW and peak separation of
    a line whose center is 'lbc'.
    """
    
    print("obtaining_flux_ew_PS, for ",lbc)
    
    auxilinflux=[]
    auxiew=[]
    auxiPS=[]
    auxiFWHM=[]
    for incs in range(0,len(contents[1])):
        xlp=contents[2][incs,:,0]   ### lambda [microns]
        ylp=contents[2][incs,:,1]   ### HDUST's flux [microns^-1]
  
        vels = (xlp - lbc) / lbc * phc.c.cgs * 1e-5 ### Doppler vels 
                                                    ### in km/s
        ### line flux in erg/s cm^2
        linflux=spt.absLineCalc(vels, ylp, vw=hwidth)
        linflux=linflux * contents[4][3]*phc.Lsun.cgs/4./np.pi/\
                    (dist_std*phc.pc.cgs)**2.*1e5*lbc/phc.c.cgs
        auxilinflux.append(linflux)

        ### Returning an array of velocities 
        ### and an array with the normalized flux (which will not be 
        ### used in the calculations below)
        xplot,yplot=spt.lineProf(xlp, ylp, lbc, hwidth=hwidth)
        ### Equivalent width [Angstroms]
        ew=spt.EWcalc(xplot, yplot, vw=hwidth)
        ew = ew*lbc/phc.c.cgs*1e9
        auxiew.append(ew)
            
        ### Try to calculate the peak separation in [km/s]
        try:
            v1,v2=spt.PScalc(xplot, yplot, vc=0.0, ssize=0.05, \
                            gaussfit=True)
        except:
            v1=np.nan; v2=np.nan
        auxiPS.append(v2-v1)
        
        ### My trial of calculating the FWHM: A gaussian is ajusted to 
        ### the absolute value of the line profile. The FWHM of this 
        ### gaussian is extracted.
        try:
            popt, pcov = curve_fit(gaussian_fit, xplot, abs(yplot-1.),p0=[1000.,0.])
            fwhm = np.sqrt(8.*popt[0]*np.log(2))
            area = popt[1]
            #plt.plot(xplot,abs(yplot-1.),color="black")
            #plt.plot(xplot,gaussian_fit(xplot,popt[0],popt[1]),color="red")
            #plt.show()
        except:
            fwhm = np.nan
            area = np.nan
        auxiFWHM.append([fwhm,area])
            
    auxilinflux=np.array(auxilinflux)
    auxiew=np.array(auxiew)
    auxiPS=np.array(auxiPS)
    auxiFWHM=np.array(auxiFWHM)
    #print(auxiFWHM)

    return auxilinflux,auxiew,auxiPS,auxiFWHM



for ifile in xrange(0,len(files_fullsed_new)):
    print("Obtaining line observables for file")
    print(files_fullsed_new[ifile][0])

    ### Line observables for the Humphrey's
    Humphreysauxi=[]
    for iHump in range(0,len(lbc_Humphreyoth)):
        
        hwidth=lbc_Humphreyoth[iHump][1]
        lbc=lbc_Humphreyoth[iHump][0]
        contents=fullsed_contents[ifile][1]
            
        auxilinflux,auxiew,auxiPS,auxiFWHM = obtaining_flux_ew_PS(contents,lbc,hwidth)
            
        Humphreysauxi.append([auxilinflux,auxiew,auxiPS,auxiFWHM])
    Humphreys.append(Humphreysauxi)

    ### Line observables for the Bracket's
    Bracketsauxi=[]
    for iHump in range(0,len(lbc_Bracketoth)):
        
        hwidth=lbc_Bracketoth[iHump][1]
        lbc=lbc_Bracketoth[iHump][0]
        contents=fullsed_contents[ifile][1]
            
        auxilinflux,auxiew,auxiPS,auxiFWHM = obtaining_flux_ew_PS(contents,lbc,hwidth)
            
        Bracketsauxi.append([auxilinflux,auxiew,auxiPS,auxiFWHM])
    Brackets.append(Bracketsauxi)
        

    ### Line observables for the Br gamma
    hwidth=lbc_Brgamma[1]
    lbc=lbc_Brgamma[0]
    contents=fullsed_contents[ifile][1]
            
    auxilinflux,auxiew,auxiPS,auxiFWHM = obtaining_flux_ew_PS(contents,lbc,hwidth)
            
    Brgamma.append([auxilinflux,auxiew,auxiPS,auxiFWHM])
        
    ### Line observables for the Br alpha    
    hwidth=lbc_Bralpha[1]
    lbc=lbc_Bralpha[0]
    contents=fullsed_contents[ifile][1]
            
    auxilinflux,auxiew,auxiPS,auxiFWHM = obtaining_flux_ew_PS(contents,lbc,hwidth)
            
    Bralpha.append([auxilinflux,auxiew,auxiPS,auxiFWHM])

    ### Line observables for the Pf gamma
    hwidth=lbc_Pfgamma[1]
    lbc=lbc_Pfgamma[0]
    contents=fullsed_contents[ifile][1]
            
    auxilinflux,auxiew,auxiPS,auxiFWHM = obtaining_flux_ew_PS(contents,lbc,hwidth)
            
    Pfgamma.append([auxilinflux,auxiew,auxiPS,auxiFWHM])
        
    ### Obtaining the B and R fluxes defined by Mennickent et al.
    print("Obtaining the B and R fluxes defined by Mennickentet al.")
    contents=fullsed_contents[ifile][1]
    BLfluxauxi=[]
    RLfluxauxi=[]
    for incs in range(0,len(contents[1])):
        xlp=contents[2][incs,:,0]
        ylp=contents[2][incs,:,1]
        Nnpts=50    ### this number must be >= 3. A good choice is 50.
            
        llamb=np.array([lamb1_BL+(lamb2_BL-lamb1_BL)/\
            float(Nnpts-1)*float(i) for i in range(0,Nnpts)])
        dllamb=np.array([llamb[i+1]-llamb[i] for i in range(0,Nnpts-1)])
        ylpf=np.array([lrr.interpLinND([llamb[i]],[xlp],ylp) \
            for i in range(0,Nnpts)])
        BLflux=lrr.integrate_trapezia(ylpf,dllamb)
        BLflux=BLflux * contents[4][3]*phc.Lsun.cgs/4./np.pi/\
                        (dist_std*phc.pc.cgs)**2.
        BLfluxauxi.append(BLflux)
            
        llamb=np.array([lamb1_RL+(lamb2_RL-lamb1_RL)/\
            float(Nnpts-1)*float(i) for i in range(0,Nnpts)])
        dllamb=np.array([llamb[i+1]-llamb[i] for i in range(0,Nnpts-1)])
        ylpf=np.array([lrr.interpLinND([llamb[i]],[xlp],ylp) \
            for i in range(0,Nnpts)])
        RLflux=lrr.integrate_trapezia(ylpf,dllamb)
        RLflux=RLflux * contents[4][3]*phc.Lsun.cgs/4./np.pi/\
                        (dist_std*phc.pc.cgs)**2.
        RLfluxauxi.append(RLflux)
        
    BLfluxauxi=np.array(BLfluxauxi)    
    BL.append(BLfluxauxi)
    RLfluxauxi=np.array(RLfluxauxi)    
    RL.append(RLfluxauxi)        
        








#############################
### Obtaining magnitudes


filters=[   'bess-u','bess-b','bess-v','bess-r','bess-i',\
            'bess-j','bess-h','bess-k',\
            'filt_Ha'
        ]

npts_interp=50 ### this number must be >= 3. A good choice is 50.

    
print("Obtaining zero point constants...")
zp=[]
for j in xrange(0,len(filters)):
    zp.append(lrr.obtain_pogson_zp('spct1',filters[j],\
                npts=npts_interp))
print("")


all_photflux=[]
all_Mag=[]
for ifile in xrange(0,len(files_fullsed_new)):
        
    print("Obtaining magnitudes for file")
    print(files_fullsed_new[ifile][0])
        
    fullsedtest=files_fullsed_new[ifile][1]
    sourcetest=files_source_new[ifile]


    photflux_vec=[]
    for j in xrange(0,len(filters)):
        print("Obtaining photon fluxes for filter "+str(filters[j]))
        mu,lamb,flambda,photflux=lrr.fullsed2photonflux(fullsedtest,\
            sourcetest,filters[j],npts=npts_interp,dist=10.)
        photflux_vec.append(photflux)
    all_photflux.append(photflux_vec)

    Mag_vec=[]
    for j in xrange(0,len(filters)):
        Mag_vec.append(lrr.pogson(photflux_vec[j],zp[j]))
    all_Mag.append(Mag_vec)

    print("")






#############################
### Obtaining vsini
            
all_vsini=[]
for ifile in xrange(0,len(files_fullsed_new)):
    contents=fullsed_contents[ifile][1]
    stelpars=[elem for elem in contents[4]]
    rpolenow=stelpars[1]
    massnow=stelpars[0]
    Wnow=stelpars[2]
    lixo,omeganow,lixo,Wnow=rss.rocheparams(Wnow,"W")
    veqfile=rss.cte_veq(rpolenow,massnow,omeganow,1.0)
        
    auxiifile=[]
    for iobs in xrange(0,len(mu)):
        auxiifile.append(veqfile*(1.-mu[iobs]**2.)**0.5)
    all_vsini.append(np.array(auxiifile))


#############################
### Writing in the external file
print("Writing in the external file")


f0=open(outputfile,"w")

for ifile in xrange(0,len(files_fullsed_new)):
    f0.write("MODEL "+\
            str(files_fullsed_new[ifile][0][0])+" "+\
            str(files_fullsed_new[ifile][0][1])+" "+\
            str(files_fullsed_new[ifile][0][2])+" "+\
            str(files_fullsed_new[ifile][0][3])+"\n")
            
    contents=fullsed_contents[ifile][1]
    
    f0.write("    SOURCE "+\
            str(contents[4][0])+" "+\
            str(contents[4][1])+" "+\
            str(contents[4][2])+" "+\
            str(contents[4][3])+" "+\
            str(contents[4][4])+"\n")
            
    f0.write("    TEMP_R ")
    for ii in range(0,len(contents[6][0,:])):
        f0.write(str(contents[6][0,ii])+" ")
    f0.write("\n")
    f0.write("    TEMP_T ")
    for ii in range(0,len(contents[6][1,:])):
        f0.write(str(contents[6][1,ii])+" ")
    f0.write("\n")
    
    for incs in range(0,len(contents[1])):
        f0.write("    COSI "+\
            str(contents[1][incs])+"\n")

        elements = []
        for ii in range(0,len(SNratios[ifile][incs])):
            elements.append(SNratios[ifile][incs][ii][0])
            elements.append(SNratios[ifile][incs][ii][1])
            elements.append(SNratios[ifile][incs][ii][2])
        f0.write("        SNRATIOS ")
        for elem in elements:
            f0.write(str(elem)+" ")
        f0.write("\n")
        
        f0.write("        UBVRI "+\
            str(all_Mag[ifile][0][incs])+" "+\
            str(all_Mag[ifile][1][incs])+" "+\
            str(all_Mag[ifile][2][incs])+" "+\
            str(all_Mag[ifile][3][incs])+" "+\
            str(all_Mag[ifile][4][incs])+"\n")
        f0.write("        JHK "+\
            str(all_Mag[ifile][5][incs])+" "+\
            str(all_Mag[ifile][6][incs])+" "+\
            str(all_Mag[ifile][7][incs])+"\n")
        f0.write("        HALPHA_SOAR "+\
            str(all_Mag[ifile][8][incs])+"\n")
        
        f0.write("        LINE_HALPHA "+\
            str(Brackets[ifile][0][0][incs])+" "+\
            str(Brackets[ifile][0][1][incs])+" "+\
            str(Brackets[ifile][0][2][incs])+" "+\
            str(Brackets[ifile][0][3][incs][0])+" "+\
            str(Brackets[ifile][0][3][incs][1])+"\n")
        f0.write("        LINE_HBETA "+\
            str(Brackets[ifile][1][0][incs])+" "+\
            str(Brackets[ifile][1][1][incs])+" "+\
            str(Brackets[ifile][1][2][incs])+" "+\
            str(Brackets[ifile][1][3][incs][0])+" "+\
            str(Brackets[ifile][1][3][incs][1])+"\n")
        f0.write("        LINE_HGAMMA "+\
            str(Brackets[ifile][2][0][incs])+" "+\
            str(Brackets[ifile][2][1][incs])+" "+\
            str(Brackets[ifile][2][2][incs])+" "+\
            str(Brackets[ifile][2][3][incs][0])+" "+\
            str(Brackets[ifile][2][3][incs][1])+"\n")

        f0.write("        LINE_BRGAMMA "+\
            str(Brgamma[ifile][0][incs])+" "+\
            str(Brgamma[ifile][1][incs])+" "+\
            str(Brgamma[ifile][2][incs])+" "+\
            str(Brgamma[ifile][3][incs][0])+" "+\
            str(Brgamma[ifile][3][incs][1])+"\n")

        f0.write("        LINE_BRALPHA "+\
            str(Bralpha[ifile][0][incs])+" "+\
            str(Bralpha[ifile][1][incs])+" "+\
            str(Bralpha[ifile][2][incs])+" "+\
            str(Bralpha[ifile][3][incs][0])+" "+\
            str(Bralpha[ifile][3][incs][1])+"\n")
        f0.write("        LINE_PFGAMMA "+\
            str(Pfgamma[ifile][0][incs])+" "+\
            str(Pfgamma[ifile][1][incs])+" "+\
            str(Pfgamma[ifile][2][incs])+" "+\
            str(Pfgamma[ifile][3][incs][0])+" "+\
            str(Pfgamma[ifile][3][incs][1])+"\n")
        f0.write("        LINE_HUMPHREY14 "+\
            str(Humphreys[ifile][0][0][incs])+" "+\
            str(Humphreys[ifile][0][1][incs])+" "+\
            str(Humphreys[ifile][0][2][incs])+" "+\
            str(Humphreys[ifile][0][3][incs][0])+" "+\
            str(Humphreys[ifile][0][3][incs][1])+"\n")
        f0.write("        LINE_HUMPHREY15 "+\
            str(Humphreys[ifile][1][0][incs])+" "+\
            str(Humphreys[ifile][1][1][incs])+" "+\
            str(Humphreys[ifile][1][2][incs])+" "+\
            str(Humphreys[ifile][1][3][incs][0])+" "+\
            str(Humphreys[ifile][1][3][incs][1])+"\n")
        f0.write("        LINE_HUMPHREY16 "+\
            str(Humphreys[ifile][2][0][incs])+" "+\
            str(Humphreys[ifile][2][1][incs])+" "+\
            str(Humphreys[ifile][2][2][incs])+" "+\
            str(Humphreys[ifile][2][3][incs][0])+" "+\
            str(Humphreys[ifile][2][3][incs][1])+"\n")
        f0.write("        LINE_HUMPHREY18 "+\
            str(Humphreys[ifile][3][0][incs])+" "+\
            str(Humphreys[ifile][3][1][incs])+" "+\
            str(Humphreys[ifile][3][2][incs])+" "+\
            str(Humphreys[ifile][3][3][incs][0])+" "+\
            str(Humphreys[ifile][3][3][incs][1])+"\n")
        f0.write("        LINE_HUMPHREY19 "+\
            str(Humphreys[ifile][4][0][incs])+" "+\
            str(Humphreys[ifile][4][1][incs])+" "+\
            str(Humphreys[ifile][4][2][incs])+" "+\
            str(Humphreys[ifile][4][3][incs][0])+" "+\
            str(Humphreys[ifile][4][3][incs][1])+"\n")
        f0.write("        LINE_HUMPHREY20 "+\
            str(Humphreys[ifile][5][0][incs])+" "+\
            str(Humphreys[ifile][5][1][incs])+" "+\
            str(Humphreys[ifile][5][2][incs])+" "+\
            str(Humphreys[ifile][5][3][incs][0])+" "+\
            str(Humphreys[ifile][5][3][incs][1])+"\n")
        f0.write("        LINE_HUMPHREY21 "+\
            str(Humphreys[ifile][6][0][incs])+" "+\
            str(Humphreys[ifile][6][1][incs])+" "+\
            str(Humphreys[ifile][6][2][incs])+" "+\
            str(Humphreys[ifile][6][3][incs][0])+" "+\
            str(Humphreys[ifile][6][3][incs][1])+"\n")
        f0.write("        LINE_HUMPHREY22 "+\
            str(Humphreys[ifile][7][0][incs])+" "+\
            str(Humphreys[ifile][7][1][incs])+" "+\
            str(Humphreys[ifile][7][2][incs])+" "+\
            str(Humphreys[ifile][7][3][incs][0])+" "+\
            str(Humphreys[ifile][7][3][incs][1])+"\n")
        f0.write("        LINE_HUMPHREY23 "+\
            str(Humphreys[ifile][8][0][incs])+" "+\
            str(Humphreys[ifile][8][1][incs])+" "+\
            str(Humphreys[ifile][8][2][incs])+" "+\
            str(Humphreys[ifile][8][3][incs][0])+" "+\
            str(Humphreys[ifile][8][3][incs][1])+"\n")
        f0.write("        LINE_HUMPHREY24 "+\
            str(Humphreys[ifile][9][0][incs])+" "+\
            str(Humphreys[ifile][9][1][incs])+" "+\
            str(Humphreys[ifile][9][2][incs])+" "+\
            str(Humphreys[ifile][9][3][incs][0])+" "+\
            str(Humphreys[ifile][9][3][incs][1])+"\n")
        f0.write("        LINE_HUMPHREY25 "+\
            str(Humphreys[ifile][10][0][incs])+" "+\
            str(Humphreys[ifile][10][1][incs])+" "+\
            str(Humphreys[ifile][10][2][incs])+" "+\
            str(Humphreys[ifile][10][3][incs][0])+" "+\
            str(Humphreys[ifile][10][3][incs][1])+"\n")
        f0.write("        BL_FLUX "+\
            str(BL[ifile][incs])+" "+\
            str(lamb1_BL)+" "+\
            str(lamb2_BL)+"\n")
        f0.write("        RL_FLUX "+\
            str(RL[ifile][incs])+" "+\
            str(lamb1_RL)+" "+\
            str(lamb2_RL)+"\n")
        
        f0.write("    END_COSI \n")
    f0.write("END_MODEL \n")

f0.close()





        
        
        
        
