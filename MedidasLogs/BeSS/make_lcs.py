import itertools as it
import glob as glob
import numpy as np
import matplotlib.pyplot as plt
import pyhdust.spectools as spt
import pyhdust.phc as phc
from scipy.optimize import curve_fit
import pyhdust.lrr.jdutil as jdu
import pyhdust.lrr as lrr


def gaussian_fit(x,sigma2,A):
    """
    This is a gaussian function centered on x = 0 and with area 'A' and 
    variance 'sigma2'.
    """
    return 0.+A/np.sqrt(2.*np.pi*sigma2)*np.exp(-0.5*(x-0.)**2./sigma2)
    
def linear_fit(x,A,B):
    """
    This is a simple linear function.
    """
    return A+B*x


folders_data = ["HD144/","HD4180/","HD5394/","HD6811/","HD11606/","HD20336/",\
        "HD23302/","HD23480/","HD23630/","HD23862/","HD187811/","HD191610/",\
        "HD193009/","HD194335/","HD194883/","HD195907/","HD197419/","HD200310/",\
        "HD202904/","HD204722/","HD208057/","HD210129/","HD217675/","HD217891/"]

Starnames =     ["10 Cas","$o$ Cas","$\\gamma$ Cas","$\\phi$ And",\
                "V777 Cas","Bk Cam","17 Tau","23 Tau","25 Tau",\
                "28 Tau","12 Vul","28 Cyg","V2113 Cyg","V2119 Cyg",\
                "V2120 Cyg","V2123 Cyg","V568 Cyg","60 Cyg",\
                "$\\upsilon$ Cyg","V2162 Cyg","16 Peg","25 Peg",\
                "$o$ And","$\\beta$ Psc"]

### These are the days, months and years to be used in the loop below
days = ["01","02","03","04","05","06","07","08","09","10",\
        "11","12","13","14","15","16","17","18","19","20",
        "21","22","23","24","25","26","27","28","29","30",
        "31"]
months = ["01","02","03","04","05","06","07","08","09","10","11","12"]
years = ["2007","2008","2009","2010","2011","2012","2013","2014",\
        "2015","2016","2017","2018"]

### 
dates = [years,months,days]
dates = [elem for elem in it.product(*dates)]

### 
dateLband = ["2017","10","01"]
dateWISE = ["2012","03","14"]

### Loop over the folders of data
for folder_data in folders_data:

    ### Obtaining files' names for this specific 'folder_data'
    files_list = glob.glob(folder_data+"*")

    ### 'files_contents' will receive all the data.
    files_contents = []
    ### Loop over all possible dates
    for idate in range(0,len(dates)):
        current_date = dates[idate][0]+dates[idate][1]+dates[idate][2]
        ### Loop over the files
        for ifile in range(0,len(files_list)):
            current_file = files_list[ifile]
            ### 
            if current_date in current_file:
                ### Julian date
                JD = jdu.date_to_jd(int(dates[idate][0]),\
                        int(dates[idate][1]),int(dates[idate][2]))
                ### Reading the file
                f0 = open(current_file,"r")
                linesf0 = f0.readlines()
                f0.close()
            
                ### Saving the table in 'array_data':
                array_data = np.zeros((len(linesf0),2))
                array_data[:,:] = np.nan
                for iline in range(0,len(linesf0)):
                    array_data[iline,0] = float(linesf0[iline].split()[0])
                    array_data[iline,1] = float(linesf0[iline].split()[1])

                ### Wavelength of Halpha [microns]
                lbc = spt.hydrogenlinewl(3, 2)*1e6

                ### Selecting region for evaluating "noise"
                if len(array_data[:,0]) >= 2 and \
                        np.nanmin(array_data[:,0]) <= lbc*1e4-40. and \
                        np.nanmax(array_data[:,0]) >= lbc*1e4+40.:
                    forerr_lambs = []
                    forerr_data = []
                    for itest in range(0,len(array_data[:,0])):
                        if lbc*1e4+30. <= array_data[itest,0] <= lbc*1e4+40.:
                            forerr_lambs.append(array_data[itest,0])
                            forerr_data.append(array_data[itest,1])

                ### If the data contains a region around the Halpha line 
                ### and a few lines for evaluating the "noise":
                if len(array_data[:,0]) >= 2 and \
                        np.nanmin(array_data[:,0]) <= lbc*1e4-40. and \
                        np.nanmax(array_data[:,0]) >= lbc*1e4+40. and \
                        len(forerr_lambs) >= 10:
                            
                            
                            
                    

                    files_contents.append([
                        current_file,
                        JD,
                        array_data,
                        "",  ### [3] (EW, errEW) [Angstroms]
                        "",  ### [4] (PS, errPS) [km/s]
                        ""]) ### [5] [(FWHM, errFWHM) [km/s], (A, errA) [km/s]]
            
                    #plt.plot(array_data[:,0],array_data[:,1])
                    #plt.show()


                    hwidth = 2000.
    
                    xlp=array_data[:,0]*1e-4 ### lambda [microns]
                    ylp=array_data[:,1]      ### BeSS's flux [unit = ???]



                    ### Returning an array of velocities [km/s]: 'xplott' 
                    ### and an array with the normalized flux: 'yplott'                    
                    xplott,yplott = spt.lineProf(np.array(forerr_lambs)*1e-4, \
                            np.array(forerr_data), lbc, hwidth=hwidth)
                    
                    popt, pcov = curve_fit(linear_fit, xplott, yplott)
                    ### Variance of spectra
                    var = 0.
                    for ifit in range(0,len(xplott)):
                        var += (yplott[ifit]-linear_fit(xplott[ifit],\
                                popt[0],popt[1]))*\
                                (yplott[ifit]-linear_fit(xplott[ifit],\
                                popt[0],popt[1]))
                    var = var/float(len(xplott))



                    ### Returning an array of velocities [km/s]: 'xplot' 
                    ### and an array with the normalized flux: 'yplot'
                    xplot,yplot = spt.lineProf(xlp, ylp, lbc, hwidth=hwidth)
                    ### Equivalent width [Angstroms]
                    ew_ = spt.EWcalc(xplot, yplot, vw=hwidth)
                    ew = ew_*lbc/phc.c.cgs*1e9
                    
                    err_ew = 2.*hwidth*lbc/phc.c.cgs*1e9*2.*np.sqrt(var)
                    
                    
                    files_contents[-1][3] = (ew,err_ew)
            
                    ### Try to calculate the peak separation in [km/s]
                    try:
                        v1,v2 = spt.PScalc(xplot, yplot, vc=0.0, ssize=0.05, \
                                        gaussfit=True)
                    except:
                        v1 = np.nan; v2 = np.nan
                    
                    files_contents[-1][4] = (v2-v1,np.sqrt(var)*abs(v2-v1))
        
                    ### Trial of calculating the FWHM: A gaussian is ajusted to 
                    ### the absolute value of the line profile. The FWHM of this 
                    ### gaussian is extracted [km/s]. 
                    ### Also, the area [km/s] is extracted.
                    try:
                        popt, pcov = curve_fit(gaussian_fit, xplot, abs(yplot-1.),\
                                p0=[1000.,0.])
                        fwhm = np.sqrt(8.*popt[0]*np.log(2))
                        area = popt[1]
                        #plt.plot(xplot,abs(yplot-1.),color="black")
                        #plt.plot(xplot,gaussian_fit(xplot,popt[0],popt[1]),\
                        #        color="red")
                        #plt.show()
                    except:
                        fwhm = np.nan
                        area = np.nan
                    
                    files_contents[-1][5] = [(fwhm,np.sqrt(var)*abs(fwhm)),\
                            (area,np.sqrt(var)*abs(area))]

    if 1==1:

        ### 
        yearlabel = ["2007","2009","2011","2013","2015","2017","2019"]
        yearpos = [jdu.date_to_jd(int(yearlabel[iy]),1,1) \
                for iy in range(0,len(yearlabel))]




        plt.figure(figsize=(11,2))
        ### 
        plt.subplot(111)
        xlc = [files_contents[ifile][1] for ifile \
                in range(0,len(files_contents))]
        ylc = [files_contents[ifile][3][0] for ifile \
                in range(0,len(files_contents))]
        errylc = [files_contents[ifile][3][1] for ifile \
                in range(0,len(files_contents))]
        
        xLband = jdu.date_to_jd(int(dateLband[0]),\
                int(dateLband[1]),int(dateLband[2]))
        
        xWISE = jdu.date_to_jd(int(dateWISE[0]),int(dateWISE[1]),\
                int(dateWISE[2]))
        
        ylims = [10.,-60.]
        plt.plot([xLband,xLband],ylims,color = "brown")
        plt.plot([xWISE,xWISE],ylims,color = "green")
        
        plt.plot(xlc,ylc,color = "black", linewidth = 0.5)
        for ilc in range(0,len(xlc)):
            plt.errorbar(xlc[ilc],ylc[ilc],yerr = errylc[ilc],color = "black")
        
        plt.xticks(yearpos,yearlabel)
        plt.ylabel("$EW\,[\mathrm{A}]$")
        plt.ylim(ylims)
        idx_fd = folders_data.index(folder_data)
        plt.title(folder_data.replace("/","").replace("HD","HD ")+\
                " ("+Starnames[idx_fd]+")")
    

        #plt.subplot(212)
        #xlc = [files_contents[ifile][1] for ifile \
        #        in range(0,len(files_contents))]
        #ylc = [files_contents[ifile][4][0] for ifile \
        #        in range(0,len(files_contents))]
        #errylc = [files_contents[ifile][4][1] for ifile \
        #        in range(0,len(files_contents))]
        #plt.plot(xlc,ylc,color = "red", linewidth = 0.5)
        #for ilc in range(0,len(xlc)):
        #    plt.errorbar(xlc[ilc],ylc[ilc],yerr = errylc[ilc],color = "red")
        #
        #
        #xllc = [files_contents[ifile][1] for ifile \
        #        in range(0,len(files_contents))]
        #yllc = [files_contents[ifile][5][0][0] for ifile \
        #        in range(0,len(files_contents))]
        #ylims = [0.,2500.]
        #erryllc = [files_contents[ifile][5][0][1] for ifile \
        #        in range(0,len(files_contents))]
        #plt.plot([xLband,xLband],ylims,color = "brown")
        #plt.plot([xWISE,xWISE],ylims,color = "green")
        #plt.plot(xllc,yllc,color = "blue", linewidth = 0.5)
        #for ilc in range(0,len(xllc)):
        #    plt.errorbar(xllc[ilc],yllc[ilc],yerr = erryllc[ilc],color = "blue")
        #plt.xticks(yearpos,yearlabel)
        #plt.ylabel("PS / FWHM [$\mathrm{km\,s^{-1}}$]")
        #plt.ylim(ylims)


        ### 
        plt.tight_layout()
        plt.savefig(folder_data.replace("/","")+".png")




    if 1==1:
    
    
        
        times = [files_contents[idate][1] for idate \
                in range(0,len(files_contents))]
                
        EWs = [files_contents[idate][3][0] for idate \
                in range(0,len(files_contents))]
        errEWs = [files_contents[idate][3][1] for idate \
                in range(0,len(files_contents))]        

        PSs = [files_contents[idate][4][0] for idate \
                in range(0,len(files_contents))]
        errPSs = [files_contents[idate][4][1] for idate \
                in range(0,len(files_contents))]                

        FWHMs = [files_contents[idate][5][0][0] for idate \
                in range(0,len(files_contents))]
        errFWHMs = [files_contents[idate][5][0][1] for idate \
                in range(0,len(files_contents))]                        
        

        xLband = jdu.date_to_jd(int(dateLband[0]),\
                int(dateLband[1]),int(dateLband[2]))
        
        xWISE = jdu.date_to_jd(int(dateWISE[0]),int(dateWISE[1]),\
                int(dateWISE[2]))
    
    
        EW_tLb = lrr.interpLinND([xLband],[times],EWs,tp="linear")
        errEW_tLb = lrr.interpLinND([xLband],[times],errEWs,tp="linear")
        
        
        f0 = open(folder_data.replace("/","")+".dat","w")
        f0.write("### [0] times [JD]\n")
        f0.write("### [1] EW [Angstrom]\n")
        f0.write("### [2] errEW [Angstrom]\n")
        f0.write("### [3] Peak separation [km/s]\n")
        f0.write("### [4] error Peak separation [km/s]\n")
        f0.write("### [5] FWHM [km/s]\n")
        f0.write("### [6] error FWHM [km/s]\n")
        for idate in range(0,len(times)):
            f0.write(str(times[idate])+" "+\
                    str(EWs[idate])+" "+\
                    str(errEWs[idate])+" "+\
                    str(PSs[idate])+" "+\
                    str(errPSs[idate])+" "+\
                    str(FWHMs[idate])+" "+\
                    str(errFWHMs[idate])+"\n")
        f0.write("##### time in L-band [JD] ; EW [Angstrom] ; errEW [Angstrom] \n")
        f0.write(str(xLband)+" "+str(EW_tLb)+" "+str(errEW_tLb)+"\n")
        f0.close()
    
    



        f1 = open(folder_data.replace("/","")+"_LC.dat","w")
        
        for idate in range(0,len(times)):
            ### 
            f1.write(str(times[idate])+"\n")
            ### 
            lineswrite = []
            [lineswrite.append(el) for el in files_contents[idate][2][:,0]]
            [f1.write(str(el)+" ") for el in lineswrite]
            f1.write("\n")
            ### 
            lineswrite = []
            [lineswrite.append(el) for el in files_contents[idate][2][:,1]]
            [f1.write(str(el)+" ") for el in lineswrite]
            f1.write("\n")
        
        f1.close()
 
 
 
 
 
 
        
    
    
    
    
    

