"""
This is an under development program.
"""

import numpy as np
import pyhdust.phc as phc
import pyhdust.spectools as spt
import pyhdust.lrr as lrr
import matplotlib.pyplot as plt


#############################
### 

### full name of the file containing the observables generated from 
### HDUST's output.
fileobservables = "observables_BeAtlas.txt"

### The domain of OldBeAtlas: 
npar=['3.0','3.5','4.0','4.5']
sigpar=['0.00','0.02','0.05','0.12','0.28','0.68','1.65','4.00']
Mpar=['03.80','04.20','04.80','05.50','06.40','07.70','08.60',\
    '09.60','10.80','12.50','14.60']
obpar=['1.10','1.20','1.30','1.40','1.45']

cosipar=[   '-4.3711e-08','0.11147','0.22155','0.33381','0.44464',\
            '0.55484','0.66653','0.77824','0.88862','1.0']


### Reading the file with the observables generated from HDUST's output.
g0 = open(fileobservables,'r')
g0linhas = g0.readlines()
g0.close()
for ilinha in range(0,len(g0linhas)):
    g0linhas[ilinha] = g0linhas[ilinha].split()

### Declaring the lists that will receive the different observables. 
### This is a preliminary step, before making the arrays for each observable.
TEMP_read = []
SOURCE_read = []
SNRATIOS_read = []
UBVRI_read = []
JHK_read = []
HALPHA_SOAR_read = []
LINE_HALPHA_read = []
LINE_HBETA_read = []
LINE_HGAMMA_read = []
LINE_BRGAMMA_read = []
LINE_BRALPHA_read = []
LINE_PFGAMMA_read = []
LINE_HUMPHREY14_read = []
LINE_HUMPHREY15_read = []
LINE_HUMPHREY16_read = []
LINE_HUMPHREY18_read = []
LINE_HUMPHREY19_read = []
LINE_HUMPHREY20_read = []
LINE_HUMPHREY21_read = []
LINE_HUMPHREY22_read = []
LINE_HUMPHREY23_read = []
LINE_HUMPHREY24_read = []
LINE_HUMPHREY25_read = []
BL_FLUX_read = []
RL_FLUX_read = []


### Filling the lists with the data from the file
MODEL_key = 0
COSI_key = 0
for ilinha in range(0,len(g0linhas)):
    
    if g0linhas[ilinha][0] == "MODEL" and MODEL_key == 0:
        ### Entering a new model. Reading the model's parameters.
        current_MODEL = [g0linhas[ilinha][ii] for ii in range(1,5)]
        MODEL_key = 1
    
    if g0linhas[ilinha][0] == "SOURCE" and MODEL_key == 1:
        ### 'SOURCE_read' will receive a new element composed of 2 lists:
        ### 0: model parameters (not including cosi)
        ### 1: elements of 'SOURCE_read'
        current_SOURCE = [g0linhas[ilinha][ii] for ii in range(1,6)]
        SOURCE_read.append([current_MODEL,current_SOURCE])

    if g0linhas[ilinha][0] == "TEMP_R" and MODEL_key == 1:
        ### 'TEMP_read' will receive a new element composed of 3 lists:
        ### 0: model parameters (not including cosi)
        ### 1: radial distances in the plane of the disk (R/Req)
        ### 2: temperatures of the disk [K]
        current_TEMP_R = [g0linhas[ilinha][ii] \
            for ii in range(1,len(g0linhas[ilinha]))]
        current_TEMP_T = [g0linhas[ilinha+1][ii] \
            for ii in range(1,len(g0linhas[ilinha+1]))]
        TEMP_read.append([current_MODEL,current_TEMP_R,current_TEMP_T])
                        
    if g0linhas[ilinha][0] == "COSI" and MODEL_key == 1 and COSI_key == 0:
        ### For a certain model, entering a new inclination (cosi)
        current_COSI = g0linhas[ilinha][1]
        COSI_key = 1
        
    
    def reading_procedure(lista,linename,Nelems):

        if g0linhas[ilinha][0] == linename and COSI_key == 1:
            ### 'lista' will receive a new element composed of 2 lists:
            ### 0: model parameters (including cosi)
            ### 1: elements of 'lista'
            auxi = [current_MODEL[ii] for ii in range(0,len(current_MODEL))]
            auxi.append(current_COSI)
            auxi2 = [g0linhas[ilinha][ii] for ii in range(1,Nelems)]
            lista.append([auxi,auxi2])
        return

    reading_procedure(SNRATIOS_read,"SNRATIOS",len(g0linhas[ilinha]))
    
    reading_procedure(UBVRI_read,"UBVRI",6)
    reading_procedure(JHK_read,"JHK",4)
    reading_procedure(HALPHA_SOAR_read,"HALPHA_SOAR",2)
    
    reading_procedure(LINE_HALPHA_read,"LINE_HALPHA",6)
    reading_procedure(LINE_HBETA_read,"LINE_HBETA",6)
    reading_procedure(LINE_HGAMMA_read,"LINE_HGAMMA",6)    

    reading_procedure(LINE_BRGAMMA_read,"LINE_BRGAMMA",6)    
    
    reading_procedure(LINE_BRALPHA_read,"LINE_BRALPHA",6)    

    reading_procedure(LINE_PFGAMMA_read,"LINE_PFGAMMA",6)    

    reading_procedure(LINE_HUMPHREY14_read,"LINE_HUMPHREY14",6)    
    reading_procedure(LINE_HUMPHREY15_read,"LINE_HUMPHREY15",6)    
    reading_procedure(LINE_HUMPHREY16_read,"LINE_HUMPHREY16",6)    
    reading_procedure(LINE_HUMPHREY18_read,"LINE_HUMPHREY18",6)    
    reading_procedure(LINE_HUMPHREY19_read,"LINE_HUMPHREY19",6)    
    reading_procedure(LINE_HUMPHREY20_read,"LINE_HUMPHREY20",6)    
    reading_procedure(LINE_HUMPHREY21_read,"LINE_HUMPHREY21",6)    
    reading_procedure(LINE_HUMPHREY22_read,"LINE_HUMPHREY22",6)    
    reading_procedure(LINE_HUMPHREY23_read,"LINE_HUMPHREY23",6)    
    reading_procedure(LINE_HUMPHREY24_read,"LINE_HUMPHREY24",6)    
    reading_procedure(LINE_HUMPHREY25_read,"LINE_HUMPHREY25",6)    

    reading_procedure(BL_FLUX_read,"BL_FLUX",4)    
    reading_procedure(RL_FLUX_read,"RL_FLUX",4)    


    if g0linhas[ilinha][0] == "END_COSI" and COSI_key == 1:
        COSI_key = 0
        
    if g0linhas[ilinha][0] == "END_MODEL" and MODEL_key == 1:
        MODEL_key = 0







def attribution_procedure(lista_read,Nelems):
    
    array = np.zeros((len(npar),len(sigpar),len(Mpar),len(obpar),\
            len(cosipar),Nelems))
    array[:] = np.nan

    for i in range(0,len(lista_read)):
        idxnpar = npar.index(lista_read[i][0][0])
        idxsigpar = sigpar.index(lista_read[i][0][1])
        idxMpar = Mpar.index(lista_read[i][0][2])
        idxobpar = obpar.index(lista_read[i][0][3])
        idxcosipar = cosipar.index(lista_read[i][0][4])
    
        for j in range(0,Nelems):
            array[idxnpar,idxsigpar,idxMpar,idxobpar,idxcosipar,j] = \
                float(lista_read[i][1][j])
    return array

SNRATIOS = attribution_procedure(SNRATIOS_read,len(SNRATIOS_read[0][1]))

UBVRI = attribution_procedure(UBVRI_read,5)
JHK = attribution_procedure(JHK_read,3)
HALPHA_SOAR = attribution_procedure(HALPHA_SOAR_read,1)

LINE_HALPHA = attribution_procedure(LINE_HALPHA_read,5)
LINE_HBETA = attribution_procedure(LINE_HBETA_read,5)
LINE_HGAMMA = attribution_procedure(LINE_HGAMMA_read,5)

LINE_BRGAMMA = attribution_procedure(LINE_BRGAMMA_read,5)

LINE_BRALPHA = attribution_procedure(LINE_BRALPHA_read,5)

LINE_PFGAMMA = attribution_procedure(LINE_PFGAMMA_read,5)

LINE_HUMPHREY14 = attribution_procedure(LINE_HUMPHREY14_read,5)
LINE_HUMPHREY15 = attribution_procedure(LINE_HUMPHREY15_read,5)
LINE_HUMPHREY16 = attribution_procedure(LINE_HUMPHREY16_read,5)
LINE_HUMPHREY18 = attribution_procedure(LINE_HUMPHREY18_read,5)
LINE_HUMPHREY19 = attribution_procedure(LINE_HUMPHREY19_read,5)
LINE_HUMPHREY20 = attribution_procedure(LINE_HUMPHREY20_read,5)
LINE_HUMPHREY21 = attribution_procedure(LINE_HUMPHREY21_read,5)
LINE_HUMPHREY22 = attribution_procedure(LINE_HUMPHREY22_read,5)
LINE_HUMPHREY23 = attribution_procedure(LINE_HUMPHREY23_read,5)
LINE_HUMPHREY24 = attribution_procedure(LINE_HUMPHREY24_read,5)
LINE_HUMPHREY25 = attribution_procedure(LINE_HUMPHREY25_read,5)

BL_FLUX = attribution_procedure(BL_FLUX_read,3)
RL_FLUX = attribution_procedure(RL_FLUX_read,3)





#############################
### Central wavelength of each line [microns]
LINE_HALPHA_lbd = spt.hydrogenlinewl(3, 2)*1e6
LINE_HBETA_lbd = spt.hydrogenlinewl(4, 2)*1e6
LINE_HGAMMA_lbd = spt.hydrogenlinewl(5, 2)*1e6
LINE_BRGAMMA_lbd = spt.hydrogenlinewl(7, 4)*1e6
LINE_BRALPHA_lbd = spt.hydrogenlinewl(5, 4)*1e6
LINE_PFGAMMA_lbd = spt.hydrogenlinewl(8, 5)*1e6
LINE_HUMPHREY14_lbd = spt.hydrogenlinewl(14, 6)*1e6
LINE_HUMPHREY15_lbd = spt.hydrogenlinewl(15, 6)*1e6
LINE_HUMPHREY16_lbd = spt.hydrogenlinewl(16, 6)*1e6
LINE_HUMPHREY18_lbd = spt.hydrogenlinewl(18, 6)*1e6
LINE_HUMPHREY19_lbd = spt.hydrogenlinewl(19, 6)*1e6
LINE_HUMPHREY20_lbd = spt.hydrogenlinewl(20, 6)*1e6
LINE_HUMPHREY21_lbd = spt.hydrogenlinewl(21, 6)*1e6
LINE_HUMPHREY22_lbd = spt.hydrogenlinewl(22, 6)*1e6
LINE_HUMPHREY23_lbd = spt.hydrogenlinewl(23, 6)*1e6
LINE_HUMPHREY24_lbd = spt.hydrogenlinewl(24, 6)*1e6
LINE_HUMPHREY25_lbd = spt.hydrogenlinewl(25, 6)*1e6







##########################################################
##########################################################
### Now, comes the part 1 of the analysis: making lots of plots!


#############################
### TODO (Fredy): plots of the SEDs and line profiles
if 1==2:

    import read_everything
    
    ### Reading fullsed, source and temperature files
    files_fullsed_new, files_source_new, files_temps_new, fullsed_contents, \
            fullsed_path, source_path, temps_path, dist_std = \
            read_everything.read_everything()




#############################
### Plotting Lenorzer Diagrams
if 1==1:

    ### Parameters for the double arcsinh scaling:
    up1 = 1.
    up2 = 100.  ### since most of the measured flux ratios are between 
                ### 0.1 and 1, the choice of up2 = 100 garantees the 
                ### nearly logarithmic behaviour of the scale for them.
    down1 = 0.2 * up1   ### This makes a "compression" of the negative axis.
    down2 = 100.

    
    XL=[]
    YL=[]
    H14BL=[]
    figname=[]
    figtitle=[]
    annotate_vec=[]
    i_n = npar.index("4.5")
    i_sig = sigpar.index("1.65")
    i_M = Mpar.index("14.60")
    i_ob = obpar.index("1.40")
    i_cosi = cosipar.index("1.0")
    
    ###
    #kinit = 1
    #colorvec=[  "gray","black","brown","red","orange","green","blue","purple"]
    #for iplot1 in range(0,len(sigpar)):
    #    auxi_XL=[]
    #    auxi_YL=[]
    #    auxi_H14BL=[]
    #    auxi_annotate_vec=[]
    #    for iplot2 in range(0,len(cosipar)):
    #        auxi_XL.append(LINE_HUMPHREY14[i_n,iplot1,i_M,i_ob,iplot2,0]/\
    #            LINE_PFGAMMA[i_n,iplot1,i_M,i_ob,iplot2,0])
    #        auxi_YL.append(LINE_HUMPHREY14[i_n,iplot1,i_M,i_ob,iplot2,0]/\
    #            LINE_BRALPHA[i_n,iplot1,i_M,i_ob,iplot2,0])
    #        auxi_H14BL.append(LINE_HUMPHREY14[i_n,iplot1,i_M,i_ob,iplot2,0]/\
    #            BL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,0])
    #        auxi_annotate_vec.append("$\cos i$ = "+\
    #                str(round(abs(float(cosipar[iplot2])),2)))
    #    XL.append(auxi_XL)
    #    YL.append(auxi_YL)
    #    H14BL.append(auxi_H14BL)
    #    annotate_vec.append(auxi_annotate_vec)
    #    figname.append("n{0}_sig{1}_M{2}_ob{3}_cosi_X.png".format(\
    #        npar[i_n],sigpar[iplot1],Mpar[i_M],obpar[i_ob]))
    #    figtitle.append("$n="+npar[i_n]+"$, $\Sigma_0="+sigpar[iplot1]+\
    #        "\,\mathrm{g\,cm^{-2}}$")

    ###
    kinit = 0
    colorvec=["black","black","black","black","black","black","black",\
                    "black","black","black"]
    for iplot2 in range(0,len(cosipar)):
        auxi_XL=[]
        auxi_YL=[]
        auxi_H14BL=[]
        auxi_annotate_vec=[]
        for iplot1 in range(0,len(sigpar)):
            auxi_XL.append(LINE_HUMPHREY14[i_n,iplot1,i_M,i_ob,iplot2,0]/\
                LINE_PFGAMMA[i_n,iplot1,i_M,i_ob,iplot2,0])
            auxi_YL.append(LINE_HUMPHREY14[i_n,iplot1,i_M,i_ob,iplot2,0]/\
                LINE_BRALPHA[i_n,iplot1,i_M,i_ob,iplot2,0])
            auxi_H14BL.append(LINE_HUMPHREY14[i_n,iplot1,i_M,i_ob,iplot2,0]/\
                BL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,0])
            auxi_annotate_vec.append("$\Sigma_0$ = "+sigpar[iplot1])
        XL.append(auxi_XL)
        YL.append(auxi_YL)
        H14BL.append(auxi_H14BL)
        annotate_vec.append(auxi_annotate_vec)
        figname.append("n{0}_sigX.XX_M{1}_ob{2}_cosi_{3}.png".format(\
            npar[i_n],Mpar[i_M],obpar[i_ob],str(round(abs(float(cosipar[iplot2])),2)) ))
        figtitle.append("$n="+npar[i_n]+"$, $\cos i = "+\
                    str(round(abs(float(cosipar[iplot2])),2))+"$")


    ### 
    XL=[np.array(XL[k]) for k in range(0,len(XL))]
    YL=[np.array(YL[k]) for k in range(0,len(XL))]
    H14BL=[np.array(H14BL[k]) for k in range(0,len(XL))]
        
    ### Defining the labels of the axis
    axisvalsy=np.array([5.,2.,1.,0.5,0.2,0.1,0.05,0.02,0.01,0.,\
                            -0.05,-0.5,-5.,-50.,500.])
    axisvalsx=np.array([5.,2.,1.,0.5,0.2,0.1,0.05,0.01,0.,\
                            -0.5,-50.,500.])
    ### Obtaing the positions of the above defined labels of the axis
    transf_axisvalsy=np.array([lrr.scale_two_arcsinh(axisvalsy[i],\
            up1,up2,down1,down2) for i in xrange(0,len(axisvalsy))])
    transf_axisvalsx=np.array([lrr.scale_two_arcsinh(axisvalsx[i],\
            up1,up2,down1,down2) for i in xrange(0,len(axisvalsx))])

    ### Converting to the double arcsinh scale, for plotting:
    XLplot=[np.array([lrr.scale_two_arcsinh(XL[k][i],up1,up2,down1,down2) \
        for i in range(0,len(XL[k]))]) for k in range(0,len(XL))]
    YLplot=[np.array([lrr.scale_two_arcsinh(YL[k][i],up1,up2,down1,down2) \
        for i in range(0,len(YL[k]))]) for k in range(0,len(YL))]


    for k in range(kinit,len(XL)):
        
        plt.figure(k,figsize=(5.5, 5.5), dpi=100)

        plt.plot([0.,0.],[-1e32,1e32],linestyle=":",color="black",\
            linewidth=0.6)  ### Vertical dotted line
        plt.plot([-1e32,1e32],[0.,0.],linestyle=":",color="black",\
            linewidth=0.6)  ### Horizontal dotted line
    
        type1_x = -0.1
        type1_y = -0.2
        plt.plot([lrr.scale_two_arcsinh(10.**type1_x,\
                up1,up2,down1,down2),lrr.scale_two_arcsinh(10.**type1_x,\
                up1,up2,down1,down2)],\
                    [lrr.scale_two_arcsinh(10.**type1_y,\
                up1,up2,down1,down2),1e32],linestyle=":",color="blue",\
            linewidth=0.5)  ### Vertical line delimiting "type 1 region"
        plt.plot([lrr.scale_two_arcsinh(10.**type1_x,\
                up1,up2,down1,down2),1e32],[lrr.scale_two_arcsinh(10.**type1_y,\
                up1,up2,down1,down2),lrr.scale_two_arcsinh(10.**type1_y,\
                up1,up2,down1,down2)],linestyle=":",color="blue",\
            linewidth=0.5)  ### Horizontal line delimiting "type 1 region"
    
    
        plt.plot(XLplot[k],YLplot[k],linewidth=0.3,color=colorvec[k])
        for i in xrange(0,len(XLplot[k])):
            if H14BL[k][i] >= 0.:
                markk="o"
            else:
                markk="^"
            plt.scatter([XLplot[k][i]],[YLplot[k][i]],marker=markk,\
                s=1e4*np.abs(H14BL[k][i]),\
                color=colorvec[k],facecolors="None")
                #alpha=0.5+0.5*float(i)/float(len(XLplot[k])-1))
                #color="black",facecolors="none")
        
        i_init = -1; key_i_init = 0
        i_final = len(XLplot[k]); key_i_final = 0
        while i_init < len(XLplot[k])-1 and key_i_init == 0:
            i_init += 1
            if (XLplot[k][i_init] is not None) and \
                        (YLplot[k][i_init] is not None):
                if (~np.isnan(XLplot[k][i_init])) and \
                        (~np.isnan(YLplot[k][i_init])):
                    key_i_init = 1
        while i_final > 0 and key_i_final == 0:
            i_final -= 1
            if (XLplot[k][i_final] is not None) and \
                        (YLplot[k][i_final] is not None):
                if (~np.isnan(XLplot[k][i_final])) and \
                        (~np.isnan(YLplot[k][i_final])):
                    key_i_final = 1
        
        if (XLplot[k][i_init] is not None) and \
                    (YLplot[k][i_init] is not None):
            if (~np.isnan(XLplot[k][i_init])) and \
                    (~np.isnan(YLplot[k][i_init])):
                plt.text(XLplot[k][i_init],YLplot[k][i_init], \
                        annotate_vec[k][i_init])
        if (XLplot[k][i_final] is not None) and \
                    (YLplot[k][i_final] is not None):
            if (~np.isnan(XLplot[k][i_final])) and \
                    (~np.isnan(YLplot[k][i_final])):
                plt.text(XLplot[k][i_final],YLplot[k][i_final], \
                        annotate_vec[k][i_final])
            
    
        plt.scatter([1e32], [1e32],marker="o",s=1e4*0.01,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = 0.01$")
        plt.scatter([1e32], [1e32],marker="o",s=1e4*0.003,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = 0.003$")
        plt.scatter([1e32], [1e32],marker="o",s=1e4*0.001,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = 0.001$")
        plt.scatter([1e32], [1e32],marker="o",s=1e4*0.0003,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = 0.0003$")
        plt.scatter([1e32], [1e32],marker="^",s=1e4*0.0003,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = -0.0003$")
        plt.scatter([1e32], [1e32],marker="^",s=1e4*0.001,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = -0.001$")
        plt.scatter([1e32], [1e32],marker="^",s=1e4*0.003,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = -0.003$")
        plt.scatter([1e32], [1e32],marker="^",s=1e4*0.01,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = -0.01$")
        ### Add a legend
        plt.legend()

        plt.title(figtitle[k])
        plt.xscale('linear')
        plt.yscale('linear')
        plt.xlabel("$F(\mathrm{H}_{14})/F(\mathrm{Pf}_{\gamma})$")
        plt.ylabel("$F(\mathrm{H}_{14})/F(\mathrm{Br}_{\\alpha})$")
        plt.xticks(transf_axisvalsx,axisvalsx)
        plt.yticks(transf_axisvalsy,axisvalsy)
        plt.xlim([lrr.scale_two_arcsinh(-80.,\
                up1,up2,down1,down2),lrr.scale_two_arcsinh(8.,\
                up1,up2,down1,down2)])
        plt.ylim([lrr.scale_two_arcsinh(-80.,\
                up1,up2,down1,down2),lrr.scale_two_arcsinh(8.,\
                up1,up2,down1,down2)])
        
        plt.savefig("LENORZER_"+figname[k])






#############################
### Plotting Humphrey's diagrams
if 1==2:
    
    ### Parameters for the double arcsinh scaling:
    up1 = 1.
    up2 = 5.    ### since most of the measured flux ratios are between 
                ### 0.1 and 1, the choice of up2 = 100 garantees the 
                ### nearly logarithmic behaviour of the scale for them.
    down1 = 0.2 * up1   ### This makes a "compression" of the negative axis.
    down2 = up2
    
    H14H19=[]; EWL_H14=[]
    H15H19=[]; EWL_H15=[]
    H16H19=[]; EWL_H16=[]
    H18H19=[]; EWL_H18=[]
    H19H19=[]; EWL_H19=[]
    H20H19=[]; EWL_H20=[]
    H21H19=[]; EWL_H21=[]
    H22H19=[]; EWL_H22=[]
    H23H19=[]; EWL_H23=[]
    H24H19=[]; EWL_H24=[]
    H25H19=[]; EWL_H25=[]
    H14BL=[]
    figname=[]
    figtitle=[]
    annotate_vec=[]
    
    i_n = npar.index("3.0")
    i_sig = sigpar.index("1.65")
    i_M = Mpar.index("14.60")
    i_ob = obpar.index("1.40")
    i_cosi = cosipar.index("1.0")


    ### 
    kinit = 1
    colorvec=[  "gray","black","brown","red","orange","green","blue","purple"]
    for iplot1 in range(0,len(sigpar)):
        auxi_H14H19=[]; auxi_EWL_H14=[]
        auxi_H15H19=[]; auxi_EWL_H15=[]
        auxi_H16H19=[]; auxi_EWL_H16=[]
        auxi_H18H19=[]; auxi_EWL_H18=[]
        auxi_H19H19=[]; auxi_EWL_H19=[]
        auxi_H20H19=[]; auxi_EWL_H20=[]
        auxi_H21H19=[]; auxi_EWL_H21=[]
        auxi_H22H19=[]; auxi_EWL_H22=[]
        auxi_H23H19=[]; auxi_EWL_H23=[]
        auxi_H24H19=[]; auxi_EWL_H24=[]
        auxi_H25H19=[]; auxi_EWL_H25=[]
        auxi_H14BL=[]
        auxi_annotate_vec=[]
        for iplot2 in range(0,len(cosipar)):
            
            def making_fluxratio(auxilist,NUM,DEN):
                auxilist.append(\
                NUM[i_n,iplot1,i_M,i_ob,iplot2,0]/\
                DEN[i_n,iplot1,i_M,i_ob,iplot2,0])
                return
            
            making_fluxratio(auxi_H14H19,LINE_HUMPHREY14,LINE_HUMPHREY19)
            making_fluxratio(auxi_H15H19,LINE_HUMPHREY15,LINE_HUMPHREY19)
            making_fluxratio(auxi_H16H19,LINE_HUMPHREY16,LINE_HUMPHREY19)
            making_fluxratio(auxi_H18H19,LINE_HUMPHREY18,LINE_HUMPHREY19)
            making_fluxratio(auxi_H19H19,LINE_HUMPHREY19,LINE_HUMPHREY19)
            making_fluxratio(auxi_H20H19,LINE_HUMPHREY20,LINE_HUMPHREY19)
            making_fluxratio(auxi_H21H19,LINE_HUMPHREY21,LINE_HUMPHREY19)
            making_fluxratio(auxi_H22H19,LINE_HUMPHREY22,LINE_HUMPHREY19)
            making_fluxratio(auxi_H23H19,LINE_HUMPHREY23,LINE_HUMPHREY19)
            making_fluxratio(auxi_H24H19,LINE_HUMPHREY24,LINE_HUMPHREY19)
            making_fluxratio(auxi_H25H19,LINE_HUMPHREY25,LINE_HUMPHREY19)

            making_fluxratio(auxi_H14BL,LINE_HUMPHREY14,BL_FLUX)            

            def making_ew(auxilist,NUM,DEN):
                auxilist.append(\
                NUM[i_n,iplot1,i_M,i_ob,iplot2,1]/\
                DEN*1e-4)                
                return
            
            making_ew(auxi_EWL_H14,LINE_HUMPHREY14,LINE_HUMPHREY14_lbd)
            making_ew(auxi_EWL_H15,LINE_HUMPHREY15,LINE_HUMPHREY15_lbd)
            making_ew(auxi_EWL_H16,LINE_HUMPHREY16,LINE_HUMPHREY16_lbd)
            making_ew(auxi_EWL_H18,LINE_HUMPHREY18,LINE_HUMPHREY18_lbd)
            making_ew(auxi_EWL_H19,LINE_HUMPHREY19,LINE_HUMPHREY19_lbd)
            making_ew(auxi_EWL_H20,LINE_HUMPHREY20,LINE_HUMPHREY20_lbd)
            making_ew(auxi_EWL_H21,LINE_HUMPHREY21,LINE_HUMPHREY21_lbd)
            making_ew(auxi_EWL_H22,LINE_HUMPHREY22,LINE_HUMPHREY22_lbd)
            making_ew(auxi_EWL_H23,LINE_HUMPHREY23,LINE_HUMPHREY23_lbd)
            making_ew(auxi_EWL_H24,LINE_HUMPHREY24,LINE_HUMPHREY24_lbd)
            making_ew(auxi_EWL_H25,LINE_HUMPHREY25,LINE_HUMPHREY25_lbd)
                        
            auxi_annotate_vec.append("$\cos i$ = "+\
                    str(round(abs(float(cosipar[iplot2])),2)))
                    
        H14H19.append(auxi_H14H19); EWL_H14.append(auxi_EWL_H14)
        H15H19.append(auxi_H15H19); EWL_H15.append(auxi_EWL_H15)
        H16H19.append(auxi_H16H19); EWL_H16.append(auxi_EWL_H16)
        H18H19.append(auxi_H18H19); EWL_H18.append(auxi_EWL_H18)
        H19H19.append(auxi_H19H19); EWL_H19.append(auxi_EWL_H19)
        H20H19.append(auxi_H20H19); EWL_H20.append(auxi_EWL_H20)
        H21H19.append(auxi_H21H19); EWL_H21.append(auxi_EWL_H21)
        H22H19.append(auxi_H22H19); EWL_H22.append(auxi_EWL_H22)
        H23H19.append(auxi_H23H19); EWL_H23.append(auxi_EWL_H23)
        H24H19.append(auxi_H24H19); EWL_H24.append(auxi_EWL_H24)
        H25H19.append(auxi_H25H19); EWL_H25.append(auxi_EWL_H25)
        H14BL.append(auxi_H14BL)
        annotate_vec.append(auxi_annotate_vec)
        figname.append("n{0}_sig{1}_M{2}_ob{3}_cosi_X.png".format(\
            npar[i_n],sigpar[iplot1],Mpar[i_M],obpar[i_ob]))
        figtitle.append("$n="+npar[i_n]+"$, $\Sigma_0="+sigpar[iplot1]+\
            "\,\mathrm{g\,cm^{-2}}$")





    ### 
    H14H19=[np.array(H14H19[k]) for k in range(0,len(H14H19))]
    H15H19=[np.array(H15H19[k]) for k in range(0,len(H15H19))]
    H16H19=[np.array(H16H19[k]) for k in range(0,len(H16H19))]
    H18H19=[np.array(H18H19[k]) for k in range(0,len(H18H19))]
    H19H19=[np.array(H19H19[k]) for k in range(0,len(H19H19))]
    H20H19=[np.array(H20H19[k]) for k in range(0,len(H20H19))]
    H21H19=[np.array(H21H19[k]) for k in range(0,len(H21H19))]
    H22H19=[np.array(H22H19[k]) for k in range(0,len(H22H19))]
    H23H19=[np.array(H23H19[k]) for k in range(0,len(H23H19))]
    H24H19=[np.array(H24H19[k]) for k in range(0,len(H24H19))]
    H25H19=[np.array(H25H19[k]) for k in range(0,len(H25H19))]
    H14BL=[np.array(H14BL[k]) for k in range(0,len(H14BL))]
    EWL_H14=[np.array(EWL_H14[k]) for k in range(0,len(EWL_H14))]
    EWL_H15=[np.array(EWL_H15[k]) for k in range(0,len(EWL_H15))]
    EWL_H16=[np.array(EWL_H16[k]) for k in range(0,len(EWL_H16))]
    EWL_H18=[np.array(EWL_H18[k]) for k in range(0,len(EWL_H18))]
    EWL_H19=[np.array(EWL_H19[k]) for k in range(0,len(EWL_H19))]
    EWL_H20=[np.array(EWL_H20[k]) for k in range(0,len(EWL_H20))]
    EWL_H21=[np.array(EWL_H21[k]) for k in range(0,len(EWL_H21))]
    EWL_H22=[np.array(EWL_H22[k]) for k in range(0,len(EWL_H22))]
    EWL_H23=[np.array(EWL_H23[k]) for k in range(0,len(EWL_H23))]
    EWL_H24=[np.array(EWL_H24[k]) for k in range(0,len(EWL_H24))]
    EWL_H25=[np.array(EWL_H25[k]) for k in range(0,len(EWL_H25))]
        
    ### Defining the labels of the axis
    axisvalsy=np.array([10.,5.,3.,2.,1.,0.5,0.2,0.1,0.,\
                            -0.5,-5.])
    ### Obtaing the positions of the above defined labels of the axis
    transf_axisvalsy=np.array([lrr.scale_two_arcsinh(axisvalsy[i],\
            up1,up2,down1,down2) for i in xrange(0,len(axisvalsy))])

    for k in range(kinit,len(H14H19)):

        Dlamb = 0.03
        
        def making_lambdas(lbdcenter,array,Dlamb):
            return np.array([lbdcenter - \
                0.5*Dlamb + Dlamb*float(i)/float(len(array)-1) \
                for i in xrange(0,len(array))])
                
        H14lamb = making_lambdas(LINE_HUMPHREY14_lbd,H14H19[k],Dlamb)
        H15lamb = making_lambdas(LINE_HUMPHREY15_lbd,H15H19[k],Dlamb)
        H16lamb = making_lambdas(LINE_HUMPHREY16_lbd,H16H19[k],Dlamb)
        H18lamb = making_lambdas(LINE_HUMPHREY18_lbd,H18H19[k],Dlamb)
        H19lamb = making_lambdas(LINE_HUMPHREY19_lbd,H19H19[k],Dlamb)
        H20lamb = making_lambdas(LINE_HUMPHREY20_lbd,H20H19[k],Dlamb)
        H21lamb = making_lambdas(LINE_HUMPHREY21_lbd,H21H19[k],Dlamb)
        H22lamb = making_lambdas(LINE_HUMPHREY22_lbd,H22H19[k],Dlamb)
        H23lamb = making_lambdas(LINE_HUMPHREY23_lbd,H23H19[k],Dlamb)
        H24lamb = making_lambdas(LINE_HUMPHREY24_lbd,H24H19[k],Dlamb)
        H25lamb = making_lambdas(LINE_HUMPHREY25_lbd,H25H19[k],Dlamb)
        
        plt.figure(k,figsize=(11, 7), dpi=100)

        plt.subplot(211)
            
        ymax = 15
        horlines = np.array([float(-15+ii) for ii in range(0,2*ymax+2)])
        for elem in horlines:
            plt.plot([-1e32,1e32],
            [lrr.scale_two_arcsinh(elem,up1,up2,down1,down2),\
            lrr.scale_two_arcsinh(elem,up1,up2,down1,down2)],
            linestyle=":",color="black",\
                linewidth=0.6)  ### Horizontal dotted lines
        
        colorvec=["red","green","blue"]
        kcolor=0
        
        def plot_ratio_line(x,y):
            plt.plot(x,np.array([\
                lrr.scale_two_arcsinh(y[j],up1,up2,down1,down2) \
                for j in range(0,len(y))]),\
                linewidth=0.3,color=colorvec[kcolor%3])
            return
        
        plot_ratio_line(H14lamb,H14H19[k]); kcolor+=1
        plot_ratio_line(H15lamb,H15H19[k]); kcolor+=1
        plot_ratio_line(H16lamb,H16H19[k]); kcolor+=1
        plot_ratio_line(H18lamb,H18H19[k]); kcolor+=1
        plot_ratio_line(H19lamb,H19H19[k]); kcolor+=1
        plot_ratio_line(H20lamb,H20H19[k]); kcolor+=1
        plot_ratio_line(H21lamb,H21H19[k]); kcolor+=1
        plot_ratio_line(H22lamb,H22H19[k]); kcolor+=1
        plot_ratio_line(H23lamb,H23H19[k]); kcolor+=1
        plot_ratio_line(H24lamb,H24H19[k]); kcolor+=1
        plot_ratio_line(H25lamb,H25H19[k]); kcolor+=1
        
                
        colorvec=["red","green","blue"]
        kcolor=0
        
        def plot_ratio_scatter(x,y):
            for i in xrange(0,len(y)):
                if H14BL[k][i] >= 0.:
                    markk="o"
                else:
                    markk="^"
                plt.scatter([x[i]],\
                    [lrr.scale_two_arcsinh(y[i],up1,up2,down1,down2)],\
                    marker=markk,\
                    s=1e4*np.abs(H14BL[k][i]),\
                    color=colorvec[kcolor%3],facecolors="None")
            return
            
        plot_ratio_scatter(H14lamb,H14H19[k]); kcolor+=1
        plot_ratio_scatter(H15lamb,H15H19[k]); kcolor+=1
        plot_ratio_scatter(H16lamb,H16H19[k]); kcolor+=1
        plot_ratio_scatter(H18lamb,H18H19[k]); kcolor+=1
        plot_ratio_scatter(H19lamb,H19H19[k]); kcolor+=1
        plot_ratio_scatter(H20lamb,H20H19[k]); kcolor+=1
        plot_ratio_scatter(H21lamb,H21H19[k]); kcolor+=1
        plot_ratio_scatter(H22lamb,H22H19[k]); kcolor+=1
        plot_ratio_scatter(H23lamb,H23H19[k]); kcolor+=1
        plot_ratio_scatter(H24lamb,H24H19[k]); kcolor+=1
        plot_ratio_scatter(H25lamb,H25H19[k]); kcolor+=1
            
            
    
        plt.scatter([1e32], [1e32],marker="o",s=1e4*0.01,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = 0.01$")
        plt.scatter([1e32], [1e32],marker="o",s=1e4*0.003,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = 0.003$")
        plt.scatter([1e32], [1e32],marker="o",s=1e4*0.001,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = 0.001$")
        plt.scatter([1e32], [1e32],marker="o",s=1e4*0.0003,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = 0.0003$")
        plt.scatter([1e32], [1e32],marker="^",s=1e4*0.0003,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = -0.0003$")
        plt.scatter([1e32], [1e32],marker="^",s=1e4*0.001,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = -0.001$")
        plt.scatter([1e32], [1e32],marker="^",s=1e4*0.003,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = -0.003$")
        plt.scatter([1e32], [1e32],marker="^",s=1e4*0.01,\
            color="black",facecolors="none", \
            label="$F(\mathrm{H}_{14})/B_L = -0.01$")
        ### Add a legend
        #plt.legend()

        plt.subplot(212)

        plt.plot([-1e32,1e32],[0.,0.],linestyle=":",color="black",\
            linewidth=0.6)  ### Horizontal dotted line
            
        colorvec=["red","green","blue"]
        kcolor=0

        def plot_ew_line(x,y):
            plt.plot(x,y*1e4,linewidth=0.3,color=colorvec[kcolor%3])
            return

        plot_ew_line(H14lamb,EWL_H14[k]); kcolor+=1
        plot_ew_line(H15lamb,EWL_H15[k]); kcolor+=1
        plot_ew_line(H16lamb,EWL_H16[k]); kcolor+=1
        plot_ew_line(H18lamb,EWL_H18[k]); kcolor+=1
        plot_ew_line(H19lamb,EWL_H19[k]); kcolor+=1
        plot_ew_line(H20lamb,EWL_H20[k]); kcolor+=1
        plot_ew_line(H21lamb,EWL_H21[k]); kcolor+=1
        plot_ew_line(H22lamb,EWL_H22[k]); kcolor+=1
        plot_ew_line(H23lamb,EWL_H23[k]); kcolor+=1
        plot_ew_line(H24lamb,EWL_H24[k]); kcolor+=1
        plot_ew_line(H25lamb,EWL_H25[k]); kcolor+=1
        
            
        colorvec=["red","green","blue"]
        kcolor=0
        
        def plot_ew_scatter(x,y):
            for i in xrange(0,len(y)):
                if H14BL[k][i] >= 0.:
                    markk="o"
                else:
                    markk="^"
                plt.scatter([x[i]],[y[i]*1e4],marker=markk,\
                    s=1e4*np.abs(H14BL[k][i]),\
                    color=colorvec[kcolor%3],facecolors="None")
            return        
        
        plot_ew_scatter(H14lamb,EWL_H14[k]); kcolor+=1
        plot_ew_scatter(H15lamb,EWL_H15[k]); kcolor+=1
        plot_ew_scatter(H16lamb,EWL_H16[k]); kcolor+=1
        plot_ew_scatter(H18lamb,EWL_H18[k]); kcolor+=1
        plot_ew_scatter(H19lamb,EWL_H19[k]); kcolor+=1
        plot_ew_scatter(H20lamb,EWL_H20[k]); kcolor+=1
        plot_ew_scatter(H21lamb,EWL_H21[k]); kcolor+=1
        plot_ew_scatter(H22lamb,EWL_H22[k]); kcolor+=1
        plot_ew_scatter(H23lamb,EWL_H23[k]); kcolor+=1
        plot_ew_scatter(H24lamb,EWL_H24[k]); kcolor+=1
        plot_ew_scatter(H25lamb,EWL_H25[k]); kcolor+=1
        



        plt.subplot(211)
        plt.title(figtitle[k])
        plt.xscale('linear')
        plt.yscale('linear')
        plt.yticks(transf_axisvalsy,axisvalsy)
        plt.ylabel("$F(\mathrm{H}_{n})/F(\mathrm{H}_{19})$")
        plt.xlabel("$\lambda\,[\mathrm{\mu m}]$")
        plt.xlim([3.4,4.1])
        plt.ylim([lrr.scale_two_arcsinh(-15.,up1,up2,down1,down2),\
            lrr.scale_two_arcsinh(15.,up1,up2,down1,down2)])

        plt.subplot(212)
        plt.xscale('linear')
        plt.yscale('linear')
        plt.ylabel("$EW/\lambda\,[\\times 10^{4}]$")
        plt.xlabel("$\lambda\,[\mathrm{\mu m}]$")
        plt.xlim([3.4,4.1])
        plt.ylim([2.5,-2.5])
        
        plt.savefig("HUMPHREYS_"+figname[k])




#############################
### Plotting Mennickent-based CMD
if 1==2:

    XL=[]
    YL=[]
    H14BL=[]
    figname=[]
    figtitle=[]
    annotate_vec=[]
    i_n = npar.index("4.5")
    i_sig = sigpar.index("1.65")
    i_M = Mpar.index("04.20")
    i_ob = obpar.index("1.40")
    i_cosi = cosipar.index("1.0")


    ###
    lamb1_BL=3.41 ; lamb2_BL=3.47
    lamb1_RL=3.93 ; lamb2_RL=4.00
    Nnpts=50
    xlp,ylp = lrr.VEGA_spct("spct1")
            
    llamb=np.array([lamb1_BL+(lamb2_BL-lamb1_BL)/\
        float(Nnpts-1)*float(i) for i in range(0,Nnpts)])*1e4
    dllamb=np.array([llamb[i+1]-llamb[i] for i in range(0,Nnpts-1)])
    ylpf=np.array([lrr.interpLinND([llamb[i]],[xlp],ylp) \
        for i in range(0,Nnpts)])
    B_Vega=lrr.integrate_trapezia(ylpf,dllamb)
    #print(B_Vega,B_Vega*0.5*(lamb1_BL+lamb2_BL)*1e-4/phc.c.cgs/phc.h.cgs)
    ###
    kinit = 0
    colorvec=[  "gray","black","brown","red","orange","green","blue","purple"]
    figname_now="n{0}_sigX.XX_M{1}_ob{2}_cosi_X.png".format(\
            npar[i_n],Mpar[i_M],obpar[i_ob])
    for iplot1 in range(0,len(sigpar)):
        auxi_XL=[]
        auxi_YL=[]
        auxi_H14BL=[]
        auxi_annotate_vec=[]
        for iplot2 in range(0,len(cosipar)):
            auxi_XL.append(1.-np.log10(
                (BL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,0]/
                (BL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,2]-\
                BL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,1]))/\
                (RL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,0]/
                (RL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,2]-\
                RL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,1]))
                )/np.log10(\
                (RL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,1]+\
                RL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,2])/\
                (BL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,1]+\
                BL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,2])
                ))
            auxi_YL.append(-2.5*np.log10(
                BL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,0]/\
                B_Vega))
            auxi_H14BL.append(LINE_HUMPHREY14[i_n,iplot1,i_M,i_ob,iplot2,0]/\
                BL_FLUX[i_n,iplot1,i_M,i_ob,iplot2,0])
            auxi_annotate_vec.append("$\cos i$ = "+\
                    str(round(abs(float(cosipar[iplot2])),2)))
        XL.append(auxi_XL)
        YL.append(auxi_YL)
        H14BL.append(auxi_H14BL)
        annotate_vec.append(auxi_annotate_vec)
        figname.append("n{0}_sig{1}_M{2}_ob{3}_cosi_X.png".format(\
            npar[i_n],sigpar[iplot1],Mpar[i_M],obpar[i_ob]))
        figtitle.append("$n="+npar[i_n]+"$, $\Sigma_0="+sigpar[iplot1]+\
            "\,\mathrm{g\,cm^{-2}}$")



    ### 
    XL=[np.array(XL[k]) for k in range(0,len(XL))]
    YL=[np.array(YL[k]) for k in range(0,len(XL))]
    YLmin=np.nanmin([np.nanmin(YL[k]) for k in range(0,len(YL))])
    if Mpar[i_M] == "14.60": YLmin = -6.3; YLmax = -2.8
    if Mpar[i_M] == "04.20": YLmin = -2.0; YLmax = -0.5
    H14BL=[np.array(H14BL[k]) for k in range(0,len(XL))]

    plt.figure(0,figsize=(5.5, 5.5), dpi=100)

    #plt.plot([1.-np.log10(2.1*(lamb2_RL-lamb1_RL)/\
    #                (lamb2_BL-lamb1_BL))/np.log10((lamb2_RL+lamb1_RL)/\
    #                (lamb2_BL+lamb1_BL)),
    #                1.-np.log10(2.1*(lamb2_RL-lamb1_RL)/\
    #                (lamb2_BL-lamb1_BL))/np.log10((lamb2_RL+lamb1_RL)/\
    #                (lamb2_BL+lamb1_BL))],[-1e32,1e32])

    for k in range(kinit,len(XL)):
        

        #plt.plot([-1e32,1e32],[0.,0.],linestyle=":",color="black",\
        #    linewidth=0.6)  ### Horizontal dotted line
    
    
        plt.plot(XL[k],YL[k],linewidth=0.3,color=colorvec[k])
        for i in xrange(0,len(XL[k])):
            if H14BL[k][i] >= 0.:
                markk="o"
            else:
                markk="^"
            plt.scatter([XL[k][i]],[YL[k][i]],marker=markk,\
                s=1e4*np.abs(H14BL[k][i]),\
                color=colorvec[k],facecolors="None")
        


    #plt.title(figtitle[k])
    plt.xscale('linear')
    plt.yscale('linear')
    plt.xlabel("$\\alpha_L$")
    plt.ylabel("$M_{B_L}-m_{B_L}(\mathrm{Vega})$")
    plt.xlim([-3.2,-1.0])
    plt.ylim([YLmax,YLmin])
        
    plt.savefig("CMD_"+figname_now)

















if 1==2: 
    
    
    ### Parameters for the double arcsinh scaling:
    up1 = 1.
    up2 = 100.  ### since most of the measured flux ratios are between 
                ### 0.1 and 1, the choice of up2 = 100 garantees the 
                ### nearly logarithmic behaviour of the scale for them.
    down1 = 0.2 * up1   ### This makes a "compression" of the negative axis.
    down2 = 100.

    i_n = npar.index("4.5")
    i_sig = sigpar.index("1.65")
    i_M = Mpar.index("14.60")
    i_ob = obpar.index("1.40")
    i_cosi = cosipar.index("1.0")


    xx = np.array([float(cosipar[k]) for k in range(0,len(cosipar))])
    yy = np.array([float(sigpar[k]) for k in range(1,len(sigpar))])
    newyy = np.array([lrr.scale_two_arcsinh(yy[k],\
            1.,100.,0.2,100.) for k in range(0,len(yy))])



    values_H14_Pfg=[]
    values_H14_Bra=[]
    for i in range(0,len(cosipar)):
        for j in range(1,len(sigpar)):
            auxi_H14_Pfg = lrr.scale_two_arcsinh(\
                        LINE_HUMPHREY14[i_n,j,i_M,i_ob,i,0]/\
                        LINE_PFGAMMA[i_n,j,i_M,i_ob,i,0],\
                        up1,up2,down1,down2)
            auxi_H14_Bra = lrr.scale_two_arcsinh(\
                        LINE_HUMPHREY14[i_n,j,i_M,i_ob,i,0]/\
                        LINE_BRALPHA[i_n,j,i_M,i_ob,i,0],\
                        up1,up2,down1,down2)
            values_H14_Pfg.append(auxi_H14_Pfg)
            values_H14_Bra.append(auxi_H14_Bra)

    axis=[xx,newyy]


    linear_cosipar = np.array([xx[0]+(xx[-1]-xx[0])/(100.-1.)*float(i)\
                        for i in range(0,100)])
    linear_newsigpar = np.array([newyy[0]+(newyy[-1]-newyy[0])/(100.-1.)*\
                        float(i) for i in range(0,100)])

    XL = np.zeros((len(linear_cosipar),len(linear_newsigpar)))
    XL[:] = np.nan
    YL = np.zeros((len(linear_cosipar),len(linear_newsigpar)))
    YL[:] = np.nan


    XL_up = lrr.scale_two_arcsinh(1.0,up1,up2,down1,down2)
    XL_low = lrr.scale_two_arcsinh(0.01,up1,up2,down1,down2)
    YL_up = lrr.scale_two_arcsinh(1.0,up1,up2,down1,down2)
    YL_low = lrr.scale_two_arcsinh(0.01,up1,up2,down1,down2)
    for i in range(0,len(linear_cosipar)):
        for j in range(0,len(linear_newsigpar)):
            point = [linear_cosipar[i],linear_newsigpar[j]]
            auxi_XL = lrr.interpLinND(point,axis,values_H14_Pfg)
            auxi_YL = lrr.interpLinND(point,axis,values_H14_Bra)

            if auxi_XL > XL_up or auxi_XL < XL_low or \
                    auxi_YL > YL_up or auxi_YL < YL_low:
                auxi_XL = XL_low
                auxi_YL = YL_low

            #if auxi_XL > XL_up: auxi_XL = XL_up
            #if auxi_XL < XL_low: auxi_XL = XL_low

            #if auxi_YL > YL_up: auxi_YL = YL_up
            #if auxi_YL < YL_low: auxi_YL = YL_low
            XL[i,j] = auxi_XL
            YL[i,j] = auxi_YL 
    XL[0,0] = XL_low; XL[-1,-1] = XL_up
    YL[0,0] = YL_low; YL[-1,-1] = YL_up


    plt.imshow(XL.T,extent=[xx.min(), xx.max(),newyy.min(), newyy.max()], \
        aspect="auto",interpolation="bilinear",origin="lower")
    plt.colorbar()
    plt.xticks(xx,np.array([abs(round(xx[k],3)) for k in range(0,len(xx))]))
    plt.yticks(newyy,yy)
    plt.show()

    plt.imshow(YL.T,extent=[xx.min(), xx.max(),newyy.min(), newyy.max()], \
        aspect="auto",interpolation="bilinear",origin="lower")
    plt.colorbar()
    plt.xticks(xx,np.array([abs(round(xx[k],3)) for k in range(0,len(xx))]))
    plt.yticks(newyy,yy)
    plt.show()    
    
    








##########################################################
##########################################################
### Now, comes the part 2 of the analysis: MCMC fitting for comparison 
### of models and observations.








