
import glob as glob
import numpy as np
import pyhdust as hdt
import pyhdust.phc as phc
import pyhdust.lrr as lrr
import pyhdust.lrr.roche_singlestar as rss


def read_everything():

    fullsed_path = '../OldBeAtlas/fullsed_v2/'
    source_path = '../OldBeAtlas/source/'
    temps_path = '../OldBeAtlas/temperatures/'
    dist_std = 10.  ### assumed distance [parsecs] for the calculations


    ###########################
    
    ### The domain of OldBeAtlas:
    npar=['3.0','3.5','4.0','4.5']
    sigpar=['0.00','0.02','0.05','0.12','0.28','0.68','1.65','4.00']
    Mpar=['03.80','04.20','04.80','05.50','06.40','07.70','08.60',\
        '09.60','10.80','12.50','14.60']
    obpar=['1.10','1.20','1.30','1.40','1.45']

    filepars=[npar,sigpar,Mpar,obpar]


    print("Reading the files...")
    print("")

    files_fullsed=sorted(glob.glob(fullsed_path+'*'))	
    files_source=sorted(glob.glob(source_path+'*'))
    files_temps=sorted(glob.glob(temps_path+'*'))

    files_fullsed_new=[]    ### will receive the names of the fullsed
                            ### files to be opened.

    ### It is assumed that the names of the fullsed files are of the form:
    ### fullsed_mod191_PLn4.0_sig0.05_h072_Rd050.0_Be_M04.80_ob1.10_H0.30_Z0.014_bE_Ell.sed2
    ### or
    ### fullsed_mod01_PLn3.5_sig0.00_h060_Rd050.0_Be_M03.80_ob1.20_H0.77_Z0.014_bE_Ell.sed2
    for i in range(0,len(npar)):
        for j in range(0,len(sigpar)):
            for k in range(0,len(Mpar)):
                for l in range(0,len(obpar)):
                    for ifile in xrange(0,len(files_fullsed)):
                        if ('PLn{0}_sig{1}_h072_Rd050.0_Be_'\
                            .format(filepars[0][i],filepars[1][j])+\
                            'M{0}_ob{1}_H0.30_Z0.014_bE_Ell'\
                            .format(filepars[2][k],filepars[3][l]) in \
                            files_fullsed[ifile]) \
                            or ('PLn{0}_sig{1}_h060_Rd050.0_Be_'\
                            .format(filepars[0][i],filepars[1][j])+\
                            'M{0}_ob{1}_H0.30_Z0.014_bE_Ell'\
                            .format(filepars[2][k],filepars[3][l]) in \
                            files_fullsed[ifile]):
                                files_fullsed_new.append([[ filepars[0][i],\
                                                            filepars[1][j],\
                                                            filepars[2][k],\
                                                            filepars[3][l]],\
                                                        files_fullsed[ifile]]) 


    ### It is assumed that the names of the source files are of the form:
    ### Be_M03.40_ob1.45_H0.54_Z0.014_bE_Ell.txt
    ### (Notice that the it is contained in the name of the fullsed file.)
    files_source_new=[] ### will receive the names of the source
                        ### files to be opened.
    for iffn in xrange(0,len(files_fullsed_new)):
        for ifs in xrange(0,len(files_source)):
            if files_source[ifs].replace(source_path,'').replace('.txt','')\
                        in files_fullsed_new[iffn][1]:
                files_source_new.append(files_source[ifs])
                #print(files_source_new[-1])

    ### It is assumed that the names of the temperature files are of the form:
    ### mod126_PLn3.5_sig0.28_h072_Rd050.0_Be_M09.60_ob1.20_H0.30_Z0.014_bE_Ell30_avg.temp
    ### (Notice that the it is contained in the name of the fullsed file.)
    files_temps_new=[]  ### will receive the names of the temperature
                        ### files to be opened.
    for iffn in xrange(0,len(files_fullsed_new)):
        achei=0 ### Some fullsed files may not have correspondent temp files,
                ### like the ones of purely photospherical models. 
        for ifs in xrange(0,len(files_temps)):
            if files_temps[ifs].replace(temps_path,'').replace(\
                    '30_avg.temp','')\
                in files_fullsed_new[iffn][1]:
                files_temps_new.append(files_temps[ifs])
                #print(files_temps_new[-1])
                achei=1
        if achei == 0:
            files_temps_new.append('EMPTY')
            #print(files_temps_new[-1])


    ### 

    fullsed_contents=[] ### This list will receive the important contents
                        ### of all the files
    for ifile in xrange(0,len(files_fullsed_new)):

        ### Reading the fullsed, source and temperature files:
    
        fullsedtest=files_fullsed_new[ifile][1]
        f0=open(fullsedtest,'r')
        f0linhas=f0.readlines()
        f0.close()

        sourcetest=files_source_new[ifile]
        f1=open(sourcetest,'r')
        f1linhas=f1.readlines()
        f1.close()    

        tempstest=files_temps_new[ifile]
        if tempstest != 'EMPTY':
            ### OBS: This pyhdust procedure will print 
            ### "'FILE' completely read!"
            ncr, ncmu, ncphi, nLTE, nNLTE, Rstarz, Raz, betaz, dataz, \
                pcr, pcmu, pcphi = hdt.readtemp(tempstest)
            abttemp=[
                    [dataz[0,i,ncmu/2,0]/Rstarz for i in \
                            xrange(0,len(dataz[0,:,ncmu/2,0]))],
                    [dataz[3,i,ncmu/2,0] for i in \
                            xrange(0,len(dataz[3,:,ncmu/2,0]))]
                    ]
        else:
            abttemp=[
                    [np.nan,np.nan],
                    [np.nan,np.nan]
                    ]


        ### Obtaining each element of the 'fullsed_contents' list

        nobs=int(f0linhas[3].split()[1])
        nlbd=int(f0linhas[3].split()[0])
        contents=[    
            fullsedtest,                    ### 0: Name of fullsed file
            np.zeros(nobs),                 ### 1: will receive the cosi's
            np.zeros((nobs,nlbd,3)),        ### 2: will receive the SED
            sourcetest,                     ### 3: Name of source file
            np.zeros(5),                    ### 4: will receive the 
                                            ###    parameters of the star 
                                            ###    (source)
            tempstest,                      ### 5: Name of temperature file
            np.zeros((2,len(abttemp[0]))),  ### 6: will receive the temp 
                                            ###    profile
            [[],[]]
            ]
        contents[1][:] = np.nan
        contents[2][:] = np.nan
        contents[4][:] = np.nan
        contents[6][:] = np.nan


        ### Receiving cosi and SED ("1" and "2")
        for iobs in xrange(0,nobs):
            mu = float(f0linhas[5+iobs*nlbd].split()[0])
            contents[1][iobs] = mu
            for ilbd in xrange(0, nlbd):
                auxi = f0linhas[5+iobs*nlbd+ilbd].split()
                contents[2][iobs, ilbd, 0] = float(auxi[2])
                contents[2][iobs, ilbd, 1] = float(auxi[3])
                contents[2][iobs, ilbd, 2] = float(auxi[7])


        ### Receiving parameters of the star (source) ("4")
        contents[4][0] = float(f1linhas[3].split()[2]) ### M
        contents[4][1] = float(f1linhas[4].split()[2]) ### R_pole
        contents[4][2] = float(f1linhas[5].split()[2]) ### W
        contents[4][3] = float(f1linhas[6].split()[2]) ### L
        contents[4][4] = float(f1linhas[7].split()[2]) ### Beta_GD
    
        ### Receiving the temperature profile ("6")
        for i in xrange(0,len(contents[6][0,:])):
            contents[6][0,i] = abttemp[0][i]
            contents[6][1,i] = abttemp[1][i]
        

        fullsed_contents.append([files_fullsed_new[ifile][0],contents])
        print(fullsed_contents[-1][0])

    print("")

    return files_fullsed_new, files_source_new, files_temps_new, fullsed_contents, \
        fullsed_path, source_path, temps_path, dist_std




