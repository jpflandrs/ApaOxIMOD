"""
                                     ApaOxIMOD
                    Modeling the evolution of Apatite Oxygen Isotopes
                                        A
Monte-Carlo simulation of the alteration of the original oxygen isotope composition of biogenic apatites

                        JP Flandrois* and C Lecuyer**
                                    2023  mk03

* Author and maintainer
** Equations and original idea
----
See the corresponding publication:
"Mitigation of the diagenesis risk in biological apatite δ18O interpretation"
Christophe Lécuyer(1) and Jean-Pierre Flandrois(2)
(1)LGL-TPE, CNRS UMR5276, ENSL, Univ Lyon, Univ Lyon 1, Institut Universitaire de France, 43 bd du 11 Novembre 1918, 69622 Villeurbanne, France.
(2)LBBE, CNRS UMR5558, Univ Lyon, Univ Lyon 1, 43 bd du 11 Novembre 1918, 69622 Villeurbanne, France.
----
The program is made available under the [CeCILL2.1](http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.txt) licence.
----
How it works :
The simulation tool ApaOxIMOD is basically a Monte Carlo process simulation, based on the creation of a huge array of A/W values taken from an uniform distribution between given threshold values. A moving A/W window may be set within the threshold values to study the influence of modifications of the A/W environment. This represents the diversity of the deposit situation and physical state of the biological apatite. Typically, 500,000 situations are simulated. For each situation equation (4) is used to compute the δ18O in biological apatite δ18O(A)f at the equilibrium on the whole array given the temperature T, initial δ18O in biological apatite (δ18O(A)i) and the δ18O in water δ18O(W)i. The four parameters (T, d18OWi, d18OAi) may be set and their variations may be explored by steps for each simulated A/W environment to analyze the relative influence of each parameter. The output is a family of histograms summarizing the possible fate of δ18O(A)i in the population after the diagenetic episode.

How it is written :
    The program is entirely relying on the array structure of NumPy to take advantage of the speed and efficiency of NumPy in manipulating large arrays. The statistics are also computing with NumPy function and python3 is minimally used to connect Nympy functions and for reporting. The 15 possible parameters and options are taken from a pilot file using the YAML format that requires only a plain-text editor to be quickly changed.
    Note that python3.11 is 10% more rapid than the previous versions, you may consider using it or future versions (>3.11).
"""

import argparse, pickle, os,sys
import yaml #yaml format for the parameters
import numpy as np
import matplotlib.pyplot as plt
from numpy import mean
from numpy import std
import time
import matplotlib.pyplot as plt
import numpy as np




def coreProgram(outputDir,arraysize,T,Tstep,Tnbsteps,d18OWi,d18OWistep,d18OWinbsteps,d18OAi,d18OAistep,d18OAinbsteps,WAlow,WAhigh,WA_window,WA_window_slices,analyticalProcessSimulation,sigmaLab,sigmaLabstep,sigmaLabnbsteps,threadsCount,instructionAndBrowserPath,figsize,bins,userTitle,d19OAf_range,typeOfFile,jobdescription):
     

    if args.Model:
        data = arraysFiller(T,d18OWi,d18OAi,arraysize,WAlow,WAhigh,analyticalProcessSimulation,sigmaLab)
        sentenceC=""
        sentenceD=""
        plt.subplots(nrows=1, ncols=1, figsize=(5,5),constrained_layout=True)
        plt.hist(data,bins=bins)
        plt.xlabel('δ18OAf')
        plt.ylabel('count')
        plt.title('Density of δ18OA final values by classes', fontsize=10)
        plt.suptitle("Fixed parameters δ18OAf values\nC. Lecuyer and JP. Flandrois 2023\n"+userTitle, fontsize=14)
        fileIs=os.path.join(outputDir,"ApaOxIS_model"+"~"+jobdescription+"."+typeOfFile)
        #plt.show()
        plt.savefig(fileIs)
    if args.TwinAW:
        INIT1=d18OAi
        STEP1=d18OAistep
        NB1=d18OAinbsteps
        sentenceA="  "
        sentenceB="  "
        sentenceC= "δ18OAi and δ18OWi "
        sentenceD=""
        summaryTag=""
        VAR="δ18OAi "
        print ("APATITE+WATER")
        INIT2=d18OWi
        STEP2=d18OWistep
        NB2=d18OWinbsteps
    if args.TwinTW: #two parameters varying + A/W by ranges
        INIT1=T
        STEP1=Tstep
        NB1=Tnbsteps
        sentenceA="  "
        sentenceB="  "
        VAR=" T "
        sentenceD=""
        summaryTag=""
        sentenceC= "T and δ18OWi "
        print ("TEMP+WATER")
        INIT2=d18OWi
        STEP2=d18OWistep
        NB2=d18OWinbsteps
        
    if args.Temperature:
        INIT=T
        STEP=Tstep
        NB=Tnbsteps
        sentenceA=" with T= "
        sentenceB=" and T= "
        sentenceC=" T "
        sentenceD="Temperature"
        summaryTag="T"
    if args.Water:
        print ('WATER')
        INIT=d18OWi
        STEP=d18OWistep
        NB=d18OWinbsteps
        sentenceA=" with δ18OWi= "
        sentenceB=" and δ18OWi= "
        sentenceC=" δ18OWi "
        sentenceD="δ18OWi"
        summaryTag="d18OWi"
    if args.Sigma:#sigmaLab,sigmaLabstep,sigmaLabnbsteps,
        INIT=sigmaLab
        STEP=sigmaLabstep
        NB=sigmaLabnbsteps
        sentenceA=" with σ = "
        sentenceB=" and σ = "
        sentenceC=" σ "
        sentenceD=" σ "
        summaryTag="sigma"
    if args.Apatite:#sigmaLab,sigmaLabstep,sigmaLabnbsteps,
        INIT=d18OAi
        STEP=d18OAistep
        NB=d18OAinbsteps
        sentenceA=" with δ18OAi= "
        sentenceB=" and δ18OAi= "
        sentenceC=" δ18OAi "
        sentenceD="δ18OAi"
        summaryTag="d18OAi"
    if args.Temperature or args.Sigma or args.Water or args.Apatite:
        PROCESS=INIT
        print (PROCESS,INIT+(NB*STEP))
        
        ncols = 2
        nrows = NB+1
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize,constrained_layout=True)
        #print (fig,axes)#Figure(1000x1000) [<AxesSubplot:> <AxesSubplot:> <AxesSubplot:>]
        #fig.tight_layout()
        linescpt=0
        kolmncpt=0
        while abs(PROCESS) <= abs(INIT+(NB*STEP)):
            print ("EVOLUTION #",PROCESS)
            #now WA
            if WA_window=="moving" or WA_window=="downtop" or WA_window=="topdown":
                jump=(WAhigh-WAlow)/WA_window_slices
                #jumpBegin= WAlow #initialization
                #jumpEnd=(WAlow+jump)
                jumpBegin= WAlow #initialization
                jumpEnd=(WAlow+jump)
                #initial
                steps=-1
#                print (WAlow,WAhigh,jump,jumpBegin,"end",jumpEnd)
                while steps < WA_window_slices:
                    #print ("STEPS",steps,"    ",WAlow,WAhigh,jump,jumpBegin,"end",jumpEnd)
                    if steps==-1 :
                        jumpEnd=WAhigh
                        #data= arraysFiller(T,d18OWi,d18OAi,arraysize,jumpBegin,jumpEnd,analyticalProcessSimulation,sigmaLab)
                        if args.Temperature: data = arraysFiller(PROCESS,d18OWi,d18OAi,arraysize,jumpBegin,jumpEnd,analyticalProcessSimulation,sigmaLab)
                        if args.Sigma: data = arraysFiller(T,d18OWi,d18OAi,arraysize,jumpBegin,jumpEnd,analyticalProcessSimulation,PROCESS)
                        if args.Water : data= arraysFiller(T,PROCESS,d18OAi,arraysize,jumpBegin,jumpEnd,analyticalProcessSimulation,sigmaLab)
                        if args.Apatite : data=arraysFiller(T,d18OWi,PROCESS,arraysize,jumpBegin,jumpEnd,analyticalProcessSimulation,sigmaLab)
                        axes[linescpt][kolmncpt].hist(data, bins=bins, range=d19OAf_range,alpha=0.5, color='#A52A2A',density=False,label=str(jumpBegin)[:4]+"-"+str(jumpEnd)[:4],animated=True)
                        axes[linescpt][kolmncpt].set_xlabel(str(PROCESS), size=14)
                        leg = axes[linescpt][kolmncpt].legend(loc='upper left')
                        leg.draw_frame(False)
                        axes[linescpt][kolmncpt].set_xlabel("δ18OAf", size=12)
                        axes[linescpt][kolmncpt].set_ylabel("Count", size=12)
                        leg = axes[linescpt][kolmncpt].legend(loc='upper left')
                        axes[linescpt][kolmncpt].set_title("A/W ="+str(jumpBegin)[:4]+"-"+str(jumpEnd)[:4]+sentenceA+ str(PROCESS), fontsize=14)
                    else:
                        if args.Temperature: data = arraysFiller(PROCESS,d18OWi,d18OAi,arraysize,jumpBegin,jumpEnd,analyticalProcessSimulation,sigmaLab)
                        if args.Sigma: data = arraysFiller(T,d18OWi,d18OAi,arraysize,jumpBegin,jumpEnd,analyticalProcessSimulation,PROCESS)
                        if args.Water : data= arraysFiller(T,PROCESS,d18OAi,arraysize,jumpBegin,jumpEnd,analyticalProcessSimulation,sigmaLab)
                        if args.Apatite : data=arraysFiller(T,d18OWi,PROCESS,arraysize,jumpBegin,jumpEnd,analyticalProcessSimulation,sigmaLab)
                        axes[linescpt][1].hist(data, bins=bins, range=d19OAf_range,alpha=0.5, density=False, label=str(jumpBegin)[:4]+"-"+str(jumpEnd)[:4])
                        axes[linescpt][1].set_xlabel(str(PROCESS), size=14)
                        leg = axes[linescpt][1].legend(loc='upper left')
                        leg.draw_frame(False)
                        axes[linescpt][1].set_xlabel("δ18OAf", size=12)
                        axes[linescpt][1].set_ylabel("Count", size=12)
                        axes[linescpt][1].set_title("A/W varying (see legend) "+sentenceB+ str(PROCESS), fontsize=14)
                    if WA_window=="moving":
                        if steps==-1 :
                            jumpEnd=WAlow+jump
                            jumpBegin=WAlow
                        else:
                            jumpBegin=jumpEnd
                            jumpEnd=jumpEnd+jump
                    elif WA_window=="downtop":                    #plt.title('δ18OA final values by classes, variations of δ18OWi', fontsize=10)
                        if steps==-1 :
                            jumpEnd=WAhigh
                            jumpBegin=WAlow
                        else:
                            jumpBegin=jumpBegin+jump
                            jumpEnd=WAhigh
                    elif WA_window=="topdown":                    #plt.title('δ18OA final values by classes, variations of δ18OWi', fontsize=10)
                        if steps==-1 :
                            jumpEnd=WAhigh
                            jumpBegin=WAlow
                        else:
                            jumpBegin=WAlow #jumpBegin+jump
                            jumpEnd=jumpEnd-jump
                    else:
                        jumpEnd=WAhigh
                        jumpBegin=WAlow
                    steps+=1
            PROCESS+=STEP
            linescpt+=1
            #kolmncpt+=1
            steps=0
            #linescpt=0 if kolmncpt+=1 > ncols else kolmncpt+=1
            
        plt.suptitle("Effect of Variation of"+sentenceC+" in multiple A/W ratio ranges and for multiple "+sentenceD+" values\nC. Lecuyer and JP. Flandrois 2023\n"+userTitle, fontsize=16)
        #plt.title("Effect of Variation of δ18OWi in multiple A/W ratio ranges")
        fileIs=os.path.join(outputDir,"ApaOxIS_model"+"~"+jobdescription+"."+typeOfFile)
        #plt.show()
        plt.savefig(fileIs)
    
    
    if args.TwinTW or args.TwinAW: #only moving window on A/W range
        colors = ['#A52A2A','#800000','#8B4513', '#a0522d', '#D2691E', '#cd853f','#b8860b','#ffa07a', '#daa520' ,'#BC8F8F',  ]#["b", "r", "m", "w", "k", "g", "c", "y"]#ffdead #A52A2A #8B4513
        
        
        if WA_window_slices >4 : WA_window_slices=4
        jump=(WAhigh-WAlow)/(WA_window_slices)
        slicesWA=[[WAlow,WAhigh]]
        WABegin= WAlow #initialization
        WAEnd=(WAlow+jump)
        for j in range(0,WA_window_slices):
            jumpBegin= WAlow+j*jump #initialization
            jumpEnd=jumpBegin+jump
            slicesWA.append([jumpBegin,jumpEnd])
        print (jump,slicesWA,len(slicesWA))
        
        PROCESS1=INIT1
        PROCESS2=INIT2
        print (INIT1,INIT1+(NB1*STEP1))
        print (INIT2,INIT2+(NB2*STEP2))
        print ('---')
        
        jumpBegin= WAlow #initialization
        jumpEnd=WAhigh
        fig, axes = plt.subplots(nrows=NB1, ncols=len(slicesWA), figsize=figsize,constrained_layout=True)
        for w in range(0,len(slicesWA)):
            jumpBegin=slicesWA[w][0]
            jumpEnd=slicesWA[w][1]
            linescpt=0
            klr=0
            for u in range(NB1):
                klr=0
                for v in range(NB2):
                    
                    print (klr,colors[klr])
                    UN=INIT1+(u*STEP1)
                    DEUX=INIT2+(v*STEP2)
                    print (w,"PARAMETRES",jumpBegin,jumpEnd,UN,DEUX)
                    
#                    data= arraysFiller(UN,DEUX,d18OAi,arraysize,WAlow,WAhigh,analyticalProcessSimulation,sigmaLab)
#                    axes[linescpt][w].hist(data, bins=bins, range=d19OAf_range,alpha=0.2, density=False, label="δ18OWi "+str(DEUX)[:4])
                    
                    data= arraysFiller(UN,DEUX,d18OAi,arraysize,jumpBegin,jumpEnd,analyticalProcessSimulation,sigmaLab)
                    if  w==0:
                        axes[linescpt][w].hist(data, bins=bins, range=d19OAf_range,color = colors[klr],alpha=0.5, density=False, label="δ18OWi "+str(DEUX)[:4])
                    else:
                        axes[linescpt][w].hist(data, bins=bins, range=d19OAf_range,alpha=0.5, density=False, label="δ18OWi "+str(DEUX)[:4])
                    axes[linescpt][w].set_xlabel("δ18OAf", size=12)
                    axes[linescpt][w].set_title(VAR+str(UN)+"  A/W range:"+str(jumpBegin).split(".")[0]+"."+str(jumpBegin).split(".")[1][0:2]+"-"+str(jumpEnd).split(".")[0]+"."+str(jumpEnd).split(".")[1][0:2], size=14)
                    leg = axes[linescpt][w].legend(loc='upper left')
                    leg.draw_frame(False)
                    axes[linescpt][w].set_xlabel("δ18OAf", size=12)
                    axes[linescpt][w].set_ylabel("Count", size=12)
                    klr+=1
                linescpt+=1
        plt.suptitle("Effect of Variation of "+sentenceC+" within a range of A/W from "+(str(WAlow)+" to "+str(WAhigh))+ " \nC. Lecuyer and JP. Flandrois 2023\n"+userTitle, fontsize=16)
        fileIs=os.path.join(outputDir,"ApaOxIS_model"+"~"+jobdescription+"."+typeOfFile)
        #plt.show()
        plt.savefig(fileIs)
        
    
    
    
    
    
    
    
    EndinstructionAndBrowserPath=instructionAndBrowserPath.replace("FILEPATH",fileIs)
    os.system(EndinstructionAndBrowserPath)
                    
    plt.suptitle("Effect of Variation of"+sentenceC+" in multiple A/W ratio ranges and for multiple "+sentenceD+" values\nC. Lecuyer and JP. Flandrois 2023\n"+userTitle, fontsize=16)
        #plt.title("Effect of Variation of δ18OWi in multiple A/W ratio ranges")
    if typeOfFile =="pdf": fileIs=os.path.join(outputDir,"ApaOxIS_model"+"~"+jobdescription+".pdf")
    if typeOfFile =="svg":fileIs=os.path.join(outputDir,"ApaOxIS_model"+"~"+jobdescription+".svg")
        #plt.show()
    plt.savefig(fileIs)
    halt()
        
def arraysFiller(T,d18OWi,d18OAi,arraysize,WAlow,WAhigh,analyticalProcessSimulation,sigmaLab):
    """
    An array is constructed (WAarray) by sampling in an uniform repartition law with boundaries defined by WAlow,WAhigh
    Then a new array (d18OAf_array) is built by applying equation 4
    At the end the normality tests are done and the error ratio returned
    """
    print ("ENVOI",T,d18OWi,d18OAi,arraysize,WAlow,WAhigh,analyticalProcessSimulation,sigmaLab)
    DELTA=(117.4-T)/4.5 #equation 3 Lécuyer et al. (2013) - CL is for ChristopheLécuyer
    WAarray=np.random.uniform(WAlow,WAhigh, arraysize)
    #print ("            VERIF -> ",T,d18OWi,d18OAi)#"      ",T,d18OWi,d18OAi,sampleNB,alpha,arraysize,WAlow,WAhigh,Analytical_MC,sigmaLab)
    d18OAf_array=computeEquation(DELTA,d18OWi,d18OAi,WAarray)
    #print ("2",d18OAf_array)
    if analyticalProcessSimulation:
        d18OAf_array=np.random.normal(d18OAf_array,sigmaLab,arraysize)
        #print ("3",d18OAf_array)
    return d18OAf_array
    

def computeEquation(DELTA,d18OWi,d18OAi,WA):
    """
    This is the computation corresponding to equation (4) in our paper.
    It has been separated as a function to clarified the program.
    As the p value level comparison is not relevant, only the putative normality is returned as a boolean
    """
    Eq4result=(WA*(d18OWi+DELTA)+d18OAi)/(WA+1.0)
    return Eq4result

def readParameters(parameterYAML):
    """
    The 8 parameters are read from a YAML file O18parameters.yaml where they are stored as a python dictionaries
    And all the possibilities of MC computing are described
    This is by far the coolest way
    Reading is done at the beginning
    The structure of the dictionary parameterdictionary is shown in the function basicParameters()
    """
    #global parameterdictionary
    yaml_file = open(parameterYAML, 'r')
    parameterdictionary = yaml.load(yaml_file, Loader=yaml.SafeLoader) #this is the safe way tp prevent exploits
    try:
        modeling=parameterdictionary["Modelparameterdictionary"]["modeling"]
        print (modeling)
        arraysize=parameterdictionary["Modelparameterdictionary"]["arraysize"]
        print (arraysize)
        T=parameterdictionary["Modelparameterdictionary"]["T"]["initialvalue"]
        Tstep=parameterdictionary["Modelparameterdictionary"]["T"]["step"]
        Tnbsteps=parameterdictionary["Modelparameterdictionary"]["T"]["nb_steps"]
        print (T,Tstep,Tnbsteps)
        d18OWi=parameterdictionary["Modelparameterdictionary"]["d18OWi"]["initialvalue"]
        d18OWistep=parameterdictionary["Modelparameterdictionary"]["d18OWi"]["step"]
        d18OWinbsteps=parameterdictionary["Modelparameterdictionary"]["d18OWi"]["nb_steps"]
        print (d18OWi,d18OWistep,d18OWinbsteps)
        d18OAi=parameterdictionary["Modelparameterdictionary"]["d18OAi"]["initialvalue"]
        d18OAistep=parameterdictionary["Modelparameterdictionary"]["d18OAi"]["step"]
        d18OAinbsteps=parameterdictionary["Modelparameterdictionary"]["d18OAi"]["nb_steps"]
        print (d18OAi,d18OAistep,d18OAinbsteps)
        WAlow=parameterdictionary["Modelparameterdictionary"]["WA"]["WAlow"]
        WAhigh=parameterdictionary["Modelparameterdictionary"]["WA"]["WAhigh"]
        WA_window=parameterdictionary["Modelparameterdictionary"]["WA"]["WA_window"] #"moving" "topdown" "downtop" "fixed"
        WA_window_slices=parameterdictionary["Modelparameterdictionary"]["WA"]["WA_window_slices"]
        print (WAlow,WAhigh,WA_window,WA_window_slices)
        analyticalProcessSimulation=parameterdictionary["AnalyticalProcessSimulation"] #True or False
        sigmaLab=parameterdictionary["analyticalStd"]["sigma"]
        sigmaLabstep=parameterdictionary["analyticalStd"]["step"]
        sigmaLabnbsteps=parameterdictionary["analyticalStd"]["nb_steps"]
        print (analyticalProcessSimulation,sigmaLab,sigmaLabstep,sigmaLabnbsteps)
        threadsCount=parameterdictionary["threadsCount"] #auto or a number
        instructionAndBrowserPath=parameterdictionary["instructionAndBrowserPath"]
        usersize=parameterdictionary["Rendering"]["figSize"].split(',')
        print(threadsCount,instructionAndBrowserPath,usersize)
        figsize=(int(usersize[0]),int(usersize[1]))
        print(figsize)
        userTitle=parameterdictionary["Rendering"]["userTitle"]
        bins=parameterdictionary["Rendering"]["bins"]
        print (userTitle)
        print (bins)
        d19OAf_range=(float(parameterdictionary["Rendering"]["d18OAf_range"].split(',')[0]),float(parameterdictionary["Rendering"]["d18OAf_range"].split(',')[1]))
        typeOfFile=parameterdictionary["Rendering"]["outputAs"]
        print (typeOfFile)
        print ("All is OK")
        print (parameterdictionary)
        
        return arraysize,T,Tstep,Tnbsteps,d18OWi,d18OWistep,d18OWinbsteps,d18OAi,d18OAistep,d18OAinbsteps,WAlow,WAhigh,WA_window,WA_window_slices,analyticalProcessSimulation,sigmaLab,sigmaLabstep,sigmaLabnbsteps,threadsCount,instructionAndBrowserPath,figsize,bins,userTitle,d19OAf_range,typeOfFile
    except :
        stop()
def stop(valeur='\n----------------------------------\nMODELING IS CANCELED\n\nSEE YOU SOON !\n----------------------------------\n'):
    print('----------------------------------\nPARAMETERS FILE IS NOT WELL FITTED\n----------------------------------\n'+str(valeur))
    exit(1)
def halt(valeur='\n----------------------------------\nMODELING IS FINISHED\n\nSEE YOU SOON !\n----------------------------------\n'):
    exit(1)
def createFolder(directory):
    """
    From BIBIengine and TEALS (jp Flandrois 2020), extremely basic and I suppose that it can be found elsewhere
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
        #else:
            #os.system('rm -r '+directory)
            #os.makedirs(directory)
    except OSError:
        print(('Error: Creating directory. ' +  directory))
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="ApaOxIS",description='''
    APAtite OXygen Isotope Simulation (modeling).
    ApaOxIS has been designed to simulate the fate of biogenic apatite δ18O during diagenesis in given conditions. Note that the .yaml file contains all the parameters, including parameters for the rendering and visiualisation.
    ''',
    epilog="""
    place our references here (paper, site...)
    """)

    #print("\n\n#####################################  ApaOxIS ###################################################\n")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-W", "--Water", action='store_true', help="Cross variations of δ18OWi and A/W")
    group.add_argument("-A", "--Apatite", action='store_true', help="Cross variations of δ18OAi and A/W")
    group.add_argument("-T", "--Temperature", action='store_true', help="Cross variations of T and A/W")
    group.add_argument("-S", "--Sigma", action='store_true', help="Cross variations of analytical σ and A/W")
    group.add_argument("-m", "--Model", action='store_true', help="Only one histogram for given fixed conditions")
    group.add_argument("-TW", "--TwinTW", action='store_true', help="Both variations of T and δ18OWi within 1-4 ranges of A/W : a collection of histograms")
    group.add_argument("-AW", "--TwinAW", action='store_true', help="Both variations of δ18OAi and δ18OWi within 1-4 ranges of A/W : a collection of histograms")
    parser.add_argument("i", type=str, help="the parameter file (yaml file, WARNING: structure is mandatory)")
    parser.add_argument("o", type=str, help="The _directory_ that will gather the results ")
    args = parser.parse_args()
    #print (args)
    userTitle=""
    createFolder(args.o)
    
    
    if args.Model:
        print("\n\n#####################################  ApaOxIS ###################################################\n")
        arraysize,T,Tstep,Tnbsteps,d18OWi,d18OWistep,d18OWinbsteps,d18OAi,d18OAistep,d18OAinbsteps,WAlow,WAhigh,WA_window,WA_window_slices,analyticalProcessSimulation,sigmaLab,sigmaLabstep,sigmaLabnbsteps,threadsCount,instructionAndBrowserPath,figsize,bins,userTitle,d19OAf_range,typeOfFile=readParameters(args.i)
        jobdescription=""
        for u in ("m",arraysize,T,Tstep,Tnbsteps,d18OWi,d18OWistep,d18OWinbsteps,d18OAi,d18OAistep,d18OAinbsteps,WAlow,WAhigh,WA_window,WA_window_slices,analyticalProcessSimulation,sigmaLab):
            jobdescription+=str(u)+'~'
        
#        DELTA=(117.4-T)/4.5
#        print ("T",T,"A",computeEquation(DELTA,-8,20,0.5))
#        halt(computeEquation(DELTA,-8,20,0.5))
        coreProgram(args.o,arraysize,T,Tstep,Tnbsteps,d18OWi,d18OWistep,d18OWinbsteps,d18OAi,d18OAistep,d18OAinbsteps,WAlow,WAhigh,WA_window,WA_window_slices,analyticalProcessSimulation,sigmaLab,sigmaLabstep,sigmaLabnbsteps,threadsCount,instructionAndBrowserPath,figsize,bins,userTitle,d19OAf_range,typeOfFile,jobdescription.strip('~'))
    if args.Water or args.Temperature or args.Sigma or args.Apatite or args.TwinTW or args.TwinAW:
        print("\n\n#####################################  ApaOxIS ###################################################\n")
        arraysize,T,Tstep,Tnbsteps,d18OWi,d18OWistep,d18OWinbsteps,d18OAi,d18OAistep,d18OAinbsteps,WAlow,WAhigh,WA_window,WA_window_slices,analyticalProcessSimulation,sigmaLab,sigmaLabstep,sigmaLabnbsteps,threadsCount,instructionAndBrowserPath,figsize,bins,userTitle,d19OAf_range,typeOfFile=readParameters(args.i)
        jobdescription=""
        for u in (sys.argv[1],arraysize,T,Tstep,Tnbsteps,d18OWi,d18OWistep,d18OWinbsteps,d18OAi,d18OAistep,d18OAinbsteps,WAlow,WAhigh,WA_window,WA_window_slices,analyticalProcessSimulation,sigmaLab):
            jobdescription+=str(u)+'~'
        coreProgram(args.o,arraysize,T,Tstep,Tnbsteps,d18OWi,d18OWistep,d18OWinbsteps,d18OAi,d18OAistep,d18OAinbsteps,WAlow,WAhigh,WA_window,WA_window_slices,analyticalProcessSimulation,sigmaLab,sigmaLabstep,sigmaLabnbsteps,threadsCount,instructionAndBrowserPath,figsize,bins,userTitle,d19OAf_range,typeOfFile,jobdescription.strip('~'))
