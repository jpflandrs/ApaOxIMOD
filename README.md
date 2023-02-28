# ApaOxIMOD
 
# Apatite Oxygen Isotopes Modeling

## Aims

ApaOxMod simulate the evolution of biological apatite δ18O (essentially in thanatocenoses) submitted to a diagenesis process in various environments after its creation.

## Version 
The current version is mk03 2023  

    Note: that the "testing" version open the possibility to draw only a part of the output to simplify the graphs for a publication. This version is not usefull and may be misleading during the exploration the phenomenon as some output are masked.  

## Authors 

### JP Flandrois 
jean-pierre.Flandrois@univ-lyon1.fr

Author and maintainer
### C Lecuyer
christophe.lecuyer@univ-lyon1.fr

Original idea and design, equations

## User Licence
The program is made available under the [CeCILL2.1](http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.txt) licence.

## Rapid description

### Hypothesis

In condition of sedimentary deposits in short time and limited space the distribution of biological apatite δ18O data is expected to be _issued from a population following a Normal Law_. 

If the samples have been submitted to a diagenesis process the original apatite δ18O signal is blurred and the resulting distribution is certainly non-normal.


### Model parameters
* It is possible to compute (under some hypothesis defined in the companion "yaml" file) one or a family of histograms depicting the result of alteration.
* The array size (not less than 250.000) corresponds to the number of simulated simulations, ideally 1.000.000.
* The initial values of the parameters are set: δ18OAi (apatite), δ18OWi (water), T (temperature) and W/A (apatite/water ratio)
* The δ18OAf (final apatite δ18O) is computed for each W/A condition according to the values of the other parameters that are allowed to vary.
    * The conditions of the diagenesis process are given: 
        * by describing a window (WAlow, WAhigh)  and a number of slices and the way to explore the W/A window 

                This is the window : [         ]
                "fixed" the window is identical to the W/A space [[xxxxxxx]] 
                "topdown" schreading from top [[xxxxxxx]]-> [[xxxxxx] ]-> [[xxxxx]  ]-> [[xxxx]   ]-> [[xxx]    ]
                "downtop" schreading from bottom [[xxxxxxx]]-> [ [xxxxxx]]-> [  [xxxxx]]-> [   [xxxx]]-> [    [xxx]]
                "moving" [[xxx]    ]-> [ [xxx]   ]-> [  [xxx]  ]-> [   [xxx] ]-> [    [xxx]]
                Note that the full range is always studied and shown in the first column of the histograms.
    * by giving an interval and a number of steps of increase (T, δ18OAi, δ18OWi).         
* There is also the possibility to introduce the analytical process uncertainty.
* Do not ask too much :

        As 4 parameters may be set and explored by steps (T, d18OWi, d18OAi) 
        and W/A may be constant or varying by using a range and a sub-range 
        as the way to explore the W/A effect...
        Note that 10 values for T, d18OWi and d18OAi led to 9x9x9 and 
        will generate 81 pages containing 9 histograms !
        Note that adding W/A variations may increase this to more than 1000 pages !

### Usage

``ApaOxIMOD [-h] [-W | -A | -T | -S | -m | -TW | -WA] i o``

#### Positional arguments:
*  i                  Parameter file (yaml file, WARNING: structure is mandatory) **or** user data
*  o                  The _directory_ that will gather the results

#### Optional arguments:
*  -h, --help        : Show an help message and exit
*  -m, --Model       : Only one histogram for given fixed conditions
*  -W, --Water       : Cross variations of δ18OWi and W/A
*  -A, --Apatite     : Cross variations of δ18OAi and W/A
*  -T, --Temperature : Cross variations of T and W/A
*  -S, --Sigma       : Cross variations of analytical σ and W/A
*  -TW, --TwinTW     : Both variations of T and δ18OWi within 1-4 ranges of W/A : a collection of histograms
*  -WA, --TwinWA     : Both variations of δ18OAi and δ18OWi within 1-4 ranges of W/A : a collection of histograms

#### Other explanations 
There are also explanations in the MODELO18parameters.yaml file.
To rebuilt the outputs provided in the publication, see the instructions in the specific parameter files MODELO18testing{1,2,3}.yaml.

### Outputs

The outputs are histograms organised by rows and lines. On the left row the histograms are that of the whole W/A window, on the left row each of the subwindows are discriminated by colors. 

Each histogram line descript the effect of given range of the other variable parameter (if -TW or -WA of the two variable parameters)

Setting the histogram parameters is possible in the yaml file. Output is saved and also displayed to the browser if defined. Note that all the values of the parameters are given in the file name in the order given in the yaml file:

``ApaOxIS_model~m~500000~15~5~5~-8.0~-2~5~20.0~1~9~0.05~0.95~moving~2~True~0.1.svg``

but you may also change the name of the yams file and save it with the outputs as a quality sheet.

### Conclusion

ApaOxIMOD enables: 

* The analysis of the diagenesis phenomenon of apatites is various physical and chemical conditions.  

## How it works :

Monte Carlo (MC) simulation is based on the creation of an array of W/A values taken from a uniform distribution between given threshold values of δ18O. This represents the diversity of the deposit situation and physical state of the biological apatite.  Typically 1,000,000 sets of length n are simulated. 

The second phase uses equation (4) to compute the δ18O in biological apatite δ18OAf at the equilibrium on each item of the whole array given (range of values) the temperature T, initial δ18O in biological apatite (δ18OAi) and the δ18O in water δ18OWi. 

The resulting δ18OAf array is analysed by a set of histograms.   
                         
