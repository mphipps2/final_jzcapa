#### Joint Calorimeter Prototype Analysis --- JCaPA           
#### Created by Y.Kulinich, R.Longo and C.Lantz on 12/12/2018 ####                                                                                                    
                                                                                                                                                 
Basic structure defined and discussed during the Thursday meeting on 12/13/2018 

JCaPA   
     Analysis   
          include    
          src   
          userFunctions   
     MC (to be implemented)    
     2018_Utils (empty for the moment)    

The project is cmake based, so you need a reasonably new cmake version ( version > 2.8 )   
The standalone Analysis part requires only a root installation (https://root.cern.ch)   
The MC part will be conditional since it requires additional software as Geant4 and all its dependencies   
The corresponding README part will be written once MC will be included.    
The 2018_Utils folder will be meant to contain useful files for 2018 test beam analysis (e.g. Summary of various scans etc)   

#### CMake and installation ####
To install the software using cmake will be trivial. 
In the same folder where you have JCaPA, just do

mkdir JCaPA_BUILD

mkdir JCaPA_INSTALL 

at this stage, remember to add to your environment

export JCaPA=/path/to/your/JCaPA_INSTALL

cd JCaPA_BUILD

cmake -DCMAKE_INSTALL_PREFIX=../JCaPA_INSTALL/ ../JCaPA

make -j8 

make install 

please remember to re-make & make install every time you change the source code 

#### Analysis ####
Each user can implement his/her own analysis creating a new userFunction.cpp in Analysis/userFunctions folder. 
Please check AnalysisExample.cpp if you're looking for a basic template. 

Two main classes are provided at the moment: 
- DataReader 
- WFAnalysis (inherits from Analysis.h)
They are well commented by Yakov for each available method.   
A doxygen documentation can also be created following the instruction below.     

#### Monte Carlo ####
Monte Carlo folder added on 12/19/2018.    
What's there at the moment is exactly what has been done so far by Mike Phipps (michael.william.phipps@cern.ch). 
The only modifications implemented were to ensure the compatibility with newer Geant4 versions (>= 10.4.3).   
Such changes have been implemented in order to allow also for backward compatibility with older Geant4 versions. 

Please note that MonteCarlo support is *DISABLED* by default. This choice is meant to avoid people not interested in MC to install Geant4 and the other dependencies. 
In order to enable it, add the option 

-DJCaPA_ENABLE_MC=YES 

while cmaking. Please note that you need the Geant4 toolkit (and the corresponding dependencies) to successfully enable the MC support. 
More details will come in the future. 

#### Doxygen documentation ####
First, check that doxygen is installed on your machine.   
If it's not the case, just check it out using   
By default a folder "doxygen" will also be installed in your JCaPA_INSTALL folder.   
To obtain the documentation, just execute   

doxygen JCaPA_doxy.cnf   

This will generate for you  $JCaPA/html and $JCaPA/latex folder.  
If you start $JCaPA/html/index.html with your browser, you will get your doxygen docs (surfable).   
Alternatively, you can compile the latex static documentation in $JCaPA/latex 


 

