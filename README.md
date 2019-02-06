#### Joint Calorimeter Prototype Analysis --- JCaPA           
#### Created by Y.Kulinich, R.Longo and C.Lantz on 12/12/2018 ####                                                                                                    
                                                                                                                                                 
Basic structure defined and discussed during the Thursday meeting on 12/13/2018 

JCaPA   
     Analysis   
          include    
          src   
          userFunctions   
     MonteCarlo    
     2018_Utils (empty for the moment)    

The project is cmake based, so you need a reasonably new cmake version ( version > 2.8 )   
The standalone Analysis part requires only a root installation (https://root.cern.ch) and a xerces-c installation (http://xerces.apache.org). 
Please note that xerces-c is usually available via your package installer (so easy to get installed). 
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


#### Structure of processed data ####
TTree structure of the processed data is illustrated in the following.

Tree name: "tree".
Inside one can find a variety of branches.
In 2018 we had N = 20 channels read. 5 DRS4 Boards x 4 Channels each. Each channel is configured to output 1024 samples.

Raw data is stored as a vector<float> in branches called 'RawCn' and 'Cn', where n is an element { 0, ..., N-1 }.
Each branch RawCn and Cn will be of size 1024 (corresponding to the number of samples).
Therefore for each event we have N vectors of size 1024.
RawCn corresponds to the "raw" output, what one would see on the DRS scope.
Cn corresponds to processed output, as determined by rdcaqAnalysis.
One can take the RawCn and reprocess it as they wish using JCaPA analysis methods (or implementing new ones). 

Also, in the tree, there are branches that come from rcdaqAnalysis which have information on the output of each channel.
These are branches such as MaxCharge. Lets say there are M of them.
These are vectors of size N, you will have one piece of information for each channel. Be it max charge, or pedestal rms, etc.
Here you have information from rcdaqAnalysis depending on how it was configured. The user can use this and or do their own processing from the raw data.
Depending on how many pieces of processed info are saved from rcdaqAnalysis, you will have M vectors of size N.

Summarizing, there are
N branches of size 1024 for the channel waveforms.
and
M branches of size N for the processed data from rcdaqAnalysis. M is an arbitrary number.
 
#### XML Setting file support - Beta version #### 
Since 29th January you will need to install also xercesc as a dependency to compile JCaPA. This is rather straight forward on machines where you have root access. 
You've just to run 

your_package_installer install xerces-c 

If you are on a cluster and xerces-c is not installed/available by default, is straightforward to checkout the source from http://xerces.apache.org/xerces-c/download.cgi and compile it following the instructions available on the webpage. 
If you are going for the source installation, please note that, for an easier detection by cmake, it is adviced to define the variables $XERCESC_DIR, $XERCESC_INCLUDE_DIR and $XERCESC_LIBRARY
