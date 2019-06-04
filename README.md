#### Joint Zerodegree Calorimeter Prototype Analysis --- JZCaPA           
#### Created by Y.Kulinich, R.Longo and C.Lantz on 12/12/2018 #### 
#### Migrated from CERN to UIUC Gitlab on 14th February 2019 ####
                                                                                                                                                 
Basic structure defined and discussed during the Thursday meeting on 12/13/2018 

JZCaPA   
     Analysis   
          include    
          src   
          userFunctions   
     MonteCarlo    
     2018_Utils 

The project is cmake based, so you need a reasonably new cmake version ( version > 2.8 )   
The standalone Analysis part requires only a root installation (https://root.cern.ch) and a xerces-c installation (http://xerces.apache.org). 
Please note that xerces-c is usually available via your package installer (so easy to get installed). 
The MC part will be conditional since it requires additional software as Geant4 and all its dependencies   
The corresponding README part will be written once MC will be included.    
The 2018_Utils folder contains the .xml settings file encoding the useful information on 2018 beam test (automatically loaded in JZCaPA).   

#### CMake and installation ####
To install the software using cmake will be trivial. 
In the same folder where you have JZCaPA, just do
```bash
mkdir JZCaPA_BUILD

mkdir JZCaPA_INSTALL 
```
at this stage, remember to add to your environment
```bash
export JZCaPA=/path/to/your/JZCaPA_INSTALL

cd JZCaPA_BUILD

cmake -DCMAKE_INSTALL_PREFIX=../JZCaPA_INSTALL/ ../JZCaPA

make -j8 

make install 
```
please remember to re-make & make install every time you change the source code 

#### Analysis ####
Each user can implement his/her own analysis creating a new userFunction.cpp in Analysis/userFunctions folder. 
Please check AnalysisExample.cpp if you're looking for a basic template. 

The code has been updated and improved considerably in the last period, so further description will be provided afterwards. 

They are well commented by Yakov for each available method.     
A doxygen documentation can also be created following the instruction below.     

#### Monte Carlo ####
Monte Carlo folder added on 12/19/2018.    
What's there at the moment is exactly what has been done so far by Mike Phipps (michael.william.phipps@cern.ch). 
The only modifications implemented were to ensure the compatibility with newer Geant4 versions (>= 10.4.3).   
Such changes have been implemented in order to allow also for backward compatibility with older Geant4 versions. 

Please note that MonteCarlo support is *DISABLED* by default. This choice is meant to avoid people not interested in MC to install Geant4 and the other dependencies. 
In order to enable it, add the option 
```bash
-DJZCaPA_ENABLE_MC=YES 
```
while cmaking. Please note that you need the Geant4 toolkit (and the corresponding dependencies) to successfully enable the MC support. 
More details will come in the future. 

#### Doxygen documentation ####
First, check that doxygen is installed on your machine.   
If it's not the case, just check it out using   
By default a folder "doxygen" will also be installed in your JZCaPA_INSTALL folder.   
To obtain the documentation, just execute   
```bash
doxygen JZCaPA_doxy.cnf   
```
This will generate for you  $JZCaPA/html and $JZCaPA/latex folder.  
If you start $JCaPA/html/index.html with your browser, you will get your doxygen docs (surfable).   
Alternatively, you can compile the latex static documentation in $JZCaPA/latex 


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
 
#### XML Setting file support #### 
Since 29th January you will need to install also xercesc as a dependency to compile JZCaPA. This is rather straight forward on machines where you have root access. 
You've just to run 
```bash
your_package_installer install xerces-c 
```
If you are on a cluster and xerces-c is not installed/available by default, is straightforward to checkout the source from http://xerces.apache.org/xerces-c/download.cgi and compile it following the instructions available on the webpage. 
If you are going for the source installation, please note that, for an easier detection by cmake, it is adviced to define the variables $XERCESC_DIR, $XERCESC_INCLUDE_DIR and $XERCESC_LIBRARY
Thanks to the .xml settings file support (hopefully) all the 2018 test beam settings will be loaded automatically when starting the analysis and associated to detectors objects (ZDC1, ZDC2 or RPD). 
More information can be added if needed, just formulate a request or do the implementation yourself!  

#### Install JZCaPA on LXPLUS #### 
Given the forthcoming migration to centos 7 of the whole lxplus system (2nd April 2019), the configuration was established directly for this OS. 
The alias of lxplus will switch automatically to centos 7 from 2nd April 2019. Until this time, to profit of this configuration do 

```
ssh your_username@lxplus7.cern.ch 
```

If you are using bash, after logging in, add ${JZCaPA}/Utils/jcapa_env_lxplus.sh to your .bashrc
Please be sure to comment out the source of other softwares/compilers. This may conflict with the view that is being loaded by the script. 
After this modification, log-out and in again. You will have all the dependencies needed (ROOT >= 6.08, Geant4 and xerces-c) loaded in your environment.
Go to the JZCaPA_BUILD directory and cmake using the flag 

```
-DJZCaPA_ON_LXPLUS=YES
```

and compile the software (with or w/o MC support, depending on your needings). The software should now compile fully.

#### Install JZCaPA on UIUC #### 
The environment for UIUC is built in Centos 7 and uses an HTCondor cluster, similar to what is available on lxplus. However, a VPN is needed
to access these resources. Instructions for download and setup can be found at https://techservices.illinois.edu/services/virtual-private-networking-vpn/download-and-set-up-the-vpn-client
When you are connected to the VPN, you can access machines by

```
ssh your_username@jzcapa-vm#.physics.illinois.edu (where # is 1-6)
or
ssh your_username@htc-login.campuscluster.illinois.edu
```

The first time you log in, enter

```
source /path/to/JZCaPA/source/Utils/jzcapa_env_uiuc.sh

echo 'source ${JZCaPA}/Utils/jzcapa_env_uiuc.sh' >> ~/.bashrc
```

The first line sources the required software so you can install JZCaPA, and the second line sources it for subsequent use.
You will have all the dependencies needed (ROOT 6.16, Geant4 and xerces-c) loaded in your environment.
Now you can follow the install instructions from the top and compile the software (with or w/o MC support, depending on your needs). 
The software should now compile fully.
