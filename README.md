#### Joint Calorimeter Prototype Analysis --- JCaPA           
#### Created by Y.Kulinich, R.Longo and C.Lantz on 12/12/2018 ####                                                                                                    
                                                                                                                                                 
Basic structure defined and discussed during the Thursday meeting on 12/13/2018 

JCaPA
#	Analysis
#		include 
#		src
#		userFunctions
#	MC (to be implemented) 
#	2018_Utils (to be added) 

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
A doxygen documentation can also be created following the instruction below 

#### Doxygen documentation ####
First, check that doxygen is installed on your machine. 
If it's not the case, just check it out using 
By default a folder "doxygen" will also be installed in your JCaPA_INSTALL folder. 


 

