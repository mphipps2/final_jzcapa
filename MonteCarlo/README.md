 Instructions from Mike Phipps (taken from his repo) 

 This simulation is intended to allow standalone (non-ATHENA) simulation of the current ZDC geometry, as well as configurable geometries that could be used for prototype studies.
	
 1- GEOMETRY DEFINITION
	
   The world is defined and the individual modules created in DetectorConstruction.cc. Depending what type of modules you've specified in the config file, different module construction classes will be called. For the EM Pixel module, use ModType1; for the Hadronic Pixel module, use ModType2; for the Hadronic Non-Pixel modules, use ModType3; and for the generic configurable modules (eg prototype designs), use ModTypeCustom. Parameters within mod types 1-3 should be kept static (in general -- since these geometries reflect the current ZDC geometry and that used in the ATHENA simulation). Parameters for custom modules should be chosen carefully within the config file.   

   For Modules 1 and 2, the geometry is highly "pixelized". In other words, it is composed of many geometry segments that must be stitched together. For example, each layer of tungsten plate is composed of 9 "lateral plates" and 8 "longitudinal plates" that are stitched together. The lateral plates are in the plane of the vertical strips and have a depth equal to the absorber gap depth. The longitudinal plates, on the other hand, run the length of the module and are located in the gaps between vertical strip groups. This scheme is done to match the geometries found in these modules in which pixels run the length of the detector in the gaps between the vertical plates, with air gaps in the vertical space between the pixels in each radiator gap. This may seem excessive, but within Geant, this sort of stitched-together, pixelized scheme is preferable to one with many boolean cut-outs. The processing time for this scheme should be similar to the other modules and an order of magnitude faster than the boolean cut-out option. For a more detailed explanation of this scheme and a diagram of the inheritance, see the ZdcSimulation section of the ATLAS ZDC twiki: https://twiki.cern.ch/twiki/bin/view/Atlas/ZdcSimulation

   For Modules 1 and 2, the pixel rods are currently not used as active elements of the detector. In other words, hits within these volumes are not retained.

   Module 4 allows the following configurations:
   - As an absorber block that could be used for test beams. To insert this, turn off the cladding, set NStripsPerGap, NRadiators, RadiatorGapLength and CladdingThickness to 0, and then set the Absorber size appropriately. Note: do not change the CoreDiameter to 0
   - As a prototype module with cladding (but no buffer)
   - As a prototype module with different sampling ratios
   - As a prototype module that uses lead or tungsten absorber
		
 2- PHYSICS LIST
 
   The standard ATLAS hadronic physics lists is FTFP_Bert (Fritiof string model (>~5 GeV); Bertini-style cascade (<~10 GeV)) and that is the default list in this simulation. Although this can be changed in the config file. The other hadronic list currently configured is QGSP_Bert.

   By default, FTFP_Bert comes with the following EM and low energy lists also turned on: G4EmStandard, G4GAmmaLeptoNuclearPhys, Decay, hElasticWEL_CHIPS, hInelastic FTFP_BERT, stopping, ionInelasticFTFP_BIC, neutronTrackingCut.

   Pion and muon decay have also been turned on within the PhysicsList.cc module
   
 3- ACTION INITALIZATION

   Action initialization begins within ActionInitialization.cc, where RunAction, EventAction and SteppingAction are all declared.

   The RunAction module is just used to initialize the output TTree at the start of each run.

   The EventAction module is used primarily at the end of each event to fill the TTree.

   The SteppingAction module is currently just a skeleton.
  	 
 4- PRIMARY GENERATOR
  
   The PrimaryGeneratorAction module uses G4ParticleGun to initialize the default positioning and configuration of the particle generator. By default the gun is set at the edge of the first module, centered in the (x,y) transverse plane.
     
 5- DETECTOR RESPONSE

   The detector response is handled in the QuartzSD module. Each quartz rod (and notice, these are the vertical quartz rods, not the pixel rods) is set as a sensitive detector and hits within these volumes are parsed in the Quartz Sensitive Detector class (QuartzSD.cc) and QuartzHits are created and filled in Collection vectors. This allows retention of detailed information about the hit including the particular strip number it was registered in, the module, the radiator gap number, the trackID, the particleID, as well as kinematic information. This is also the stage where individual hits are turned into Cherenkov radiation and the number of Cherenkov photons retained in the acceptance of the quartz rod is calculated. This calculation finds the average number of photons created by each charged particle above threshold, shifts it by a Poisson smearing and then uses the appropriate emittance angle to determine whether the photon would have been captured by the quartz. All photons directed vertically downward in the detector are assumed to be lost, while all photons captured and directed upward are considered to be retained. This calculation also assumes Cherenkov photons are created in the quartz core and not the cladding.
   
   The Cherenkov sampling window can be set from the config file. The min and max wavelengths are taken as inputs. The Frank-Tamm relationship is used to calculate the number of photons created across the step length and the wavelength range specified. Note, this relationship is sensitive to the charge^2, 1/(beta^2), and 1/(n^2), with the number of photons falling off linearly with wavelength.
   

 A- VISUALISATION

   The visualization manager is set via the G4VisExecutive class
   in the main() function in exampleB1.cc.    
   The initialisation of the drawing is done via a set of /vis/ commands
   in the macro vis.mac. This macro is automatically read from
   the main function when the program is used in interactive running mode.

   By default, vis.mac opens an OpenGL viewer (/vis/open OGL).
   The user can change the initial viewer by commenting out this line
   and instead uncommenting one of the other /vis/open statements, such as
   HepRepFile or DAWNFILE (which produce files that can be viewed with the
   HepRApp and DAWN viewers, respectively).  Note that one can always
   open new viewers at any time from the command line.  For example, if
   you already have a view in, say, an OpenGL window with a name
   "viewer-0", then
      /vis/open DAWNFILE
   then to get the same view
      /vis/viewer/copyView viewer-0
   or to get the same view *plus* scene-modifications
      /vis/viewer/set/all viewer-0
   then to see the result
      /vis/viewer/flush

   The DAWNFILE, HepRepFile drivers are always available
   (since they require no external libraries), but the OGL driver requires
   that the Geant4 libraries have been built with the OpenGL option.

   From Release 9.6 the vis.mac macro has additional commands
   that demonstrate additional functionality of the vis system, such as
   displaying text, axes, scales, date, logo and shows how to change
   viewpoint and style.  Consider copying these to other examples or
   your application.  To see even more commands use help or
   ls or browse the available UI commands in the Application
   Developers Guide, Section 7.1.

   For more information on visualization, including information on how to
   install and run DAWN, OpenGL and HepRApp, see the visualization tutorials,
   for example,
   http://geant4.slac.stanford.edu/Presentations/vis/G4[VIS]Tutorial/G4[VIS]Tutorial.html
   (where [VIS] can be replaced by DAWN, OpenGL and HepRApp)

   The tracks are automatically drawn at the end of each event, accumulated
   for all events and erased at the beginning of the next run.

 B- USER INTERFACES
 
   The user command interface is set via the G4UIExecutive class
   in the main() function in zdc.cc. This script is the driver for your simulation 
    
 
 C- HOW TO RUN

    - Execute zdc in the 'interactive mode' with visualization (must be run from the atlasZDC directory):
        % ./build/zdc
      and type in the commands from run1.mac line by line:  
        Idle> /run/beamOn 1 
        Idle> ...
        Idle> exit
      or
        Idle> /control/execute run1.mac
        ....
        Idle> exit

    - Execute zdc in the 'batch' mode from macro files 
      (without visualization)
        % ./build/zdc run.mac results/resultsFile.root

    - Execute zdc MANY TIMES in the 'batch' mode from macro files 
      (without visualization) -- This is intended if you events broken up into separate runs/ttrees -- if you use this you should modify the arguments inside scan.sh 
        % source scan.sh
       
