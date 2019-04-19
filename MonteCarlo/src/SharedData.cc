/** @file SharedData.cxx
 *  @brief Implementation of SharedData.
 *
 *  SharedData class holds data common to analysis packages.
 *  This includes histograms, trees, event counter.
 *  Analysis packages will have a pointer to this class passed from
 *  the manager.
 * 
 *  @author Yakov Kulinich
 *  @bug No known bugs.
 */
#include "SharedData.hh"
#include <iostream>

/** @brief Default Constructor for SharedData.
 */
SharedData :: SharedData ()
{
  m_eventCounter = 0;     
  m_outputFileName = "myOut.root";
  m_configFileName = std::getenv("JZCaPA");
  m_configFileName.replace(m_configFileName.length()-15,15,"/JZCaPA/MonteCarlo/config/config.cfg");
 // std::cout << "* current config file path * = " << m_configFileName << std::endl;
 // m_configFileName = "config/config.cfg";

  m_fout = NULL;
  m_tree = NULL;
  m_config = NULL;

}
/** @brief Constructor for SharedData. 
 *
 *  @param1 string with output filename
 *  @param2 string with name of config file
 */
SharedData :: SharedData ( const std::string& outputFileName,
			   const std::string& configFileName )
{
  m_eventCounter = 0;
  m_outputFileName = outputFileName;
  m_configFileName = configFileName;
 // std::cout << "* current config file path * = " << m_configFileName << std::endl;
 
  m_fout = NULL ;
  m_tree = NULL;
  m_config = NULL;
}

/** @brief Destructor for SharedData.
 *
 *  Cleans up an SharedData object.
 */
SharedData :: ~SharedData() 
{
  delete m_tree;
  delete m_fout;
  delete m_config;
}


/** @brief Function to add an event store.
 *
 *  Initialize TFile, Tree, TEnv (config) 
 *
 *  @return void
 */

void SharedData :: Initialize()
{
  m_fout         = new TFile( m_outputFileName.c_str(), "RECREATE" );
  m_tree         = new TTree( "tree"                  , "tree"     );

  m_config       = new TEnv ();
  int success;
  success = m_config->ReadFile( m_configFileName.c_str(), EEnvLevel(0));
  std::cout << " Config File path = " <<  m_configFileName.c_str() << " Check if read successfully (0 = success) " << success << std::endl;
}

/** @brief Function to add an output histo.
 *
 *  @param1 Pointer to histogram
 *
 *  @return void
 */
void SharedData :: AddOutputHistogram( TH1* h )
{
  m_v_hists.push_back( h );
}
/** @brief End of event
 *
 *  Fill the tree, increment event counter 
 *
 *  @return void 
 */
void SharedData :: EndOfEvent()
{
  //std::cout << " filling tree " << std::endl;
  //  m_tree->Fill();
  //    std::cout << " filling tree 2" << std::endl;
  m_eventCounter++;
}

/** @brief Function to check if printing necessary
 *
 *  Print every event until 10, then every 10,
 *  then every 100, then every 1000, etc.
 *
 *  @return bool to print or not
 */
bool SharedData :: DoPrint() 
{
  int statSize=1;
  if( m_eventCounter != 0){
    double power=std::floor(log10(m_eventCounter));
    statSize=(int)std::pow(10.,power);
  }
  if(m_eventCounter%statSize==0) return true;
  return false;
}

/** @brief Finalize Shared Data
    
    Writes histos and closes the TFile

    @return void
*/
void SharedData :: Finalize() 
{
  m_fout->cd();
  m_tree->Write();
  for( auto& h : m_v_hists ) { h->Write(); }
  
  m_fout->Close();
}

/*
void SharedData :: AddDoubleToTree( const std::string& name, std::vector<double> *pObj ){
  std::cout << "Adding " << name << std::endl;
  m_tree->Branch( name.c_str(), "std::vector<double>", &pObj );
}

void SharedData :: AddIntToTree( const std::string& name, std::vector<int> *pObj ){
  std::cout << "Adding " << name << std::endl;
  m_tree->Branch( name.c_str(), "std::vector<int>", &pObj );
}
*/


//***********************************************************************************
// READ Survey_2018.xml and Alignment_2018.xml
//***********************************************************************************
/**
 * @brief Reads the .xml configuration file and load characteristics for all the alignments, immediately sorted into alignment objects
 * @param _inFile
 */
 
void SharedData::LoadConfigurationFile( int m_runNumber, std::string _inFile  ){

    m_XMLparser = new XMLSettingsReader();
	
	 
    if (!m_XMLparser->parseFile(_inFile)) {
            std::cerr << " Data Reader could not parse file : " << _inFile << std::endl;
            return;
    }
	
    std::cout << "Loading .xml Configuration File..." << std::endl;
    std::cout << "Found " << m_XMLparser->getBaseNodeCount("Survey") << " survey entries " << std::endl;

    int first_run, last_run;

    for (int i = 0; i < m_XMLparser->getBaseNodeCount("Survey"); i++) { //this was unsigned int i = 0
        
        m_XMLparser->getChildValue("Survey",i,"start_run",first_run);
        m_XMLparser->getChildValue("Survey",i,"end_run",last_run);

        //Discard entries for any alignment that does not apply to our run
        if(m_runNumber < first_run || m_runNumber > last_run) continue;
		m_survey = new Survey();

        //If the entry applies, we store it in the vector
        m_XMLparser->getChildValue("Survey",i,"detector",m_survey->detector);
        m_XMLparser->getChildValue("Survey",i,"x_pos",m_survey->x_pos);
        m_XMLparser->getChildValue("Survey",i,"y_pos",m_survey->y_pos);
        m_XMLparser->getChildValue("Survey",i,"z_pos",m_survey->z_pos);
        m_XMLparser->getChildValue("Survey",i,"cos_x",m_survey->cos_x);
        m_XMLparser->getChildValue("Survey",i,"cos_y",m_survey->cos_y);
        m_XMLparser->getChildValue("Survey",i,"cos_z",m_survey->cos_z);
  
		surveyEntries.push_back(m_survey);
    }

	if(surveyEntries.size() == 0) std::cout << "WARNING: SURVEY NOT FOUND!!!" << std::endl;		

	return;
}

/**
 * @brief Reads the .xml configuration file and load characteristics for all the channels, immediately sorted into detectors objects
 * @param _inFile
 */
void SharedData::LoadAlignmentFile( int m_runNumber, std::string _inFile ){

    m_XMLparser = new XMLSettingsReader();

    if (!m_XMLparser->parseFile(_inFile)) {
            std::cerr << " Data Reader could not parse file : " << _inFile << std::endl;
            return;
    }

    m_alignment = new Alignment();

    std::cout << "Loading .xml Alignment File..." << std::endl;
    std::cout << "Found " << m_XMLparser->getBaseNodeCount("Alignment") << " alignment entries " << std::endl;
    std::cout << "Retrieving the information for run " << m_runNumber << std::endl;

    int run;
    for ( int i = 0; i < m_XMLparser->getBaseNodeCount("Alignment"); i++) {
        m_XMLparser->getChildValue("Alignment",i,"run",run);
        if(run != m_runNumber) continue;
        std::cout << "Found Run Entry in Alignment file for run " << m_runNumber << std::endl;
        m_XMLparser->getChildValue("Alignment",i,"x_table",m_alignment->x_table);
        m_XMLparser->getChildValue("Alignment",i,"y_table",m_alignment->y_table);
        m_XMLparser->getChildValue("Alignment",i,"upstream_Det",m_alignment->upstream_Det);
        m_XMLparser->getChildValue("Alignment",i,"mid_Det",m_alignment->mid_Det);
        m_XMLparser->getChildValue("Alignment",i,"downstream_Det",m_alignment->downstream_Det);
        m_XMLparser->getChildValue("Alignment",i,"target_In",m_alignment->target_In);
        m_XMLparser->getChildValue("Alignment",i,"lead_In",m_alignment->lead_In);
        m_XMLparser->getChildValue("Alignment",i,"magnet_On",m_alignment->magnet_On);
    }

    if(m_alignment == NULL) std::cout << "WARNING: ALIGNMENT NOT FOUND!!!" << std::endl;
    return;
}

Survey* SharedData::GetSurvey(std::string name){
	for(unsigned int i = 0; i < surveyEntries.size(); i++){	
		if( name == surveyEntries[i]->detector ){ return surveyEntries[i]; }
	}
	Survey* empty=NULL;
	return empty;
}

Alignment* SharedData::GetAlignment(){
	return m_alignment;
}
