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
