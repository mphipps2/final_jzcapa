 /** @file DataReader.cxxs
 *  @brief Implementation of DataReader.
 *
 *  Function definitions for DataReader are provided. 
 *  This class reads a rootfile with raw waveforms 
 *  that are processed by rcdaqAnalysis running on prdf files.
 *  Then, in the event loop, analysis classes can be called.
 *
 *  @author Yakov Kulinich, Riccardo Longo
 *  @bug No known bugs.
 */

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TChain.h>

#include <iostream>

#include "DataReader.h"
#include "Analysis.h"
#include "Containers.h"
#include "RPD.h"
#include "ZDC.h"

/** @brief Default Constructor for DataReader.
 */
DataReader::DataReader() : DataReader( 0, 0, "", 0 ){

}

/** @brief Constructor for DataReader.
 *
 *  @param1 Number of channels being read
 *  @param2 Number of samples per channel
 */
DataReader::DataReader( const unsigned int nCh, const unsigned int nSamp )
  : DataReader( nCh, nSamp, "", 0 ){

  // here say we will read in list of files, maybe
  // there is better way to do it, just an idea for now.
  m_readListOfFiles = true;
}

/** @brief Constructor for DataReader.
 *
 *  @param1 Number of channels being read
 *  @param2 Number of samples per channel
 *  @param3 Input filename.
 */
DataReader::DataReader( const uint nCh, const uint nSamp,
			const std::string& fNameIn )
  : DataReader( nCh, nSamp, fNameIn, 0 ){
  
}


/** @brief Constructor for DataReader.
 *
 *  @param1 Number of channels being read
 *  @param2 Number of samples per channel
 *  @param4 Output file name (custom)
 *  @param3 Run number being used.
 
 */
DataReader::DataReader( const uint nCh, const uint nSamp,
			const std::string& fNameIn, const uint runNum )
  : m_nCh( nCh ), m_nSamp( nSamp ), m_fNameIn( fNameIn ),
    m_runNumber( runNum ), m_readListOfFiles( false ), m_fIn( NULL ){

}


/** @brief Destructor for DataReader.
 */
DataReader::~DataReader(){

  for( auto& ana : m_ana ){
    delete ana; ana = NULL;
  }
}

/** @brief Adds an analysis to vector of analysis
 *
 *  @param1 Pointer to an Analysis.
 *
 *  @return none
 */
void DataReader::AddAnalysis( Analysis* ana ){
  m_ana.push_back( ana );
}


/** @brief Enables the read from list of files option for DataReader
 *
 *  @param1 name of the list of files to be processed (with the full path if it's not in the execution folder)
 *
 *  @return none
 */
void DataReader::ReadListOfFiles( std::string listname ){

    m_readListOfFiles = true;
    m_fListOfFiles = listname;
}


/**
 * @brief Reads the .xml configuration file and load characteristics for all the channels, immediately sorted into detectors objects
 * @param _inFile
 */
void DataReader::LoadAlignmentFile(std::string _inFile ){

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
    for (unsigned int i = 0; i < m_XMLparser->getBaseNodeCount("Alignment"); i++) {
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

/**
 * @brief Reads the .xml configuration file and load characteristics for all the channels, immediately sorted into detectors objects
 * @param _inFile
 */
void DataReader::LoadConfigurationFile(std::string _inFile ){

    m_XMLparser = new XMLSettingsReader();

    if (!m_XMLparser->parseFile(_inFile)) {
            std::cerr << " Data Reader could not parse file : " << _inFile << std::endl;
            return;
    }

    std::cout << "Loading .xml Configuration File..." << std::endl;
    std::cout << "Found " << m_XMLparser->getBaseNodeCount("channel") << " channel entries " << std::endl;

    std::vector < Channel* > channelEntries;
    int first_run, last_run;

    for (unsigned int i = 0; i < m_XMLparser->getBaseNodeCount("channel"); i++) {
        Channel *buffer_ch = new Channel();
        m_XMLparser->getChildValue("channel",i,"start_run",first_run);
        m_XMLparser->getChildValue("channel",i,"end_run",last_run);

        //Discard entries for any channel that does not apply to our run
        if(m_runNumber < first_run || m_runNumber > last_run) continue;


        //If the entry applies, we store it in the vector
        m_XMLparser->getChildValue("channel",i,"detector",buffer_ch->detector);
        m_XMLparser->getChildValue("channel",i,"name",buffer_ch->name);
        m_XMLparser->getChildValue("channel",i,"mapping_row",buffer_ch->mapping_row);
        m_XMLparser->getChildValue("channel",i,"mapping_column",buffer_ch->mapping_column);
        m_XMLparser->getChildValue("channel",i,"delay",buffer_ch->delay);
        m_XMLparser->getChildValue("channel",i,"offset",buffer_ch->offset);
        m_XMLparser->getChildValue("channel",i,"HV",buffer_ch->HV);
        m_XMLparser->getChildValue("channel",i,"is_on",buffer_ch->is_on);
        m_XMLparser->getChildValue("channel",i,"Vop",buffer_ch->Vop);

        bool isNew(true);
        for( int k = 0; k < channelEntries.size(); k++){
            if(buffer_ch->name == channelEntries.at(k)->name){
                std::cout << "WARNING!!! Redundancy in your settings file for " << buffer_ch->name << ". Check it carefully. The second entry found will be skipped..." << std::endl;
                isNew = false;
            }
        }
        if(isNew) channelEntries.push_back(buffer_ch);
    }

    std::cout << "Loaded " << channelEntries.size() << " configuration entries " << std::endl;
    if( channelEntries.size() < 18 ) std::cout << "WARNING!!!! Number of Channels < 18. Seems that some entry is missed for this run in the config.xml. BE CAREFUL!" << std::endl;
    ZDC* zdc1 = new ZDC(channelEntries,1);
    ZDC* zdc2 = new ZDC(channelEntries,2);
    RPD* rpd = new RPD(channelEntries);

    m_detectors.push_back(zdc1); //Position 0 goes for ZDC1
    m_detectors.push_back(zdc2); //Position 1 goes for ZDC2
    m_detectors.push_back(rpd); //Position 2 goes for the RPD

    std::cout << "Detector configuration: loading complete! " << std::endl;

    return;
}

/**
 * @brief DataReader::GetDetector allows the user to access the detectors objects after loading them at the beginning of the execution
 * @param _detName can be ZDC1 (upstream module), ZDC2 (downstream module) or RPD.
 * @return
 */
Detector* DataReader::GetDetector( std::string _detName ){

    if(_detName != "ZDC1" && _detName != "ZDC2" && _detName != "zdc1" && _detName != "zdc2" && _detName != "RPD" && _detName != "rpd")
    {
        std::cout << "The detector you're looking for is not ZDC1, ZDC2 or RPD. Please check and correct your request " << std::endl;
        return NULL;
    }
    if( _detName == "ZDC1" || _detName == "zdc1" ) return m_detectors.at(0);
    if( _detName == "ZDC2" || _detName == "zdc2" ) return m_detectors.at(1);
    if( _detName == "RPD" || _detName == "rpd" ) return m_detectors.at(2);

    std::cout << "WARNING: detector recognition glitch. NULL pointer being returned..." << std::endl;
    return NULL;
}


/** @brief Run method for DataReader
 *
 *  Call Initialize, ProcessEvents, and Finalize
 *  Made so the driver class only has to call one method.
 *  
 *  @return none
 */
void DataReader::Run(){

  Initialize();
  ProcessEvents();
  Finalize();
}

/** @brief Initialization method for DataReader
 *
 *  Select which file(s) to read. For now just a single
 *  file, later this can be extended to read many and make
 *  chain of files for example. Also create and initialize
 *  the analysis that will be running.
 *  
 *  @return none
 */
void DataReader::Initialize(){

  if( m_readListOfFiles ){

    // Riccardo - 21/01/2019 - TChain implementation
    // The file list name is supposed to be already set
    m_fileChain = new TChain("tree");
    std::ifstream inFile;
    inFile.open(m_fListOfFiles.c_str());
    std::string reader_buff;
    while(inFile >> reader_buff){
        //Let's push all the files in the list into the TChain
        m_fileChain->Add(reader_buff.c_str());
    }
    /** TODO - Add fileChain reading below
     */
  } else {
    m_fIn = TFile::Open( m_fNameIn.c_str() );
  }

  // If we are reading a list of files, or have no run number
  // make default name output.root, otherwise make it
  // outputN.root, where N is a run number of a file.
  std::string fNameOut = m_readListOfFiles ?
    "output.root" : Form( "output%d.root", m_runNumber );

  m_fOut = new TFile( fNameOut.c_str(), "RECREATE" );
  
  for( auto& ana : m_ana ){
    ana->Initialize();
    ana->SetupHistograms();
  }
}



/** @brief Process Events method for DataReader
 *
 *  Connect to files / chain of files.
 *  Read tree event by event, and fill basic waveforms.
 *  So far (12.12.18) For each event - Read Raw data for all channels,
 *  put it into 2D vector size NxM where N = nCh, M = nSamp
 *  Also, For each event - Fill N 1D histos that have size M.
 *  Then one can send these to any Analysis function they make.
 *  Here for example (12.12.18), we send to WFAnalysis::AnalyzeEvent( .. )
 *  See below - Have fun!
 *  
 *  @return none
 */
void DataReader::ProcessEvents(){
  
  //Added for debug/demo purposes. Remove after implementation of Visualizer
  TCanvas *canvas = new TCanvas( "Diff Demo", "Diff Demo", 200, 10, 1000, 600);
  TPad *pad = new TPad("pad", "pad",0.15,0.11,0.85,0.79);
  canvas->Divide(4,2);
  
  // Processed Raw data to read in as vector of vectors size NxM
  // Where N = nCh and M = nSamples per channel.
  std::vector< std::vector< float >  >  vWF;
  std::vector< std::vector< float >* > pvWF;

  // Histograms (N of them) for the raw waveforms from each event.
  // They will go to AnalyzeEvent for processing
  std::vector< TH1* > vWFH;

  // Resize these to be of size nCh.
  vWF .resize( m_nCh );
  pvWF.resize( m_nCh );
  vWFH.resize( m_nCh );

  /** TODO : add reading for list of files
    * Please note that many of the implementations are now for a single-file treatment
    */
  TTree* tree = static_cast< TTree* >( m_fIn->Get( "tree" ) );

  //Specific pointers to each detector, if needed afterwards
  ZDC* zdc1 = static_cast< ZDC* >( GetDetector("ZDC1") );
  ZDC* zdc2 = static_cast< ZDC* >( GetDetector("ZDC2") );
  RPD* rpd =  static_cast< RPD* >( GetDetector("RPD") );

  //All the raw channels addresses set for read-out
  for( uint detID = 0; detID < (int) m_detectors.size(); detID++ ){
       m_detectors.at(detID)->SetBranches(tree);
       m_detectors.at(detID)->SetNSamples(m_nSamp);
       m_detectors.at(detID)->DeclareHistograms();
  }

  // Connect raw data to tree. For the moment, the only reading implemented is the raw data from each channel.
  // Other items should be implemented in the same way if needed.
  // Also create the histograms to fill
  for( uint ch = 0; ch < m_nCh; ch++ ){
    //Example - here we retrieve the already processed waveform (M.Phipps approach during the data production)
    pvWF[ ch ] = &vWF[ ch ];
    tree->SetBranchAddress( Form( "C%d", ch ), &pvWF[ ch ] );
    vWFH[ ch ] = new TH1D( Form( "hWF%d", ch ), Form( "hWF%d;samp;amp", ch ), m_nSamp, 0, m_nSamp );
  }

  std::cout << "File: " << m_fIn->GetName() << " has " << tree->GetEntries() << " events." << std::endl;
  
  // !! EVENT LOOP
  for( int ev = 0; ev < tree->GetEntries(); ev++ ){
    
    tree->GetEntry( ev );
    
    std::cout << "Event : " << ev << std::endl;

    // Fill the raw waveforms
    for( uint detID = 0; detID < (int) m_detectors.size(); detID++ )
        m_detectors.at(detID)->FillHistograms();

    std::cout << "Event : " << ev << std::endl;

   //Here if you're interested in already processed waveform
   for( uint ch = 0; ch < m_nCh; ch++ ) {
   // Loop over samples in each channel
    for( uint samp = 0; samp < m_nSamp; samp++ ){
      vWFH[ ch ]->SetBinContent( samp + 1, vWF[ ch ][ samp ] );
      } // End loop over samples in each channel
    } // End loop over channels

    // Now call all analysis and run their AnalyzeEvent.
    // Can either send a vector of histos, or a 2D vector
    // Of all the data, depending on what you want to do.
    // Note that at the moment none of this methods is doing anything
    for( auto& ana : m_ana ){
      //raw data analysis
      if(m_debug){    
        // Uncomment to run a few events at a time
        //if(ev==8) break;
        // Uncomment to run a single event
        if(ev!=8) continue;
        
        //ana->AnalyzeEvent( zdc1->GetChannelsVector(), canvas->cd() );
        //ana->AnalyzeEvent( zdc2->GetChannelsVector(), canvas->cd() );
        ana->AnalyzeEvent( rpd->GetChannelsVector() , canvas->cd() );
      }else{
      ana->AnalyzeEvent( zdc1->GetChannelsVector() );
      ana->AnalyzeEvent( zdc2->GetChannelsVector() );
      ana->AnalyzeEvent( rpd->GetChannelsVector()  );
      //already processed wf analysis
      //ana->AnalyzeEvent( vWFH );
      //ana->AnalyzeEvent( vWF  );
      }
    }
  } // End event loop
  
  
  //Added for debug/demo purposes. Remove after implementation of Visualizer
  if(m_debug){
    pad->Update();
    canvas->Draw();
    canvas->Print( "Output.pdf" );
  }

  for( auto& h : vWFH ){ delete h; }
}

/** @brief Finalize method for DataReader
 *
 *  Close any input files 
 *  Call analysis finalize methods
 *  
 *  @return none
 */
void DataReader::Finalize(){
  
  if( m_fIn ){
    m_fIn->Close();
  }

  // enter the output file since
  // we will be writing to it now.
  m_fOut->cd();

  for( auto& ana : m_ana ){
    ana->Finalize();
  }

  m_fOut->Close();
}

