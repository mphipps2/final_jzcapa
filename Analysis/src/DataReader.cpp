/** @file DataReader.cxx
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

#include <iostream>

#include "DataReader.h"
#include "WFAnalysis.h"

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
  : m_nCh( nCh ), m_nSamp( nSamp ),
    m_fNameIn( fNameIn ), m_runNumber( runNum ),
    m_ana( NULL ), m_fIn( NULL ){
}


/** @brief Destructor for DataReader.
 */
DataReader::~DataReader(){
  delete m_fIn;
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

  // If we are reading a list of files, or have no run number
  // make default name output.root, otherwise make it
  // outputN.root, where N is a run number of a file.
  std::string fNameOut = m_readListOfFiles ?
    "output.root" : Form( "output%d.root", m_runNumber );

  if( m_readListOfFiles ){
    // for now do nothing. Was thinking of having
    // a file called inputFiles.txt, from which to
    // read a list of files. TBD (I am in Berlin now).
  } else {
      m_fIn = TFile::Open( m_fNameIn.c_str() );
  }
  
  m_ana = new WFAnalysis( fNameOut );
  m_ana->Initialize();
  m_ana->SetupHistograms();
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

  std::vector< std::vector< float >  >  vWF;
  std::vector< std::vector< float >* > pvWF;
  std::vector< TH1* > vWFH;
  vWF .resize( m_nCh );
  pvWF.resize( m_nCh );
  vWFH.resize( m_nCh );
  
  TTree* tree = static_cast< TTree* >
    ( m_fIn->Get( "tree" ) );

  // Connect raw data to tree
  // Also create the histograms to fill
  for( uint ch = 0; ch < m_nCh; ch++ ){
    pvWF[ ch ] = &vWF[ ch ];
    tree->SetBranchAddress( Form( "RawC%d", ch ), &pvWF[ ch ] );
    vWFH[ ch ]  = new TH1D( Form( "hWF%d", ch ), Form( "hWF%d;samp;amp", ch ), m_nSamp, 0, m_nSamp );
  }

  std::cout << "File: " << m_fIn->GetName() << " has "
	    << tree->GetEntries() << " events." << std::endl;
  
  // !! EVENT LOOP
  for( int ev = 0; ev < tree->GetEntries(); ev++ ){
    tree->GetEntry( ev );

    // Fill the waveforms
    for( uint ch = 0; ch < m_nCh; ch++ ){
      vWFH[ ch ]->Reset();

    // Loop over samples in each channel
    for( uint samp = 0; samp < m_nSamp; samp++ ){
      vWFH[ ch ]->SetBinContent( samp + 1, vWF[ ch ][ samp ] );
      } // End loop over samples in each channel
    } // End loop over channels

    // Now call analysis. With the provided WFAnalysis,
    // Can either send a vector of histos, or a 2D vector
    // Of all the data, depending on what you want to do.
    m_ana->AnalyzeEvent( vWFH );
    m_ana->AnalyzeEvent( vWF  );
  } // End event loop

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

  m_ana->Finalize();
}

