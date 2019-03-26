/** @file DataReader
 *  @brief Function prototypes for DataReader
 *
 *  This contains the prototypes and members 
 *  for DataReader
 *
 *  @author Yakov Kulinich
 *  @bug No known bugs.
 */

#ifndef DATAREADER_H
#define DATAREADER_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

#include "XMLSettingsReader.h"
#include "Containers.h"
#include "Detector.h"

#include <TChain.h>

class TFile;
class Analysis;

class DataReader{

 public:
  DataReader( );
  DataReader( const unsigned int = 0,  const unsigned int = 0 );
  DataReader( const unsigned int = 0,  const unsigned int = 0,
              const std::string& = "" );
  DataReader( const unsigned int = 0,  const unsigned int = 0,
              const std::string& = "", const unsigned int = 0 );
  virtual ~DataReader();

  void AddAnalysis  ( Analysis* );
  
  void ReadListOfFiles( std::string listname );

  void LoadAlignmentFile     (std::string _inFile = std::getenv("JZCaPA") + std::string("/Utils/Alignment_2018.xml"));
  void LoadConfigurationFile (std::string _inFile = std::getenv("JZCaPA") + std::string("/Utils/ConfigFile2018.xml"));
  void SetDebugMode          ( ) { m_debug = true; }
  void SetVerbosity          ( int _level ){ m_verbose = _level; }
  void UpdateConsole         ( );

  Detector* GetDetector( std::string _detName );

  void Run();
  
  void Initialize   ( );
  void ProcessEvents( );
  void Finalize     ( );
  
 private:
  // output file
  TFile* m_fOut;
  
  // vector of all analysis
  std::vector< Analysis* > m_ana;

  //Number of channels to be read
  unsigned int m_nCh;
  //Number of samples per channel
  unsigned int m_nSamp;

  //Input file name
  std::string m_fNameIn;
  //Input list of files
  std::string m_fListOfFiles;

  //Run number
  unsigned int m_runNumber;
  
  //Boolean switch to enable the reading of a list of files
  bool m_readListOfFiles;

  //Input file (in case of a single processing)
  TFile* m_fIn;
  //TChain to accomodate many files (in case of a list of files)
  TChain* m_fileChain;

  //Vector of detectors placed in the 2018 setup (2 ZDCs, 1 RPD)
  std::vector < Detector* > m_detectors;

  //Alignment information for the given run
  Alignment* m_alignment;

  //XML parser
  XMLSettingsReader *m_XMLparser;
  
  //Current event
  int m_event;
  //Event number of last update
  int m_event_old = 0;
  //Rate of console updates in ms
  int m_update_rate = 1000;

  //DebugVariable
  bool m_debug = false;
  //Verbosity level
  int m_verbose = 0;
};

#endif
