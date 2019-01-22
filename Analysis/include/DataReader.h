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

#include <TChain.h>

class TFile;
class Analysis;

class DataReader{

 public:
  DataReader( );
  DataReader( const unsigned int = 0, const unsigned int = 0 );
  DataReader( const unsigned int = 0, const unsigned int = 0, const std::string& = "" );
  DataReader( const unsigned int = 0, const unsigned int = 0,
	      const std::string& = "", const unsigned int = 0 );
  virtual ~DataReader();

  void ReadListOfFiles(std::string listname);

  void Initialize   ( );
  void ProcessEvents( );
  void Finalize     ( );

 private:
  Analysis* m_ana;

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
};

#endif
