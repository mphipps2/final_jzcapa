/** @file Analysis
 *  @brief Function prototypes for Analysis
 *
 *  This contains the prototypes and members 
 *  for Analysis
 *
 *  @author Yakov Kulinich
 *  @bug No known bugs.
 */

#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <string>
#include <vector>

class TFile;
class TH1;
class TH2;
class TH3;
class TTree;

class Analysis{

 public:
  Analysis( ){};
  Analysis( const std::string& ){};
  virtual ~Analysis( ){};

  virtual void Initialize     ( ) = 0;
  virtual void SetupHistograms( ) = 0;
  virtual void AnalyzeEvent   ( const std::vector< TH1* >& ) = 0;
  virtual void AnalyzeEvent   ( const std::vector< std::vector< float > >& ) = 0;
  virtual void Finalize       ( ) = 0;

 protected:
  std::string m_fNameOut;
  TFile* m_fOut;

};

#endif
