/** @file WFAnalysis
 *  @brief Function prototypes for WFAnalysis
 *
 *  This contains the prototypes and members 
 *  for WFAnalysis
 *
 *  @author Yakov Kulinich
 *  @bug No known bugs.
 */

#ifndef WFANALYSIS_H
#define WFANALYSIS_H

#include "Analysis.h"

class WFAnalysis : public Analysis{

 public :
  WFAnalysis( );
  WFAnalysis( const std::string& = "" );
  virtual ~WFAnalysis( );

  virtual void Initialize     ( );
  virtual void SetupHistograms( );
  virtual void AnalyzeEvent   ( const std::vector< TH1* >& );
  virtual void AnalyzeEvent   ( const std::vector< std::vector< float > >& );
  virtual void Finalize       ( );

};

#endif
