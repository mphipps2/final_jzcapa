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
#include "Containers.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TF1.h"
#include "TLine.h"


class WFAnalysis : public Analysis{

 public :
  WFAnalysis( );
  virtual ~WFAnalysis( );

  virtual void   Initialize     ( );
  virtual void   SetupHistograms( );
  virtual void   GetDifferential( TH1D *hIN, TH1D *hOUT, int N);
  virtual double GetRMS         ( TH1D *h, int diff_window, bool debug = false) ;
  virtual void   FindHitWindow  ( Channel* ch, double threshMultiple );
  virtual void   AnalyzeEvent   ( const std::vector< TH1* >& );
  virtual void   AnalyzeEvent   ( const std::vector< std::vector< float > >& );
  virtual void   AnalyzeEvent   ( const std::vector< Channel* > vCh );
  virtual void   Finalize       ( );

};

#endif
