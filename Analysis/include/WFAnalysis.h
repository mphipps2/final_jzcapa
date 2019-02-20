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
#include "TVirtualPad.h"

class WFAnalysis : public Analysis{

 public :
  WFAnalysis( );
  virtual ~WFAnalysis( );

  virtual void Initialize     ( );
  virtual void SetupHistograms( );
  virtual TH1 *GetDifferential( TH1D *h, unsigned int ch, int window,bool debug=false);
  virtual void OverlayHistos  ( TH1D *h1, TH1D *h2 , TVirtualPad* pad, bool save);
  virtual void AnalyzeEvent   ( const std::vector< TH1* >& );
  virtual void AnalyzeEvent   ( const std::vector< std::vector< float > >& );
  virtual void AnalyzeEvent   ( const std::vector< Channel* > );
  virtual void Finalize       ( );

};

#endif
