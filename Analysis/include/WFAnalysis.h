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
  WFAnalysis( int _sensitivity, double _threshold ){ m_diffSens = _sensitivity; m_Tmultiple = _threshold; }
  virtual ~WFAnalysis( );

  virtual void   Initialize     ( );
  virtual void   SetupHistograms( );
  virtual void   GetDifferential( TH1D *hIN, TH1D *hOUT );
  virtual double GetRMS         ( TH1D *h, bool debug = false) ;
  virtual void   GetPedestal    ( Channel* ch );
  virtual void   FindHitWindow  ( Channel* ch );
  virtual void   AnalyzeEvent   ( const std::vector< TH1* >& );
  virtual void   AnalyzeEvent   ( const std::vector< std::vector< float > >& );
  virtual void   AnalyzeEvent   ( const std::vector< Channel* > vCh );
  virtual void   Finalize       ( );

 private :
  /** Sensitivity level to hits (differentiation window) */
  int m_diffSens = 25;
  /** Hit threshold multiplier */
  double m_Tmultiple = 3.5;
};

#endif
