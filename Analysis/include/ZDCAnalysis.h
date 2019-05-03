/** @file ZDCAnalysis
 *  @brief Function prototypes for ZDCAnalysis
 *
 *  This contains the prototypes and members 
 *  for ZDCAnalysis
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
 */

#ifndef ZDCANALYSIS_H
#define ZDCANALYSIS_H

#include "Analysis.h"
#include "TH2D.h"
#include "Containers.h"
#include "Detector.h"

class ZDCAnalysis : public Analysis{

 public :
  ZDCAnalysis( );
  virtual ~ZDCAnalysis( );
  
  virtual void   Initialize     ( );
  virtual void   Initialize     ( std::vector < Detector* > _vDet);
  virtual void   SetupHistograms( );
  virtual void   AnalyzeEvent   ( );
  virtual void   AnalyzeEvent   ( const std::vector< TH1* >& ){};
  virtual void   AnalyzeEvent   ( const std::vector< std::vector< float > >& ){};
  virtual void   AnalyzeEvent   ( const std::vector< Channel* > ){};
  virtual void   Finalize       ( );

  
  /** Charge correlation */
  TH2D *hCharge;
  /** Peak height correlation */
  TH2D *hPeak;
  
 private :
  /** Pointer to the first ZDC module */
  Detector *m_zdc1 = 0;
  /** Pointer to the second ZDC module */
  Detector *m_zdc2 = 0;
  /** Pointer to the first ZDC channel */
  Channel* zdc1 = 0;
  /** Pointer to the second ZDC channel */
  Channel* zdc2 = 0;
 
  
};

#endif