/** @file RPDAnalysis
 *  @brief Function prototypes for RPDAnalysis
 *
 *  This contains the prototypes and members 
 *  for RPDAnalysis
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
 */

#ifndef RPDANALYSIS_H
#define RPDANALYSIS_H

#include "Analysis.h"
#include "TH2D.h"
#include "Containers.h"
#include "RPD.h"

class RPDAnalysis : public Analysis{

 public :
  RPDAnalysis( );
  virtual ~RPDAnalysis( );
  
  virtual void   Initialize     ( ){};
  virtual void   Initialize     ( std::vector < Detector* > _vDet );
  virtual void   SetupHistograms( );
  virtual void   AnalyzeEvent   ( );
  virtual void   AnalyzeEvent   ( const std::vector< TH1* >& ){};
  virtual void   AnalyzeEvent   ( const std::vector< std::vector< float > >& ){};
  virtual void   AnalyzeEvent   ( const std::vector< Channel* > ){};
  virtual void   SetBranches    ( TTree* _tree );
  virtual void   Finalize       ( );

  
  /** Running average of charges */
  double m_charge[5][5];
  /** Running average of charges */
  double m_peak[5][5];
  /** Center of mass x */
  double xCoM;
  /** Center of mass y */
  double yCoM;
  /** Sum of all channel charges */
  double ChargeSum;
  /** Sum of all channel peak heights */
  double PeakSum;
  /** Average charge per tile */
  TH2D *hCharge;
  /** Average peak height per tile */
  TH2D *hPeak;
  /** Calculated center of mass */
  TH2D *hCenter;
  /** Center of tiles in X mm */
  double xPos[4] = {32.34,10.79,-10.79,-32.37};
  /** Gap between tiles in Y mm */
  double yPos[4] = {30.75,10.25,-10.25,-30.75}; 
  
 private :
  /** Pointer to the RPD */
  RPD *m_RPD = 0;
  /** Array of pointers to the RPD channels to avoid repeated GetElement() calls */
  Channel* rpd[5][5];
 
  
};

#endif