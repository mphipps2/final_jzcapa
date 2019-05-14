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
  /** Sum of all channel differential peak heights */
  double DiffPeakSum;
  /** Average charge per tile */
  TH2D *hCharge;
  /** Average peak height per tile */
  TH2D *hPeak;
  /** Calculated center of mass */
  TH2D *hCenter;
  /** Charge vs peak height */
  TH2D *hChgVsPk;
  /** Peak height vs differential peak height */
  TH2D *hPkVsDiffPk;
  /** Charge sum histogram */
  TH1D *hChargeSum;
  /** Peak height sum histogram */
  TH1D *hPeakSum;
  /** Differential peak height sum histogram */
  TH1D *hDiffPeakSum;
  /** Array of charge histograms. One per tile */
  TH1D *hChargeArr[5][5];
  /** Array of peak height histograms. One per tile */
  TH1D *hPeakArr[5][5];
  /** Array of differential peak histograms. One per tile */
  TH1D *hDPeakArr[5][5];
  /** Center of tiles in X mm */
  double xPos[5] = {0,32.34,10.79,-10.79,-32.37};
  /** Gap between tiles in Y mm */
  double yPos[5] = {0,30.75,10.25,-10.25,-30.75}; 
  /** Estimated position of the beam calculated form table position */
  double beamPosX, beamPosY;
  
 private :
  /** Pointer to the RPD */
  RPD *m_RPD = 0;
  /** Array of pointers to the RPD channels to avoid repeated GetElement() calls */
  Channel* rpd[5][5];
 
  
};

#endif