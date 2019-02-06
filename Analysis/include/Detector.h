/** @file Detector
 *  @brief Function prototypes for Detector
 *
 *  This contains the prototypes and members 
 *  for Detector
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
 */

#ifndef DETECTOR_H
#define DETECTOR_H

#include <vector>
#include "Containers.h"

class Detector{

 public:
  Detector( ){};
  virtual ~Detector( ){};
  
  virtual Channel GetElement  (int row, int column) = 0;
  virtual void    SetElement  (int row, int column) = 0;
  virtual double* GetPosition ( ) = 0;
  virtual void    SetPosition (double x, double y, double z) { Position[0] = x; Position[1] = y; Position[3] = z; };
  virtual void    PrintMap    ( ) = 0;
  
  private:
  std::vector< Channel > Element;
  /** Four element array with x, y, and z angle of some pre-defined point on the detector **/
  double Position[3];
  double Angle[3];
  };

#endif
