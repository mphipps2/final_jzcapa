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
  Detector( std::vector< Channel > _element ){ Element = _element}
  virtual ~Detector( ){};

  
  virtual Channel GetElement  (int row, int column) = 0;
  virtual void    SetElement  ( Channel _entry) { element.push_back(_entry); };
  virtual double* GetPosition ( ) = 0;
  virtual void    SetPosition (double x, double y, double z) { Position[0] = x; Position[1] = y; Position[3] = z; };
  virtual void    SetAngle    (double a, double b, double c) { Angle[0] = a; Angle[1] = b; Angle[3] = c; };
  virtual void    PrintMap    ( ) = 0;
  
  private:
  std::vector< Channel > Element;
  /** Four element array with x, y, and z angle of some pre-defined point on the detector **/
  double Position[3];
  double Angle[3];
  };

#endif
