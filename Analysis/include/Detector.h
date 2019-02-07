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
#include <string>
#include "Containers.h"

class Detector{

 public:
  Detector( ){};
  Detector( std::vector< Channel > _element ){ Element = _element}
  virtual ~Detector( ){};

  
  virtual Channel* GetElement  (int row, int column,int runNum, std::string type) = 0;
  virtual float*   GetPosition ( ) { return Position; };
  virtual float*   GetAngle    ( ) { return Angle; };
  virtual void     SetElement  ( Channel _entry) { element.push_back(_entry); };
  virtual void     SetPosition (double x, double y, double z) { Position[0] = x; Position[1] = y; Position[3] = z; };
  virtual void     SetAngle    (double a = 0, double b = 0, double c = 0) { Angle[0] = a; Angle[1] = b; Angle[3] = c; };
  virtual void     PrintMap    ( ) = 0;
  
  private:
  std::vector< Channel > Element;
  /** Three element array with x, y, and z of some pre-defined point on the detector **/
  double Position[3];
  /** Three element array of angle about the x, y, and z axis **/
  double Angle[3];
  };

#endif
