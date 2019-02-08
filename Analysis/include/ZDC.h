/** @file ZDC
 *  @brief Function prototypes for ZDC
 *
 *  This contains the prototypes and members 
 *  for ZDC
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
 */

#ifndef ZDC_H
#define ZDC_H

#include "Detector.h"

class ZDC : public Detector{

 public:
  ZDC( );
  virtual ~ZDC( );
  virtual ~ZDC( int _number ) { Number = _number; };
  
  virtual void PrintMap  ( );
  virtual void SetNumber ( int _number ) { Number = _number; };
 
 private:
  int Number;
};

#endif
