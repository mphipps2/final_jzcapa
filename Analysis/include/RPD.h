/** @file RPD
 *  @brief Function prototypes for RPD
 *
 *  This contains the prototypes and members 
 *  for RPD
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
 */

#ifndef RPD_H
#define RPD_H

#include "Detector.h"

class RPD : public Detector{

 public:
  RPD( );
  virtual ~RPD( );
  
  virtual void  PrintMap  ( );
  
    };

#endif
