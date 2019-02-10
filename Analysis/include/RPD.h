/** @file RPD
 *  @brief Function prototypes for RPD
 *
 *  This contains the prototypes and members 
 *  for RPD
 *
 *  @author Chad Lantz, Riccardo Longo
 *  @bug No known bugs.
 */

#ifndef RPD_H
#define RPD_H

#include "Containers.h"
#include "Detector.h"

class RPD : public Detector{

 public:
  RPD( );
  RPD( std::vector< Channel* > _readOut );
  ~RPD( );
  
  virtual void PrintMap( );
  
    };

#endif
