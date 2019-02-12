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
  
  Channel* GetElement(int row, int column);
  
  void SetnRows(int rows){nRows=rows; nElements = nRows * nColumns; };
  void SetnCols(int cols){nColumns=cols; nElements = nRows * nColumns; };
  
  virtual void PrintMap( );
  
 private:
  int nRows = 4;
  int nColumns = 4;
  int nElements = nRows * nColumns;
  
    };

#endif
