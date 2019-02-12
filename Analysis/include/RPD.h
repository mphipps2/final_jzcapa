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
  
  void SetnRows(int rows){nRows=rows; nElements = nRows * nColumns; ResizeSortedElements(); };
  void SetnCols(int cols){nColumns=cols; nElements = nRows * nColumns; ResizeSortedElements(); };
  void ResizeSortedElements();
  
  virtual void PrintMap( );
  
 private:
  /** Number of rows **/
  int nRows = 4;
  /** Number of columns **/
  int nColumns = 4;
  /** Total elements **/
  int nElements = nRows * nColumns;
  /** 2D vector of channels sorted in a [row][column] format **/
  std::vector< std::vector< Channel* > > m_SortedElements;
  
};

#endif
