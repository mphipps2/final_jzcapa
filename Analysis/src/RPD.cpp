/** @file RPD.cpp
 *  @brief Implementation of RPD.
 *
 *  Function definitions for RPD are provided. 
 *  This is a daughter class of Detector.
 *  Methods specific to RPDs are implemented here.
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
 */

#include "RPD.h"

#include <string>
#include <stdio.h>

/** @brief Default Constructor for RPD.
 */
RPD::RPD( ){

}

/** @brief Constructor that takes the whole vector of channels readout and selects and stores only the RPD ones
*/

RPD::RPD( std::vector < Channel* > _readOut){

    ResizeSortedElements();
    for(int row=0; row<nRows; row++){
        for(int column; column<nColumns; column++){
            m_SortedElements[row][column]=GetElement(row,column);
        }
    }

    std::cout << "RPD object created with " << nRows << " rows and " << nColumns << " columns" << std::endl;
}

/** @brief Destructor for RPD.
 */
RPD::~RPD( ){

}

/** @brief Get the properties of a detector element
 *
 * If m_SortedElements is populated, return the element in 
 * m_SortedElements[row][column].
 * Otherwise returns a pointer to the Channel stored in m_Element 
 * with the requested row and column.
 * If the requested element is not found, return a warning message and a NULL pointer.
 * 
 */
Channel* RPD::GetElement(int row, int column){
    
    if((int)m_SortedElements.size()==nElements){ return m_SortedElements[row][column];
    }else{
        for(int i=0; i < (int)GetChannelsVector().size(); i++){
            if(row == GetChannelsVector().at(i)->mapping_row && column == GetChannelsVector().at(i)->mapping_column){
            return GetChannelsVector().at(i);
            }
        }
     }
  std::cerr << " WARNING: Element (" << row << "," << column << ") not found! " << std::endl;
  return nullptr;
    
}

/** @brief Prints a map of the RPD to the terminal
 *
 * Prints the map of a 4x4 RPD.
 * Displays a grid of elements with row,column on the top line, DRS4
 * Channel on the second line, and if the element is functioning on the
 * third line.
 * 
 */
void RPD::PrintMap(){
  std::cout << " ___________________________________ " << std::endl;
  //RPD has 4 rows and 4 columns
  for(int row = 1; row <= 4; row++){
      Channel* c[4];
      std::string status, name;
      for(int cln = 1; cln <= 4; cln++){
          
          c[cln] = GetElement(row,cln);
          if(c[cln]->is_on) status = "ON";
              else status = "OFF";
          if (cln == 1) std::cout << "|  " << row << "," << cln << " --> " << c[cln]->name << " , " << status;
          else std::cout << "  |  " << row << "," << cln << " --> " << c[cln]->name << " , " << status;
          
        }//End of the loop over columns
      }//End of the loop over rows
}

/** @brief Resizes m_SortedElements
 *
 * Uses nRows and nColumns to resize m_SortedElements
 * 2D vectors must have columns sized individually.
 * 
 */
void RPD::ResizeSortedElements(){
    m_SortedElements.resize(nRows);
    for(int row=0; row < nRows; row++){
        m_SortedElements.resize(nColumns);
    }
}








