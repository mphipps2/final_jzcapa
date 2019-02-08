/** @file Detector.cpp
 *  @brief Implementation of Detector.
 *
 *  Function definitions for Detector are provided. 
 *  This is the mother class for detectors. 
 *  Methods common to all detectors are implemented here.
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
 */
 
#include "Detector.h"
 
 
/** @brief Default Constructor for Detector.
 */
Detector::Detector( ){

}

/** @brief Destructor for Detector.
 */
Detector::~Detector( ){

}

/** @brief Get the properties of a detector element
 *
 * Returns a pointer to the Channel stored in Elements with the
 * requested row and column.
 * If the requested element is not found, return an NULL pointer.
 * 
 */
Channel* Detector::GetElement(int row, int column){
  int entries=Element.size();
  for(int i=0; i<entries; i++){
    if(row==Element[i].mapping_row && column==Element[i].mapping_column){
      return &Element[i];
    }
  }
  return nullptr;
}












