/** @file Detector.cpp
 *  @brief Implementation of Detector.
 *
 *  Function definitions for Detector are provided. 
 *  This is the mother class for detectors. 
 *  Methods common to all detectors are implemented here.
 *
 *  @author Chad Lantz, Riccardo Longo
 *  @bug No known bugs.
 */
 
#include "Detector.h"
#include "Containers.h"
 
 
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
 * Returns a pointer to the Channel stored in m_Elements with the
 * requested row and column.
 * If the requested element is not found, return a warning message and a NULL pointer.
 * 
 */
Channel* Detector::GetElement(int row, int column){
    
    for(int i=0; i < (int)m_Element.size(); i++){
            if(row == m_Element[i]->mapping_row && column == m_Element[i]->mapping_column){
            return m_Element[i];
    }
  }
  std::cerr << " WARNING: Element (" << row << "," << column << ") not found! " << std::endl;
  return nullptr;
    
}












