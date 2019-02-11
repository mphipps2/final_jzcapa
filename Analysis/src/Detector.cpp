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
#include <vector>
 
 
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

/** @brief Get the properties of a detector element
 *
 * Returns a pointer to the Channel stored in m_Elements with the
 * requested channel name.
 * If the requested element is not found, return a warning message and a NULL pointer.
 * 
 */
Channel* Detector::GetElement(std::string _name){
    
    for(int i=0; i < (int)m_Element.size(); i++){
            if(!_name.compare(m_Element.at(i)->name)) return m_Element[i];
  }
  std::cerr << " WARNING: Element (" << _name << ") not found! " << std::endl;
  return nullptr;
    
}

void Detector::SetBranches( TTree *_dataTree ){

    std::vector< std::vector< float >* > pvWF;
    pvWF.resize(m_Element.size());

    for( uint ch = 0; ch < m_Element.size(); ch++ ){
      pvWF[ ch ] = &m_Element[ch]->WF;
      _dataTree->SetBranchAddress( ("Raw" + m_Element.at(ch)->name).c_str(), &pvWF[ ch ] );
    }

}










