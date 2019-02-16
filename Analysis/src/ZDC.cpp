/** @file ZDC.cpp
 *  @brief Implementation of ZDC.
 *
 *  Function definitions for ZDC are provided. 
 *  This is a daughter class of Detector.
 *  Methods specific to ZDCs are implemented here.
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
 */

#include "ZDC.h"


/** @brief Default Constructor for ZDC.
 */
ZDC::ZDC( ){
	
}


ZDC::ZDC( std::vector < Channel* > _readOut, int _zdcNumber){

    m_Number = _zdcNumber;
    for(int i = 0; i < (int)_readOut.size(); i++){
        if((_readOut.at(i)->mapping_row == 0) && (_readOut.at(i)->mapping_column == m_Number)){
            SetElement(_readOut.at(i));
        }
    }
    if(GetChannelsVector().size() > 1) std::cout << "WARNING : more than one entry for one ZDC module. Check the config.xml" << std::endl;
    std::cout << "ZDC object created with " << GetChannelsVector().size() << " channel entries " << std::endl;
}

/** @brief Destructor for ZDC.
 */
ZDC::~ZDC( ){

}

/** @brief Prints a map of the ZDC to the terminal
 *
 * Prints a "map" of the ZDC.
 * Displays one element with ZDC number on the top line, DRS4
 * Channel on the second line, and if the element is functioning on the
 * third line.
 * 
 */
void ZDC::PrintMap(){
        Channel *c = GetElement(0,m_Number);
	
        std::cout<<"|   "     << m_Number  <<"  |"<<std::endl;
	std::cout<<"|   "     << c->name <<"  |"<<std::endl;
	std::cout<<"|       |"<< std::endl;
	std::cout<<"|_______|"<< std::endl;

}
