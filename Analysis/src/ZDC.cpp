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
	Channel *c = GetElement(0,Number);
	
	std::cout<<"|   "     << Number  <<"  |"<<std::endl;
	std::cout<<"|   "     << c->name <<"  |"<<std::endl;
	std::cout<<"|       |"<< std::endl;
	std::cout<<"|_______|"<< std::endl;

}