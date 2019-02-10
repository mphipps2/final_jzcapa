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

    for(int i = 0; i < (int)_readOut.size(); i++){
        if(_readOut.at(i)->name == "RPD"){
            SetElement(_readOut.at(i));
        }
    }
    std::cout << "RPD object created with " << GetChannelsVector().size() << " channel entries " << std::endl;
}

/** @brief Destructor for RPD.
 */
RPD::~RPD( ){

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
