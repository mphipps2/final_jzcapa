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


/** @brief Default Constructor for RPD.
 */
RPD::RPD( ){
	
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
void PrintMap(){
  std::cout<<" _______________________________ "<<std::endl;
  for(int row=1;row<=4;row++){
    Channel *c[4];
    c[0]=GetElement(row,1,"RPD");
    c[1]=GetElement(row,2,"RPD");
    c[2]=GetElement(row,3,"RPD");
    c[3]=GetElement(row,4,"RPD");

    std::cout<<"|  "<<Form("%d,4",row);
    std::cout<<"  |  "<<Form("%d,4",row);
    std::cout<<"  |  "<<Form("%d,4",row);
    std::cout<<"  |  "<<Form("%d,4",row);
    std::cout<<"  |"<<std::endl;

    std::cout<<"|   "<<*c[3]->name;
    std::cout<<"  |   "<<*c[2]->name;
    std::cout<<"  |  "<<*c[1]->name;
    std::cout<<"   |   "<<*c[0]->name;
    std::cout<<"  |"<<std::endl;

    std::cout<<"|   "<<*c[3]->is_on;
    std::cout<<"  |   "<<*c[2]->is_on;
    std::cout<<"  |  "<<*c[1]->is_on;
    std::cout<<"   |   "<<*c[0]->is_on;
    std::cout<<"  |"<<std::endl;

    std::cout<<"|_______|_______|_______|_______|"<<std::endl;
  }	
}
