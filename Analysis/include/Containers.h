/** @file Containers.h
 *  @brief Class to define containers: child classes with only public data memebers
 *
 *
 *  @author Riccardo Longo
 *  @bug No known bugs.
 */

#ifndef CONTAINERS_H
#define CONTAINERS_H

#include "Containers.h"

#include <iostream>
#include <vector>

class Channel {

 public :
    /** Detector associated to this channel. RPD or ZDC
    *   Mapping (sitting on the impinging beam) [R,C] for ZDC:
    *   ZDC1 -> [0,1]    |
    *   ZDC2 -> [0,2]   \ /
    *   Mapping (sitting on the impinging beam) [R,C] for RPD:
    *   Jura [1,4] [1,3] [1,2] [1,1] Saleve
    *        [2,4] [2,3] [2,2] [2,1]
    *        [3,4] [3,3] [3,2] [3,1]
    *        [4,4] [4,3] [4,2] [4,1]
    */
	/** Type of detector - ZDC or RPD **/
    std::string detector;
    /** Channel name - CXX with XX in [1,20]  */
    std::string name;
    /** High voltage set for this channel - in V */
    double HV;
    /** Delay set for this channel - in ns */
    double delay;
    /** Offset set for this channel */
    double offset;
    /** Mapping in the horizontal direction [R] */
    int mapping_row;
    /** Mapping in the vertical direction [C] */
    int mapping_column;
    /** Was the channel functioning */
    bool is_on;
    /** Operating voltage of the channel**/
    int Vop;
    /** Raw waveform for a particular event **/
    std::vector < float > WF;

};

#endif
