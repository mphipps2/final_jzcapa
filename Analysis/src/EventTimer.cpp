/** @file EventTimer.cxx
 *  @brief Implementation of EventTimer.
 *
 *  Function definitions for EventTimer are provided. 
 *  This class is a reimplementation of TTimer with
 *  a new notify method for JZCaPA
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
 */

#include "EventTimer.h"

/** @brief Default constructor for EventTimer
 *
 */
EventTimer::EventTimer(Long_t milliSec, DataReader* obj, Bool_t mode)
    : TTimer(milliSec, mode){
    m_object = obj;
}

/** @brief Default destructor for EventTimer
 *
 */
EventTimer::~EventTimer(){
    
}


/** @brief Notify when timer times out.
 *
 *  Calls the update method from the object pointer
 *  The timer is always reset. To stop the timer 
 *  call TurnOff().
 */
Bool_t EventTimer::Notify(){
    Timeout();       // emit Timeout() signal
    if (m_object) m_object->UpdateConsole();
    
    Reset();
    return kTRUE;
 }