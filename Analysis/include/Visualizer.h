/** @file Visualizer.h
 *  @brief Class to have a simple visualization of plots through JZCaPA. If you want to plot something, just define here a child method and use it recursively, w/o filling the code with canvases et al. 
 *
 *  @author Riccardo Longo
 *  @bug No known bugs.
 */

#ifndef VISUALIZER_H
#define VISUALIZER_H

#include "TH1.h"
#include "TStyle.h"

#include <iostream>
#include <vector>

class Visualizer {

    public :
    Visualizer( std::string _style );

    //Styles (so fa only ATLAS)
    TStyle* AtlasStyle();
    //Methods to set a give style
    void SetAtlasStyle();

    /** @brief make the argument true to activate the debug mode. False to deactivate it*/
    void SetDebugMode( bool _isDebugon ) { m_debug = _isDebugon; }

    //Special plot treatments
    void   OverlayHistos  ( TH1D *h1, TH1D *h2 , TVirtualPad* pad);

    //Main visualization methods
    void ManyPadsPlot( std:: vector < TH1 > raw_form, std::vector < TH1 > der_form, int nx, int ny, std::string out_name, std::string treatment );


    private :
    /** String identifying the style to be applied */
    std::string m_style;
    /** Debug flag */
    bool m_debug = false;
};

#endif
