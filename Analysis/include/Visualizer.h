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
    //Methods to set a given style
    void SetAtlasStyle();

    /** @brief allow the user to define the extension of the plots once they're printed. ".pdf" by default */
    void SetPlotExtension( std::string _extension ) { m_extension = _extension; }

    /** @brief make the argument true to activate the debug mode. False to deactivate it*/
    void   SetDebugMode   ( bool _isDebugon ) { m_debug = _isDebugon; }

    //Special plot treatments
    void   OverlayHistos  ( TH1 *h1, TH1 *h2 , TVirtualPad* pad);

    //Main visualization methods
    void   ManyPadsPlot   ( std::vector< TH1* > _first_form, std::vector< TH1* > _second_form, int _ncol, int _nrow, std::string _out_name, std::string _treatment );


    private :
    /** String identifying the style to be applied */
    std::string m_style;
    /** Debug flag */
    bool m_debug = false;
    /** String defining the extension used to print the plots */
    std::string m_extension = ".pdf";
};

#endif
