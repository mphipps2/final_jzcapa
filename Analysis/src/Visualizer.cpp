/** @file Visualizer.cpp
 *  @brief Implementation of Visualizer class.
 *
 *
 *  @author Riccardo Longo
 *  @bug No known bugs.
 */

#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TPad.h>

#include <iostream>

#include "Visualizer.h"

Visualizer::Visualizer( std::string _style ){
    m_style = _style;
    if(m_style == "ATLAS" || m_style == "atlas") Visualizer::SetAtlasStyle();
}

/**
 * @brief Function inherited from rcdaq. Creates the ATLAS TStyle object
 * @return
 */
TStyle* Visualizer::AtlasStyle(){
    TStyle *atlasStyle = new TStyle("ATLAS","Atlas style");

      // use plain black on white colors
      Int_t icol=0; // WHITE
      atlasStyle->SetFrameBorderMode(icol);
      atlasStyle->SetFrameFillColor(icol);
      atlasStyle->SetCanvasBorderMode(icol);
      atlasStyle->SetCanvasColor(icol);
      atlasStyle->SetPadBorderMode(icol);
      atlasStyle->SetPadColor(icol);
      atlasStyle->SetStatColor(icol);
      //atlasStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

      // set the paper & margin sizes
      atlasStyle->SetPaperSize(20,26);

      // set margin sizes
      atlasStyle->SetPadTopMargin(0.05);
      atlasStyle->SetPadRightMargin(0.05);
      atlasStyle->SetPadBottomMargin(0.16);
      atlasStyle->SetPadLeftMargin(0.16);

      // set title offsets (for axis label)
      atlasStyle->SetTitleXOffset(1.4);
      atlasStyle->SetTitleYOffset(1.4);

      // use large fonts
      //Int_t font=72; // Helvetica italics
      Int_t font=42; // Helvetica
      Double_t tsize=0.05;
      atlasStyle->SetTextFont(font);

      atlasStyle->SetTextSize(tsize);
      atlasStyle->SetLabelFont(font,"x");
      atlasStyle->SetTitleFont(font,"x");
      atlasStyle->SetLabelFont(font,"y");
      atlasStyle->SetTitleFont(font,"y");
      atlasStyle->SetLabelFont(font,"z");
      atlasStyle->SetTitleFont(font,"z");

      atlasStyle->SetLabelSize(tsize,"x");
      atlasStyle->SetTitleSize(tsize,"x");
      atlasStyle->SetLabelSize(tsize,"y");
      atlasStyle->SetTitleSize(tsize,"y");
      atlasStyle->SetLabelSize(tsize,"z");
      atlasStyle->SetTitleSize(tsize,"z");

      // use bold lines and markers
      atlasStyle->SetMarkerStyle(20);
      atlasStyle->SetMarkerSize(1.2);
      atlasStyle->SetHistLineWidth(2.);
      atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

      // get rid of X error bars
      //atlasStyle->SetErrorX(0.001);
      // get rid of error bar caps
      atlasStyle->SetEndErrorSize(0.);

      // do not display any of the standard histogram decorations
      atlasStyle->SetOptTitle(0);
      //atlasStyle->SetOptStat(1111);
      atlasStyle->SetOptStat(0);
      //atlasStyle->SetOptFit(1111);
      atlasStyle->SetOptFit(0);

      // put tick marks on top and RHS of plots
      atlasStyle->SetPadTickX(1);
      atlasStyle->SetPadTickY(1);

      return atlasStyle;

}

/**
 * @brief Function inherited from rcdaq. Set the ATLAS style.
 */
void Visualizer::SetAtlasStyle(){

    static TStyle* atlasStyle = 0;
    std::cout << "\nApplying ATLAS style settings...\n" << std::endl ;
    if ( atlasStyle==0 ) atlasStyle = AtlasStyle();
    gROOT->SetStyle("ATLAS");
    gROOT->ForceStyle();

}


/** @brief OverlayHistos method for Visualizer
 *
 *  Plots two input histograms on the same, specified pad with
 *  separate axis. Saves the plots as PDFs with the name of the
 *  base histogram if given the option.
 *
 *  @param1 Base histogram (left y-axis)
 *  @param2 Overlayed histogram (right y-axis)
 *  @param3 Address of a pad (TPad or TCanvas) to be drawn on
 *  @param4 Save option. If true, save a .pdf
 */
void Visualizer::OverlayHistos( TH1D *h1, TH1D *h2 , TVirtualPad* pad){

    // If there is no pad or no data in the histograms, return
    if( pad == nullptr ) {std::cerr<< "WARNING: No pad to overlay histos onto" << std::endl; return;}
    if( !h1->GetMinimum() && !h1->GetMaximum()) {std::cerr << "WARNING: "<< h1->GetTitle() << " is empty. Can't overlay" << std::endl; return;}

    //Remove Stat box and double the y-axis range to include negative values
    gStyle->SetOptStat( kFALSE );
    h1->Draw();
    h1->SetAxisRange( -h1->GetMaximum()*1.1, h1->GetMaximum()*1.1, "Y");
    pad->Update();

   //scale h2 to the pad coordinates
   float rightmax = 1.1*h2->GetMaximum();
   float scale = gPad->GetUymax()/rightmax;
   h2->SetLineColor( kRed );
   h2->Scale( scale );
   h2->DrawCopy( "same" );

   //draw an axis on the right side
   TGaxis axis( gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), -rightmax, rightmax, 510, "+L");
   axis.SetLineColor( kRed );
   axis.SetLabelColor( kRed );
   axis.DrawClone();

   //Rescale h2 back to the original size
   h2->Scale( 1/scale );
   if( m_debug ) pad->Print( Form( "%s_Overlay.pdf", h1->GetTitle() ) ) ;
}

