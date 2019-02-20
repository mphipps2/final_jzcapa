/** @file WFAnalysis.cxx
 *  @brief Implementation of WFAnalysis.
 *
 *  Function definitions for WFAnalysis are provided. 
 *  This class is the main  class for the waveform analysis.
 *  It initializes the histograms used for output of processed events.
 *  Also includes methods that accept all waveforms in an event. 
 *  This is where the analysis should be done.
 *
 *  @author Yakov Kulinich
 *  @bug No known bugs.
 */

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TCanvas.h>//Sheng
#include <iostream>

#include "WFAnalysis.h"

/** @brief Default Constructor for WFAnalysis.
 */
WFAnalysis::WFAnalysis( ){

}

/** @brief Destructor for WFAnalysis.
 */
WFAnalysis::~WFAnalysis( ){

}

/** @brief Initialization method for WFAnalysis
 *
 *  Can add other things here that you would 
 *  perhaps not put into the constructor.
 *  I.e. a TTree, some tools. Etc.
 *
 *  @return none
 */
void WFAnalysis::Initialize( ){

}


/** @brief Historgam Setup method for WFAnalysis
 *
 *  Should instantiate any histograms you wish to output here.
 *
 *  @return none
 */
void WFAnalysis::SetupHistograms( ){
  
}


/** @brief Analyze Events method for WFAnalysis
 *
 *  Here is the event-based analysis code.
 *  A const 1D vector of N TH1* is received.
 *  N is the number of channels
 *
 *  @param1 1D vector of all waveforms for N ch 
 *  @return none
 */
void WFAnalysis::AnalyzeEvent( const std::vector< TH1* >& vWFH ){

  // example... you can loop through all histos as follows
  for( unsigned int ch = 0; ch < vWFH.size(); ch++ ){
    TH1* h = vWFH[ch];
    for( unsigned int samp = 0; samp < h->GetNbinsX(); samp++ ){
      // will print what each sample in each channel equals
      // std::cout << ch << " " << samp << " = " << h->GetBinContent( samp + 1 ) << std::endl;
    }
  }
}


/** @brief Analyze Events method for WFAnalysis
 *
 *  Here is the event-based analysis code.
 *  A const 2D vector of NxM is received.
 *  N is the number of channels
 *  M is the number of samples per channel
 *
 *  @param1 cosnt 2D vector for all waveforms in event 
 *  @return none
 */
void WFAnalysis::AnalyzeEvent( const std::vector< std::vector< float > >& vWF ){

  // example... you can loop through each sample in each channel
  // Can do the same in previous AnalyzeEvent( .. ) function, just
  // looking at vWFH[ch]->GetBinContent( samp + 1 );
  for( unsigned int ch = 0; ch < vWF.size(); ch++ ){
    for( unsigned int samp = 0; samp < vWF[ch].size(); samp++ ){
      // will print what each sample in each channel equals
      // std::cout << ch << " " << samp << " = " << vWF[ch][samp] << std::endl;
    }
  }
}

/**
 * @brief Analyze Event method for WF analysis
 *  Receives a  const vector of Channels. Loops over all Channels within
 *  the vector for analysis.
 *  Example of how to retrieve raw data and associated histograms are provided
 *
 * @param1 A vector of pointers to Channel objects
 */
void WFAnalysis::AnalyzeEvent( const std::vector< Channel* > vCh ){

    for( unsigned int ch = 0; ch < vCh.size(); ch++ ){
      //retrieving information for each channel
      TH1D* h = vCh.at(ch)->WF_histo;
      std::vector < float > chEntries = vCh.at(ch)->WF;
      
      TH1D *diff = (TH1D*) GetDifferential( h, 50);
      
      
  }

}

/**
 * @brief Analyze Event method for WF analysis
 *  Receives a  const vector of Channels. Loops over all Channels within
 *  the vector for analysis. Also receives a pointer to a TVirtualPad for graphical output
 *  Example of how to retrieve raw data and associated histograms are provided
 *
 * @param1 A vector of pointers to Channel objects
 * @param2 Pointer to a pad to be drawn to
 */
void WFAnalysis::AnalyzeEvent( std::vector<Channel*> vCh, TVirtualPad* pad ){

    for( unsigned int ch = 0; ch < vCh.size(); ch++ ){
      //retrieving information for each channel
      TH1D* h = vCh.at(ch)->WF_histo;
      std::vector < float > chEntries = vCh.at(ch)->WF;
      
      TH1D *diff = (TH1D*) GetDifferential( h, ch, 50);
      OverlayHistos( h, diff, pad);
      
      
  }

}


/** @brief GetDifferential method for WFAnalysis
 *
 * The derivative is calculated as the difference of the summation N points before and after the ith point
 * i.e new_sample[i] = sum(samp[i+k])-sum(samp[i-k]), k =1, 2, 3 ... sample_range.
 *
 *
 *  @param1 Histogram to be differentiated
 *  @param2 Number of points used for summation window
 *  @param3 Debug boolean, prints samples if true
 *  @return truncated (M-2*window) TH1 histogram 
 */
TH1 *WFAnalysis::GetDifferential( TH1D *h, int N, bool debug){
    TH1D *hNew = new TH1D( Form( "%s Differential", h->GetTitle() ), "diff;samp;derivative", 1024, 0, 1024);

    // Loop over histogram
    for( unsigned int bin = N; bin < h->GetNbinsX() - N ; bin++ ){
      //decalare temp variable to store the sum before and after the i data point
      double sum_before = 0;
      double sum_after = 0;
      //calculate the sum before and after the i data point
      for (int i = 0; i < N; i++  ){
        sum_before += h->GetBinContent( bin - i );
        sum_after  += h->GetBinContent( bin + i );
      }
        //set the bin to the calculated derivative value     
        hNew->SetBinContent(bin,(sum_after - sum_before));
    }//end derivative loop
    
    if (debug){
        TCanvas *c = new TCanvas( Form( "%s Differential", h->GetTitle() ), "Debug", 200, 10, 1000, 600);
        OverlayHistos(h,hNew,c->cd(),true);
    }//end if debug
    return hNew;
}

/** @brief OverlayHistos method for WFAnalysis
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
void WFAnalysis::OverlayHistos( TH1D *h1, TH1D *h2 , TVirtualPad* pad, bool save){
    if( pad == nullptr ) TCanvas pad( h1->GetTitle(), h1->GetTitle(), 200, 10, 1000, 600);
    //Remove Stat box and double the y-axis range to include negative values
    gStyle->SetOptStat( kFALSE );
    pad->cd();
    h1->Draw();
    h1->SetAxisRange( -h1->GetMaximum()*1.1, h1->GetMaximum()*1.1, "Y");
    pad->Update();
    
   //scale h2 to the pad coordinates
   float rightmax = 1.1*h2->GetMaximum();
   float scale = gPad->GetUymax()/rightmax;
   h2->SetLineColor( kRed );
   h2->Scale( scale );
   h2->Draw( "same" );

   //draw an axis on the right side
   TGaxis axis( gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), -rightmax, rightmax, 510, "+L");
   axis.SetLineColor( kRed );
   axis.SetLabelColor( kRed );
   axis.Draw();

   if( save ) pad->Print( Form( "%s_Overlay.pdf", h1->GetTitle() ) ) ;
}

/** @brief OverlayHistos method for WFAnalysis
 *
 *  Generates a new TCanvas and calls OverlayHistos on it.
 *  This allows the user to save the result of OverlayHistos without
 *  providing a pad
 *
 *  @param1 Base histogram (left y-axis)
 *  @param2 Overlayed histogram (right y-axis)
 *  @param3 Save option. If true, save a .pdf
 */
void WFAnalysis::OverlayHistos( TH1D *h1, TH1D *h2 , bool save){
    TCanvas *pad = new TCanvas( h1->GetTitle(), h1->GetTitle(), 200, 10, 1000, 600);
    OverlayHistos( h1, h2, pad, save);
    delete pad;
}

/** @brief Finalize method for WFAnalysis
 *
 *  Write histograms, TTree if it exists.
 *
 *  @return none
 */
void WFAnalysis::Finalize( ){

  // If these exist...
  // m_tree->Write();
  // m_hist->Write();
}
