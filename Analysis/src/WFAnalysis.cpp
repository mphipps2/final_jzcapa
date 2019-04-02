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
 *  @param vWF 1D vector of all waveforms for N ch 
 *
 *  Here is the event-based analysis code.
 *  A const 1D vector of N TH1* is received.
 *  N is the number of channels
 *
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
 *  @param vWF const 2D vector for all waveforms in event 
 *
 *  Here is the event-based analysis code.
 *  A const 2D vector of NxM is received.
 *  N is the number of channels
 *  M is the number of samples per channel
 *
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
 * @param vCh A vector of pointers to Channel objects
 *
 *  Receives a  const vector of Channels. Loops over all Channels within
 *  the vector for analysis.
 *  Example of how to retrieve raw data and associated histograms are provided
 *
 */
void WFAnalysis::AnalyzeEvent( const std::vector< Channel* > vCh ){
    int  ch = 0;
    int  dWindow = 25;
    int  RMSthreshold = 3.5;
    
    for( unsigned int ch = 0; ch < vCh.size(); ch++ ){
      //retrieving information for each channel as a histogram
        TH1D* h = vCh.at(ch)->WF_histo;
        TH1D* hProcessed = vCh.at(ch)->PWF_histo;
        TH1D* hDiff = vCh.at(ch)->FirstDerivative;
        
      //retrieve information as a vector of floats
        std::vector < float > chEntries = vCh.at(ch)->WF;
      
        GetDifferential( h, hDiff, dWindow );
        vCh.at(ch)->FirstDerivativeRMS = GetRMS( hDiff, dWindow );
        FindHitWindow( vCh.at(ch), RMSthreshold);
        
    }
    
}

/** @brief GetDifferential method for WFAnalysis
 *  @param hIN Histogram to be differentiated
 *  @param hOUT Output histogram
 *  @param N Number of points used for summation window 
 *
 * The derivative is calculated as the difference of the summation N points before and after the ith point
 *  @f[
 *       \delta_i(N)= \sum_{k=1}^{N} (s_i + k)  -  \sum_{k=1}^{N} (s_i - k)
 *  @f]
 *
 *
 */
void WFAnalysis::GetDifferential( TH1D *hIN, TH1D *hOUT, int N){
    
    if(!hIN)  std::cerr <<    "WARNING: Nothing to differentiate"     << std::endl;
    if(!hOUT) std::cerr << "WARNING: Nowhere to put the differential" << std::endl;

    // Loop over histogram
    for( unsigned int bin = N; bin < hIN->GetNbinsX() - N ; bin++ ){
      //decalare temp variable to store the sum before and after the i data point
      double sum_before = 0;
      double sum_after = 0;
      //calculate the sum before and after the i data point
      for (int i = 0; i < N; i++  ){
        sum_before += hIN->GetBinContent( bin - i );
        sum_after  += hIN->GetBinContent( bin + i );
      }
        //set the bin to the calculated derivative value     
        hOUT->SetBinContent(bin,(sum_after - sum_before));
        
    }//end derivative loop
}

/** @brief GetRMS method for WFAnalysis
 *  @param h Input histogram
 *  @param diff_window Averaging window used in GetDifferential
 *  @param debug Saves a PDF of the RMS histogram if true
 *  @return Width of gaussian fit excluding tails created by peaks
 *
 *  Given an input histogram, outputs an RMS value
 *  based on a gaussian fit concentrating on the center.
 *  The result can be saved to a PDF if debug is set to true.
 *
 */
double WFAnalysis::GetRMS( TH1D *h , int diff_window, bool debug){
    
    //Make a histogram with x-range to match the y-range of the input
    Double_t xmin,xmax;
    h->GetMinimumAndMaximum(xmin,xmax);
    if(xmax == 0) return 0;
    
    TH1D hRMS("RMS","RMS",5*(xmax-xmin)/diff_window,xmin,xmax);
    
    
    //Loop over the histogram excluding the window used for differentiating to fill hRMS
    Int_t nbins = h->GetNbinsX();
    for(int bin = diff_window; bin < nbins - diff_window; bin++){
        hRMS.Fill( h->GetBinContent( bin ) );
    }
    
    //Make a gaussian fit and apply it to our histogram quietly and using our range
    TF1 f("f","gaus",-250,250);
    hRMS.Fit("f","qR+");
    
    if( debug ){
        TCanvas c( "RMS" , "RMS", 200, 10, 1000, 600);
        c.cd();
        hRMS.Draw();
        f.Draw("same");
        c.Draw();
        c.Print( "RMS.pdf" );
    }
    
    //Return parameter 2. "gaus" is [0]*exp(-0.5*((x-[1])/[2])**2)
    return f.GetParameter(2);
    
}


/** @brief Defines the hit window for a given channel
 *  @param ch
 *
 *  Determines the hit window using the first derivative and first derivative RMS.
 *  Also determines the peak height and peak center using the raw waveform value 
 *  at the first derivative zero crossing. Saves the results to Channel members.
 */
void WFAnalysis::FindHitWindow( Channel* ch, double threshMultiple){
    int risingEdge = ch->FirstDerivative->GetMaximumBin();
    int fallingEdge = ch->FirstDerivative->GetMinimumBin();
    int nBins = ch->FirstDerivative->GetNbinsX();
    
    ch->Diff_max = ch->FirstDerivative->GetMaximum();
    double diffThresh = threshMultiple*ch->FirstDerivativeRMS;
    
    // Find the beginning of the hit window
    for(int bin = risingEdge; bin > 0; bin--){
        if(ch->FirstDerivative->GetBinContent(bin) < diffThresh){
            ch->hit_window.first = bin;
            break;
        }
    }
    
    // Find the peak center using the derivative
    for(int bin = risingEdge; bin < fallingEdge; bin++){
        if(ch->FirstDerivative->GetBinContent(bin) < 0){
            ch->Peak_center = bin;
            ch->Peak_max = ch->WF_histo->GetBinContent(bin) - ch->offset;
            break;
        }
    }
    
    // Find the end of the hit window
    for(int bin = fallingEdge; bin < nBins; bin++){
        if(ch->FirstDerivative->GetBinContent(bin) > -1*diffThresh){
            ch->hit_window.second = bin;
            break;
        }
    }
}

/** @brief Finalize method for WFAnalysis
 *  @return none
 *
 *  Write histograms, TTree if it exists.
 *
 */
void WFAnalysis::Finalize( ){

  // If these exist...
  // m_tree->Write();
  // m_hist->Write();
}
