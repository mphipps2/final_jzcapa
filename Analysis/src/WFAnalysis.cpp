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

/** @brief GetDifferential method for WFAnalysis
 *
 *  h const 1D vector of M dimenstion is received, M is the number of samples per channel
 *  ch is the channel number.
 *  window is the integral window for the pre-integration
 *  debug is to draw the first 5 raw signals and the derivative 
 *
 *  @return truncated (M-2*window) 1D TH1 histogram 
 */
TH1 *WFAnalysis::GetDifferential( const TH1 *h, unsigned int ch, int window, bool debug){
    TH1F *hNew = new TH1F( Form( "ch_%d", ch ), "diff;samp;derivative", 1004, 0 + window, 1024 - window );
    int new_sample = window;

    // sample loop
    for( unsigned int samp = window; samp < h->GetNbinsX() - window ; samp++ ){
      //decalare temp variable to store the sum before and after the i data point
      double sum_previous = 0;
      double sum_after = 0;
      //calculate the sum before and after the i data point
      for (int i = 0; i < window; i++  ){
        sum_previous = (h->GetBinContent( samp - i ) + sum_previous);
        sum_after = (h->GetBinContent( samp + i ) + sum_after);
      }
      //set the difference of two sum to the new histogram        
        hNew->SetBinContent(new_sample,(sum_after - sum_previous));
      new_sample++;
    }// end of sample loop
    if (debug){
      if (ch <=5){
        TCanvas *c2 = new TCanvas(Form( "ch_%d", ch ),"Canvas debug",200,10,1000,600);
        c2->Divide(1,2);
        c2->cd(1);
        h->DrawCopy();
        c2->cd(2);
        hNew->DrawCopy();
        c2->Print(Form( "ch_%d.pdf", ch ));
      }
    }//if debug end
    return hNew;
  }

/**
 * @brief Analyze Event method for WF analysis
 *  A const vector of channels is received. Can be either the vector coming from a single detector or formed using all the channels.
 *  Example of how to retrieve raw data and associated histograms are provided
 * @param vCh
 */
void WFAnalysis::AnalyzeEvent( const std::vector<Channel *> vCh ){
    /* The derivative is achieve by taking the difference of the summation N points before and after the i channel
    i.e new_sample[i] = sum(samp[i+k])-sum(samp[i-k]), k =1, 2, 3 ... sample_range.
    */
    //Create a Canvas to output raw signal and first derivative of the raw signal
    TCanvas *c1 = new TCanvas("c2","Canvas ped",200,10,1000,600);
    c1->Divide(1,2);
    // create a hitogram to store the derivative of the signal
    // discard the first and last sample_range points in the raw signal
     
    //set up the sample range for the derivative
    unsigned int sample_range = 50;
    //in order to align with the raw signal,the new histogram start at the 50th point of the raw signal
    int new_sample = sample_range;

    for( unsigned int ch = 0; ch < vCh.size(); ch++ ){
      //retrieving information for each channel
      TH1* h = vCh.at(ch)->WF_histo;
      std::vector < float > chEntries = vCh.at(ch)->WF;
      

      GetDifferential(h,ch,sample_range);

  }

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
