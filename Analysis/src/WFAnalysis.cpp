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
    
    for( unsigned int ch = 0; ch < vCh.size(); ch++ ){
      //retrieving information for each channel as a histogram
        TH1D* h = vCh.at(ch)->WF_histo;
        TH1D* hProcessed = vCh.at(ch)->PWF_histo;
        TH1D* hDiff = vCh.at(ch)->FirstDerivative;
        
      //retrieve information as a vector of floats
        std::vector < float > chEntries = vCh.at(ch)->WF;
      
        GetDifferential( h, hDiff );
        vCh.at(ch)->FirstDerivativeRMS = GetRMS( hDiff );
        FindHitWindow( vCh.at(ch) );
        
        if( vCh.at(ch)->was_hit ){
            GetPedestal( vCh.at(ch) );
            vCh.at(ch)->Charge = hProcessed->Integral( vCh.at(ch)->hit_window.first, vCh.at(ch)->hit_window.second );
        }
    }
}

/** @brief GetDifferential method for WFAnalysis
 *  @param hIN Histogram to be differentiated
 *  @param hOUT Output histogram
 *
 * The derivative is calculated as the difference of the summation N points before and after the ith point
 *  @f[
 *       \delta_i(N)= \sum_{k=1}^{N} (s_i + k)  -  \sum_{k=1}^{N} (s_i - k)
 *  @f]
 *
 * N is given by m_diffSens and is set by the constructor or SetDiffSense()
 */
void WFAnalysis::GetDifferential( TH1D *hIN, TH1D *hOUT ){
    
    if(!hIN)  std::cerr <<    "WARNING: Nothing to differentiate"     << std::endl;
    if(!hOUT) std::cerr << "WARNING: Nowhere to put the differential" << std::endl;

    // Loop over histogram
    for( unsigned int bin = m_diffSens; bin < hIN->GetNbinsX() - m_diffSens ; bin++ ){
      //decalare temp variable to store the sum before and after the i data point
      double sum_before = 0;
      double sum_after = 0;
      //calculate the sum before and after the i data point
      for (int i = 0; i < m_diffSens; i++  ){
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
double WFAnalysis::GetRMS( TH1D *h, bool debug){
    
    //Make a histogram with x-range to match the y-range of the input
    Double_t xmin,xmax;
    h->GetMinimumAndMaximum(xmin,xmax);
    if(xmax == 0) return 0;
    
    TH1D hRMS("RMS","RMS",5*(xmax-xmin)/m_diffSens,xmin,xmax);
    
    
    //Loop over the histogram excluding the window used for differentiating to fill hRMS
    Int_t nbins = h->GetNbinsX();
    for(int bin = m_diffSens; bin < nbins - m_diffSens; bin++){
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
 *  @param ch Channel to be processed
 *
 *  Determines if the channel was hit first, then determines the hit window using 
 *  the first derivative and first derivative RMS. Also determines the peak height 
 *  and peak center using the raw waveform value at the first derivative zero crossing. 
 *  Saves the results to Channel members.
 */
void WFAnalysis::FindHitWindow( Channel* ch ){
    
    double threshold = m_Tmultiple*ch->FirstDerivativeRMS;
    ch->Diff_max = ch->FirstDerivative->GetMaximum();
    
    //////////////////////////////////////////////////////////////////
    //  Hit condition likely needs to be refined
    if( ch->Diff_max <= threshold || ch->FirstDerivative->GetMinimum() >= -1*threshold ){
        if( ch->is_on && m_verbose > 1){ std::cerr << std::endl << "No hit found on " << ch->name << std::endl; }
        ch->was_hit = false;
        ch->Peak_center = ch->hit_window.first = ch->hit_window.second = -1;
        ch->Peak_max = 0.0;
        return;
    }
    
    ch->was_hit = true;
    
    int risingEdge = ch->FirstDerivative->GetMaximumBin();
    int fallingEdge = ch->FirstDerivative->GetMinimumBin();
    int nBins = ch->FirstDerivative->GetNbinsX();
    
    // Find the beginning of the hit window
    for(int bin = risingEdge; bin > 0; bin--){
        if(ch->FirstDerivative->GetBinContent(bin) < threshold){
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
        if(ch->FirstDerivative->GetBinContent(bin) > -1*threshold){
            ch->hit_window.second = bin;
            break;
        }
    }
}


/** @brief Finds the pedestal and pedestal rms of the raw waveform
 *  @param ch Channel to be processed
 *
 *  Histograms data from the raw waveform excluding the hit window.
 *  Saves the mean and rms of that histogram as PedMean and PedRMS.
 *  Also saves a pedestal subtracted, zero supressed version of WF_histo to PWF_histo.
 */
void WFAnalysis::GetPedestal( Channel* ch ){
    int nBins = ch->WF_histo->GetNbinsX();
    TH1D h("ped","ped", nBins, 0, nBins);
    
    if( ch->was_hit && ch->hit_window.first > ch->hit_window.second ){
        if(ch->is_on && m_verbose > 1){
            std::cerr << std::endl << "Bad hit window on " << ch->name << ": Cannot get pedestal" << std::endl;
        }
        return;
    }
    // Find PedMean and PedRMS
    for(int bin = 0; bin < nBins; bin++){
        // Skip the hit window
        if( bin == ch->hit_window.first ){ bin = ch->hit_window.second; }
        
        h.Fill( ch->WF_histo->GetBinContent(bin) );
    }
    
    ch->PedMean = h.GetMean();
    ch->PedRMS  = h.GetRMS();
    
    // Use PedMean and PedRMS to generate a processed waveform
    for(int bin = 0; bin < nBins; bin++){
        double content = ch->WF_histo->GetBinContent(bin) - ch->PedMean;
        if( content <= ch->PedRMS ){
            ch->PWF_histo->SetBinContent( bin, 0 );
        }else{
            ch->PWF_histo->SetBinContent( bin, content );
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
