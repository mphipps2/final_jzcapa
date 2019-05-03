/** @file ZDCAnalysis.cxx
 *  @brief Implementation of ZDCAnalysis.
 *
 *  
 *  Function definitions for ZDCAnalysis are provided. 
 *  This class is the main class for the ZDC analysis.
 *  The analysis is done on both ZDCs simultaneously
 *  
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
*/
 
 
#include "ZDCAnalysis.h"
#include "ZDC.h"
#include "Containers.h"
#include "TH2D.h"
#include "TCanvas.h"


/** @brief Default Constructor for ZDCAnalysis.
 */
ZDCAnalysis::ZDCAnalysis( ){

}

/** @brief Destructor for ZDCAnalysis.
 */
ZDCAnalysis::~ZDCAnalysis( ){

}

/** @brief Initialization method for ZDCAnalysis
 *
 *  Can add other things here that you would 
 *  perhaps not put into the constructor.
 *  I.e. a TTree, some tools. Etc.
 *
 */
void ZDCAnalysis::Initialize( ){

}

/** @brief Initialization method for ZDCAnalysis
 *
 *  Takes a vector of detectors, picks out the ZDCs
 *  and assigns them to member pointers
 *
 */
void ZDCAnalysis::Initialize( std::vector < Detector* > _vDet ){

    for( auto& det : _vDet ){
        if( det->GetChannelsVector().at(0)->detector == "ZDC" && det->GetChannelsVector()[0]->mapping_column == 1){
            m_zdc1 = (ZDC*)det;
            zdc1 = det->GetElement(0,1);
        }
        if( det->GetChannelsVector().at(0)->detector == "ZDC" && det->GetChannelsVector()[0]->mapping_column == 2){
            m_zdc2 = (ZDC*)det;
            zdc2 = det->GetElement(0,2);
        }
    }
}


/** @brief Historgam Setup method for ZDCAnalysis
 *
 *  Should instantiate any histograms you wish to output here.
 *
 *  @return none
 */
void ZDCAnalysis::SetupHistograms( ){
  
    hCharge = new TH2D("ZDC Charge Correlation", "ZDC Charge Correlation", 100, 0, 15000, 100, 0, 9500);
    hPeak   = new TH2D("ZDC Peak Correlation", "ZDC Peak Correlation", 100, 0, 1200, 100, 0, 1200);
  
}


/** @brief Analyze Events method for ZDCAnalysis
 *
 *
 */
void ZDCAnalysis::AnalyzeEvent( ){

    
    zdc1 = m_zdc1->GetElement(0,1);
    zdc2 = m_zdc2->GetElement(0,2);
    
    if(zdc1 && zdc2){
        hCharge->Fill( zdc1->Charge,   zdc2->Charge );
        hPeak->Fill(   zdc1->Peak_max, zdc2->Peak_max );
    }
    
}

/** @brief Finalize method for ZDCAnalysis
 *
 *
 */
void ZDCAnalysis::Finalize( ){
    
    std::string output =  std::getenv("JZCaPA");
    output += "/results/";
    
    TCanvas *c = new TCanvas("ZDCAnalysis","ZDCAnalysis",800,600);
    c->cd();
    hCharge->Draw("COLZ");
    c->Print( (output + "ZDC_charge.png").c_str() );
    hPeak->Draw("COLZ");
    c->Print( (output + "ZDC_peak.png").c_str() );
    delete c;
    
    
}



























