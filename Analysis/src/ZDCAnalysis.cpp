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


/** @brief Historgam Setup method for ZDCAnalysis
 *
 *  Should instantiate any histograms you wish to output here.
 *
 *  @return none
 */
void ZDCAnalysis::SetupHistograms( ){
  
    hCharge = new TH2D("ZDC Charge Correlation", "ZDC Charge Correlation", 1000, 0, 15000, 1000, 0, 15000);
    hPeak   = new TH2D("ZDC Peak Correlation", "ZDC Peak Correlation", 1000, 0, 500, 1000, 0, 500);
  
}


/** @brief Analyze Events method for ZDCAnalysis
 *
 *
 */
void ZDCAnalysis::AnalyzeEvent( std::vector < Detector* > _vDet ){

    // If the detectors haven't been filled, fill them
    if( m_zdc1 == 0 || m_zdc2 == 0 ){
        for( auto& det : _vDet ){
            if(det->GetChannelsVector().size() > 1 ){continue;}
            if(det->GetElement(0,1) != nullptr ){ 
                m_zdc1 = (ZDC*)det; 
                std::cout << "ZDC1 assigned row " << det->GetElement(0,1)->mapping_row << ", column " << det->GetElement(0,1)->mapping_column << std::endl;
            }
            if(det->GetElement(0,2) != nullptr){ 
                m_zdc2 = (ZDC*)det; 
                std::cout << "ZDC2 assigned row " << det->GetElement(0,2)->mapping_row << ", column " << det->GetElement(0,2)->mapping_column << std::endl;
            }
        }
        std::cout << " ^^^^ Don't worry about those warnings ^^^^" << std::endl;
        std::cout << "      Unless there are wanings below...  " << std::endl;
        
    }
    
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
    
    TCanvas *c = new TCanvas("ZDCAnalysis","ZDCAnalysis",800,600);
    c->cd();
    hCharge->Draw("COLZ");
    c->Print("output/ZDC_charge.png");
    hPeak->Draw("COLZ");
    c->Print("output/ZDC_peak.png");
    delete c;
    
    
}



























