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
#include "TH1.h"
#include "TCanvas.h"
#include "Visualizer.h"

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
  
    //TH1D
    hChargeRatio = new TH1D("ChargeZDC2_over_ChargeZDC1","Q_{ZDC2}/Q_{ZDC1}",100,0,1.);
    hChargeRatio->SetCanExtend(TH1::kXaxis);

    hPeakRatio = new TH1D("PeakZDC2_over_PeakZDC1","Peak_{ZDC2}/Peak_{ZDC1}",100,0,1.);
    hPeakRatio->SetCanExtend(TH1::kXaxis);

    //TH2D
    hCharge = new TH2D("ZDC_Charge_Correlation", "ZDC Charge Correlation", 100, 0, 7000, 100, 0, 7000);
    hCharge->SetCanExtend(TH1::kAllAxes);

    hPeak   = new TH2D("ZDC_Peak_Correlation", "ZDC Peak Correlation", 100, 0, 350, 100, 0, 350);
    hPeak->SetCanExtend(TH1::kAllAxes);

    hChargePeakZDC1 = new TH2D("ZDC1_ChargePeakCorrelation","Q_{ZDC1} vs Peak_{ZDC1}",50,0,10000,50,0,1000);
    //hChargePeakZDC1->SetCanExtend(TH1::kAllAxes);

    hChargePeakZDC2 = new TH2D("ZDC2_ChargePeakCorrelation","Q_{ZDC2} vs Peak_{ZDC2}",50,0,10000,50,0,1000);
    //hChargePeakZDC2->SetCanExtend(TH1::kAllAxes);

}


/** @brief Analyze Events method for ZDCAnalysis
 *
 *
 */
void ZDCAnalysis::AnalyzeEvent( ){

    
    zdc1 = m_zdc1->GetElement(0,1);
    zdc2 = m_zdc2->GetElement(0,2);
    
    if(zdc1 && zdc2){

        hChargeRatio->Fill(zdc2->Charge/zdc1->Charge);
        hPeakRatio->Fill(zdc2->Peak_max/zdc1->Peak_max);

        hCharge->Fill( zdc1->Charge,   zdc2->Charge );
        hPeak->Fill(   zdc1->Peak_max, zdc2->Peak_max );

        hChargePeakZDC1->Fill(zdc1->Charge,zdc1->Peak_max);
        hChargePeakZDC2->Fill(zdc2->Charge,zdc2->Peak_max);
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

    Visualizer* viz = new Visualizer( "ATLAS" );

    viz->DrawPlot(hChargeRatio,"Q_{ZDC2}/Q_{ZDC1}","Counts","ZDC_chargeRatio.png","");
    viz->DrawPlot(hChargeRatio,"Peak_{ZDC2}/Peak_{ZDC1}","Counts","ZDC_peakRatio.png","");

    viz->DrawPlot(hCharge,"Q_{ZDC1}","Q_{ZDC2}","ZDC_charge.png","COLZ");
    viz->DrawPlot(hPeak,"Peak_{ZDC1}","Peak_{ZDC2}","ZDC_peak.png","COLZ");
    viz->DrawPlot(hChargePeakZDC1,"Q_{ZDC1}","Peak_{ZDC1}","ZDC1_ChargePeak.png","COLZ");
    viz->DrawPlot(hChargePeakZDC2,"Q_{ZDC2}","Peak_{ZDC2}","ZDC2_ChargePeak.png","COLZ");

    delete c;
    
}



























