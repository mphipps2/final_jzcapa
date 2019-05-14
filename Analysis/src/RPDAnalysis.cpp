/** @file RPDAnalysis.cxx
 *  @brief Implementation of RPDAnalysis.
 *
 *  
 *  Function definitions for RPDAnalysis are provided. 
 *  This class is the main class for the RPD analysis.
 *  The analysis is done on both RPDs simultaneously
 *  
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
*/
 
 
#include "RPDAnalysis.h"
#include "RPD.h"
#include "Containers.h"
#include "TMarker.h"
#include "TH2D.h"
#include "TCanvas.h"


/** @brief Default Constructor for RPDAnalysis.
 */
RPDAnalysis::RPDAnalysis( ){

}

/** @brief Destructor for RPDAnalysis.
 */
RPDAnalysis::~RPDAnalysis( ){

}

/** @brief Initialization method for ZDCAnalysis
 *
 *  Takes a vector of detectors, picks out the RPD
 *  and channels, then assigns them to member pointers
 *
 */
void RPDAnalysis::Initialize( std::vector < Detector* > _vDet ){

    for( auto& det : _vDet ){
        if(det->GetChannelsVector().at(0)->detector == "RPD" ){
            m_RPD = (RPD*) det;
            m_alignment = m_RPD->GetAlignment();
            for(int row = 1; row <= 4; row++){
                for(int col = 1; col <= 4; col++){
                    rpd[row][col] = m_RPD->GetElement(row,col);
                }
            }
        }
    }
    beamPosX = -(m_alignment->x_table + 2250)/100;
    beamPosY = -(m_alignment->y_table - 500)/100;
}

/** @brief Historgam Setup method for RPDAnalysis
 *
 *  Should instantiate any histograms you wish to output here.
 *
 *  @return none
 */
void RPDAnalysis::SetupHistograms( ){
  
   hCharge     = new TH2D("RPD Charge", "average charge per tile", 4, 0, 4, 4, 0, 4);
   hPeak       = new TH2D("RPD Peak", "average peak height per tile", 4, 1, 5, 4, 1, 5);
   hCenter     = new TH2D("Shower Center", "Calculated center of mass", 200, -50, 50, 200, -50, 50);
   hChgVsPk    = new TH2D("Charge vs Peak height","Q vs Peak; Q; Peak", 300, 0, 10000, 300, 0, 1100);
   hPkVsDiffPk = new TH2D("Peak height vs Differential peak height","Peak vs Diff Peak", 200, 0, 5000, 200, 0, 300);
   hChargeSum  = new TH1D("RPD Integrated Signal","RPD Integrated Signal", 200, 0, 160000);
   hPeakSum    = new TH1D("RPD Peak Height Sum","RPD Peak Height Sum", 200, 0, 3000);
   hDiffPeakSum= new TH1D("RPD Differential Peak Height Sum","RPD Diff Peak Sum", 200, 0, 55000);
   
   for(int row = 1; row <=4; row++){
       for(int col = 1; col <= 4; col++){
            hChargeArr[row][col] = new TH1D(Form("rpd%d_%d_Charge",row,col),Form("rpd%d_%d_Charge",row,col), 200 , 0, 15000);
            hPeakArr[row][col]   = new TH1D(Form("rpd%d_%d_Peak",row,col),Form("rpd%d_%d_Peak",row,col), 200, 0, 300);
            hDPeakArr[row][col]  = new TH1D(Form("rpd%d_%d_Diff_peak",row,col),Form("rpd%d_%d_Diff_peak",row,col), 200, 0, 5000);
       }
   }
  
}

/** @brief Branch setup method for ZDCAnalysis
 *
 *  Adds branches with data created by the analysis
 *
 */
void RPDAnalysis::SetBranches( TTree* _tree ){
    m_AnalysisTree = _tree;
    
    for(int row = 1; row <= 4; row++){
        for(int col = 1; col <= 4; col++){
            m_AnalysisTree->Branch( Form("rpd%d_%d_Charge", row, col),           &rpd[row][col]->Charge,           Form("rpd%d_%d_Charge/D", row, col) );
            m_AnalysisTree->Branch( Form("rpd%d_%d_Peak_max", row, col),         &rpd[row][col]->Peak_max,         Form("rpd%d_%d_Peak_max/D", row, col) );
            m_AnalysisTree->Branch( Form("rpd%d_%d_Diff_max", row, col),         &rpd[row][col]->Diff_max,         Form("rpd%d_%d_Diff_max/D", row, col) );
            m_AnalysisTree->Branch( Form("rpd%d_%d_Peak_center", row, col),      &rpd[row][col]->Peak_center,      Form("rpd%d_%d_Peak_center/I", row, col) );
            m_AnalysisTree->Branch( Form("rpd%d_%d_Diff_Peak_center", row, col), &rpd[row][col]->Diff_Peak_center, Form("rpd%d_%d_Diff_Peak_center/I", row, col) );
        }
    }
   
    m_AnalysisTree->Branch("rpd_xCoM", &xCoM, "xCoM/D" );
    m_AnalysisTree->Branch("rpd_yCoM", &yCoM, "yCoM/D" );
    m_AnalysisTree->Branch("rpd_Charge_sum", &ChargeSum, "ChargeSum/D" );
    m_AnalysisTree->Branch("rpd_Peak_sum",   &PeakSum,   "PeakSum/D" );
    m_AnalysisTree->Branch("rpd_Diff_Peak_sum",   &DiffPeakSum,   "DiffPeakSum/D" );
    
}

/** @brief Analyze Events method for RPDAnalysis
 *
 *
 */
void RPDAnalysis::AnalyzeEvent( ){
    ChargeSum = 0;
    PeakSum = 0;
    DiffPeakSum = 0;
    
    //This loop takes the running average of charge and peak height per tile
    for(int row = 1; row <= 4; row++){
        for(int col = 1; col <= 4; col++){
            if( rpd[row][col]->is_on && rpd[row][col]->was_hit ){
                ChargeSum   += rpd[row][col]->Charge;
                PeakSum     += rpd[row][col]->Peak_max;
                DiffPeakSum += rpd[row][col]->Diff_max;
                
                hChargeArr[row][col]->Fill(rpd[row][col]->Charge);
                hPeakArr[row][col]->Fill(rpd[row][col]->Peak_max);
                hDPeakArr[row][col]->Fill(rpd[row][col]->Diff_max);
                
                m_charge[row][col] = (m_charge[row][col] + rpd[row][col]->Charge  )/2;
                m_peak[row][col] =   (m_peak[row][col]   + rpd[row][col]->Peak_max)/2;
            }
        }
    }
 
    //This loop finds the center of mass per event (should probably be it's own function)
    
    double totalCharge = 0;
    double weightedRow = 0;
    double weightedCol = 0;
    for(int i = 1; i <= 4; i++){
        for(int j = 1; j <= 4; j++){
            if( rpd[i][j]->is_on && rpd[i][j]->was_hit ){
                totalCharge += m_charge[i][j];
                
                weightedCol += m_charge[j][i] * yPos[j];
                weightedRow += m_charge[i][j] * xPos[j];
            }
        }
    }
    xCoM = weightedRow/totalCharge;
    yCoM = weightedCol/totalCharge;
    hCenter->Fill(xCoM,yCoM);
    
    hChargeSum->Fill(ChargeSum);
    hPeakSum->Fill(PeakSum);
    hDiffPeakSum->Fill(DiffPeakSum);


 
}

/** @brief Finalize method for RPDAnalysis
 *
 *
 */
void RPDAnalysis::Finalize( ){

    if(m_viz == NULL) m_viz = new Visualizer( "ATLAS" );
    
    TCanvas cCharge("charge","charge",800,600);
    TCanvas cPeak("peak","peak",800,600);
    TCanvas cDpeak("dpeak","dpeak",800,600);
    cCharge.Divide(4,4);
    cPeak.Divide(4,4);
    cDpeak.Divide(4,4);
    int pad;
    
    for(int row = 1; row <= 4; row++){
        for(int col = 1; col <= 4; col++){
            pad = row*4 - (col - 1);
            
            cCharge.cd(pad);
            hChargeArr[row][col]->Draw();
            
            cPeak.cd(pad);
            hPeakArr[row][col]->Draw();
            
            cDpeak.cd(pad);
            hDPeakArr[row][col]->Draw();
            
            hCharge->SetBinContent(5-col, 5-row, m_charge[row][col]);
            hPeak->SetBinContent(  5-col, 5-row, m_peak[row][col]  );
        }
    }
    
    std::string output =  std::getenv("JZCaPA");
    output += "/results/";
    
    TMarker *marker = new TMarker(beamPosX,beamPosY,5);
    marker->SetMarkerColor(kRed);
    marker->SetMarkerSize(4);
    
    cCharge.Print("RPD_Charge_per_tile.png");
    cPeak.Print("RPD_Peak_per_tile.png");
    cDpeak.Print("RPD_Diff_peak_per_tile.png");
    
    m_viz->DrawPlot(hChargeSum,"Q","Counts","RPD_TotalCharge.png","");
    m_viz->DrawPlot(hPeakSum,"Peak","Counts","RPD_PeakSum.png","");
    m_viz->DrawPlot(hDiffPeakSum,"$\frac{#partial V}{#partial t}_{max}$","Counts","RPD_DiffSum.png","");
    m_viz->DrawPlot(hCharge,"Row","Col","RPD_Charge.png","COLZ",marker);
    m_viz->DrawPlot(hPeak,"Row","Col","RPD_Peak.png","COLZ",marker);
    m_viz->DrawPlot(hCenter,"x (mm)","y (mm)","RPD_CoM.png","COLZ",marker);
    
}



























