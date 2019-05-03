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
            
            for(int row = 1; row <= 4; row++){
                for(int col = 1; col <= 4; col++){
                    rpd[row][col] = m_RPD->GetElement(row,col);
                }
            }
        }
    }
}

/** @brief Historgam Setup method for RPDAnalysis
 *
 *  Should instantiate any histograms you wish to output here.
 *
 *  @return none
 */
void RPDAnalysis::SetupHistograms( ){
  
   hCharge = new TH2D("RPD Charge", "average charge per tile", 4, 0, 4, 4, 0, 4);
   hPeak   = new TH2D("RPD Peak", "average peak height per tile", 4, 1, 5, 4, 1, 5);
   hCenter = new TH2D("Shower Center", "Calculated center of mass", 200, -50, 50, 200, -50, 50);
  
}


/** @brief Analyze Events method for RPDAnalysis
 *
 *
 */
void RPDAnalysis::AnalyzeEvent( ){
    
    
    //This loop takes the running average of charge and peak height per tile
    for(int row = 1; row <= 4; row++){
        for(int col = 1; col <= 4; col++){
            if( rpd[row][col]->is_on && rpd[row][col]->was_hit ){
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
    double xCM = weightedRow/totalCharge;
    double yCM = weightedCol/totalCharge;
    hCenter->Fill(xCM,yCM);


 
}

/** @brief Finalize method for RPDAnalysis
 *
 *
 */
void RPDAnalysis::Finalize( ){
    
    TCanvas c("RPD","RPD",800,600);
    for(int row = 1; row <= 4; row++){
        for(int col = 1; col <= 4; col++){
            hCharge->SetBinContent(5-col, 5-row, m_charge[row][col]);
            hPeak->SetBinContent(  5-col, 5-row, m_peak[row][col]  );
        }
    }
    
    std::string output =  std::getenv("JZCaPA");
    output += "/results/";
    
    c.cd();
    hCharge->Draw("COLZ text");
    c.Print( (output + "RPD_charge.C").c_str() );
    c.Print( (output + "RPD_charge.png").c_str() );
    hPeak->Draw("COLZ text");
    c.Print( (output + "RPD_peak.png").c_str() );
    c.Print( (output + "RPD_peak.C").c_str() );
    hCenter->Draw("COLZ");
    c.Print( (output + "RPD_CoM.png").c_str() );
    c.Print( (output + "RPD_Com.C").c_str() );
    
}



























