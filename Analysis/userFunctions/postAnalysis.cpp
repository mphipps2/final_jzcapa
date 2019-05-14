/** @file postAnalysis.cpp
 *  @brief Example of an analysis of the output of the primary JZCaPA analysis
 *
 *  Takes run number as an argument and searches $JZCaPA/results directory for the data.
 *  This is not for the analysis of raw data, but instead is for faster analyis of pre-processed data.
 *  
 *
 *
 *  @author Chad Lantz
 *  @bug No known bugs.
*/

#include <stdlib.h>
#include "Visualizer.h"
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TImage.h>
#include <TText.h>
#include <TASImage.h>
#include <TH2.h>


using namespace std;


int main(int argc, char *argv[]){
    
    //Open the file and exit if it fails
    int runNo = atoi(argv[1]);
    string Path = (string)getenv("JZCaPA") + "/results/";
    string file = Path + Form("output%d.root", runNo);
    
    TFile *f = new TFile(file.c_str());
    if(f->IsZombie()){
        cerr << "File didn't open... closing" << endl;
        exit(-1);
    }
    TTree *t=(TTree*)f->Get("AnalysisTree");
    
    //Make variables to hold data
    double zdc1_Charge, zdc1_Peak_max, zdc2_Charge, zdc2_Peak_max, rpd_xCoM, rpd_yCoM, rpd_Charge_sum, rpd_Peak_sum, rpd_Charge[5][5], rpd_Peak_max[5][5], rpd_Diff_max[5][5];
    int runNumber, evNo, zdc1_Peak_center, zdc1_Diff_Peak_center, zdc2_Peak_center, zdc2_Diff_Peak_center, rpd_Peak_center[5][5], rpd_Diff_Peak_center[5][5];
    
    //Set branch addresses
    t->SetBranchAddress("evNo", &evNo);
    t->SetBranchAddress("zdc1_Charge", &zdc1_Charge);
    t->SetBranchAddress("zdc1_Peak_max", &zdc1_Peak_max);
    t->SetBranchAddress("zdc1_Peak_center", &zdc1_Peak_center);
    t->SetBranchAddress("zdc1_Diff_Peak_center", &zdc1_Diff_Peak_center);
    t->SetBranchAddress("zdc2_Charge", &zdc2_Charge);
    t->SetBranchAddress("zdc2_Peak_max", &zdc2_Peak_max);
    t->SetBranchAddress("zdc2_Peak_center", &zdc2_Peak_center);
    t->SetBranchAddress("zdc2_Diff_Peak_center", &zdc2_Diff_Peak_center);
    t->SetBranchAddress("rpd_xCoM", &rpd_xCoM);
    t->SetBranchAddress("rpd_yCoM", &rpd_yCoM);
    t->SetBranchAddress("rpd_Charge_sum", &rpd_Charge_sum);
    t->SetBranchAddress("rpd_Peak_sum", &rpd_Peak_sum);
    for(int row = 1; row <= 4; row++){
        for(int col = 1; col <=4; col++){
            t->SetBranchAddress(Form("rpd%d_%d_Charge",row,col), &rpd_Charge[row][col]);
            t->SetBranchAddress(Form("rpd%d_%d_Peak_max",row,col), &rpd_Peak_max[row][col]);
            t->SetBranchAddress(Form("rpd%d_%d_Diff_max",row,col), &rpd_Diff_max[row][col]);
            t->SetBranchAddress(Form("rpd%d_%d_Peak_center",row,col), &rpd_Peak_center[row][col]);
            t->SetBranchAddress(Form("rpd%d_%d_Diff_Peak_center",row,col), &rpd_Diff_Peak_center[row][col]);
        }
    }
    
    //Check RPD center of mass for time dependence. Saves a gif of time slices
    int nEntries = t->GetEntries();
    int nEntriesPerFrame = 100;
    
    TImage *img = 0;
    gSystem->Unlink("anim1.gif");  // delete existing file
    TCanvas *c = new TCanvas("blah","blah",800,600);
    
    TPad *c1_1 = new TPad("c1_1", "c1_1",0.01,0.01,0.75,0.99); //Left half
    TPad *c1_2 = new TPad("c1_2", "c1_2",0.75,0.5,0.99,0.99); // Top right
    TPad *c1_3 = new TPad("c1_3", "c1_3",0.75,0.01,0.99,0.5);  // Bottom right
    
    TH2D *hCoM = new TH2D("RPD_CoM","RPD Center of Mass;x (mm);y(mm)", 200, -50, 50, 200, -50, 50);
    hCoM->SetAxisRange(0,15,"Z");
    //hCoM->SetStats(false);
    TH1D *hxCoM = new TH1D("RPD_xCoM","RPD Center of Mass x;x (mm)", 50, -25, 25 );
    TH1D *hyCoM = new TH1D("RPD_yCoM","RPD Center of Mass y;y (mm)", 50, -25, 25 );
    
    TText *text = new TText();
    
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.9);
    
    for(int ev = 0; ev < nEntries; ev++){
        t->GetEntry(ev);
        hCoM->Fill(rpd_xCoM,rpd_yCoM);
        hxCoM->Fill(rpd_xCoM);
        hyCoM->Fill(rpd_yCoM);
        
        if(ev%nEntriesPerFrame == 0){
            delete img;
            
            //Draw the figures
            c->cd();
            c1_1->Draw();
            c1_1->cd();
            hCoM->Draw("COLZ");
            text->DrawTextNDC(0.25,0.75, Form("Event %d",ev) );
            c->cd();
            c1_2->Draw();
            c1_2->cd();
            hxCoM->Draw();
            c->cd();
            c1_3->Draw();
            c1_3->cd();
            hyCoM->Draw();
            
            //Print the canvas and open the file with TImage
            c->Print( (Path + "CoM.png").c_str() );
            img = TImage::Open( (Path + "CoM.png").c_str() );
            
            //Add the image to the back of the stack. If it's the last image, set it to loop.
            if( ev < nEntries-nEntriesPerFrame ){
                img->WriteImage( (Path + "anim1.gif+10").c_str() );
            }else{// the last iamge written. "++" stands for infinite loop (Doesn't work for some reason)
                img->WriteImage( (Path + "anim1.gif++10++").c_str() );
            }
            
            //Reset all histograms
            hCoM->Reset();
            hxCoM->Reset();
            hyCoM->Reset();
        }
    }
    
    delete c;
    delete hCoM;
    delete hxCoM;
    delete hyCoM;
    delete text;
    delete f;
    
    return 0;
}



















