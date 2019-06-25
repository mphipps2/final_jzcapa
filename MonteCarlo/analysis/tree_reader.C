#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include <TTree.h>
#include <TROOT.h>
#include "TTreeReaderValue.h"
#include "TMath.h"
#include "TBox.h"
#include <vector>
#include <deque>


using namespace std;
	
	int Col(int num){//takes in tile num and returns column
		if		(num == 0 || num == 4 || num == 8 	|| num == 12) return 1;
		else if	(num == 1 || num == 5 || num == 9 	|| num == 13) return 2;
		else if	(num == 2 || num == 6 || num == 10 	|| num == 14) return 3;
		else if	(num == 3 || num == 7 || num == 11 	|| num == 15) return 4;
		else return 0;
	}
	int Row(int num){//takes in tile num and returns row
		if		(num == 0 	|| num == 1 	|| num == 2 	|| num == 3) return 1;
		else if	(num == 4	|| num == 5 	|| num == 6 	|| num == 7) return 2;
		else if	(num == 8 	|| num == 9 	|| num == 10 	|| num == 11) return 3;
		else if	(num == 12	|| num == 13 	|| num == 14 	|| num == 15) return 4;
		else return 0;
	}

 void tree_reader() {
	
	auto myFile = TFile::Open("outputname.root");
	
	TTree *RPDtree= (TTree*)myFile->Get("RPDtree");
	TTree *Fibertree= (TTree*)myFile->Get("Fibertree");
	
	//Create rpd branch pointers
	TBranch *bX = 0, 	*bY = 0, 	*bZ = 0, 	*bpX = 0, 	*bpY = 0,	*bpZ = 0, *bEnergy=0, *bVx=0, *bVy=0, *bVz=0;
	TBranch *bPID = 0, 	*bID = 0, 	*bEventNo = 0 , *bNCherenkovs = 0, *bModNb=0, *bRadNb=0, *bRodNb=0;
	
	//Create fiber branch pointers
	TBranch *b1X = 0, 	*b1Y = 0, 	*b1Z = 0, 	*b1pX = 0, 	*b1pY = 0,	*b1pZ = 0, *b1Energy=0, *b1Vx=0, *b1Vy=0, *b1Vz=0;
	TBranch *b1PID = 0, 	*b1ID = 0, 	*b1EventNo = 0 , *b1NCherenkovs = 0, *b1ModNb=0, *b1RadNb=0, *b1RodNb=0;
	
	//Create rpd vectors
	std::vector<double> *X=0, *Y=0, *Z=0, *Px=0, *Py=0, *Pz=0, *Energy=0, *NCherenkovs=0;
	std::vector<int> *PID=0, *ID=0, *EventNo=0, *ModNb=0, *RadNb=0, *RodNb=0;
	
	//Create fiber vectors
	std::vector<double> *Xf=0, *Yf=0, *Zf=0, *Pxf=0, *Pyf=0, *Pzf=0, *Energyf=0, *NCherenkovsf=0;
	std::vector<int> *PIDf=0, *IDf=0, *EventNof=0, *ModNbf=0, *RadNbf=0, *RodNbf=0;

	//Set RPD branch addresses
	RPDtree->SetBranchAddress("X",&X,&bX);//
	RPDtree->SetBranchAddress("Y",&Y,&bY);//
	RPDtree->SetBranchAddress("Z",&Z,&bZ);//
	RPDtree->SetBranchAddress("Px",&Px,&bpX);//
	RPDtree->SetBranchAddress("Py",&Py,&bpY);//
	RPDtree->SetBranchAddress("Pz",&Pz,&bpZ);//
	RPDtree->SetBranchAddress("Pid",&PID,&bPID);//
	RPDtree->SetBranchAddress("Energy",&Energy,&bEnergy);//
	RPDtree->SetBranchAddress("ID",&ID,&bID);//
	RPDtree->SetBranchAddress("EventNo",&EventNo,&bEventNo);//
	RPDtree->SetBranchAddress("ModNb",&ModNb,&bModNb);//
	RPDtree->SetBranchAddress("RadNb",&RadNb,&bRadNb);//
	RPDtree->SetBranchAddress("RodNb",&RodNb,&bRodNb);//
	RPDtree->SetBranchAddress("NCherenkovs",&NCherenkovs,&bNCherenkovs);//

	//Set Fiber branch addresses
	Fibertree->SetBranchAddress("X",&Xf,&b1X);//
	Fibertree->SetBranchAddress("Y",&Yf,&b1Y);//
	Fibertree->SetBranchAddress("Z",&Zf,&b1Z);//
	Fibertree->SetBranchAddress("Px",&Pxf,&b1pX);//
	Fibertree->SetBranchAddress("Py",&Pyf,&b1pY);//
	Fibertree->SetBranchAddress("Pz",&Pzf,&b1pZ);//
	Fibertree->SetBranchAddress("Pid",&PIDf,&b1PID);//
	Fibertree->SetBranchAddress("Energy",&Energyf,&b1Energy);//
	Fibertree->SetBranchAddress("ID",&IDf,&b1ID);//
	Fibertree->SetBranchAddress("EventNo",&EventNof,&b1EventNo);//
	Fibertree->SetBranchAddress("ModNb",&ModNbf,&b1ModNb);//
	Fibertree->SetBranchAddress("RadNb",&RadNbf,&b1RadNb);//
	Fibertree->SetBranchAddress("RodNb",&RodNbf,&b1RodNb);//
	Fibertree->SetBranchAddress("NCherenkovs",&NCherenkovsf,&b1NCherenkovs);//
	
	//START ANALYSIS VARIABLES
	Int_t rpdCol[4], rpdRow[4];
	Int_t fiberCol[4], fiberRow[4];
	Int_t fiber_contribution[16];
	
	Double_t rpdColWeighted[4], rpdRowWeighted[4];
	Double_t CoMcol, CoMrow;
	Double_t value_holder=0;
	
	int totalCol = 0, totalRow = 0;
	double rpdColSum = 0, rpdRowSum =0;
	int tile;
	
	// x_start, y_start correspond to tile0s center which can be retrieved by running your geant in interactive mode and reading the position. the current setup is for the RPD centerred on (0,0)
	float x_start, y_start, x_offset, y_offset;
	 x_start  = 32.37;
	 y_start  = 30.75;
	 x_offset = 21.58;
	 y_offset = 20.5;

	 float tile_midx[4];
	 float tile_midy[4];

		//zeroing loop and tile midpoint calculation
	for(int i=0;i<4;i++){
		rpdCol[i] = 0;
		rpdRow[i] = 0;
		fiberCol[i] = 0;
		fiberRow[i] = 0;
		rpdColWeighted[i]= 0;
		rpdRowWeighted[i]= 0;
		
		tile_midx[i]=x_start-(i*x_offset);
	 	tile_midy[i]=y_start-(i*y_offset);
	}
	for(int i=0;i<16;i++){
		fiber_contribution[i]=0;
	}
	//END ANALYSIS VARIABLES
	
	
	//START CANVAS CREATION
	TCanvas *canv1 = new TCanvas("c1" , "TEST CANVAS" , 1800, 1000);
	canv1->Divide(2,1);
	//END CANVAS CREATION
	
	
	//START HISTOGRAM CREATION
	const int rpd_highx	=  50;
	const int rpd_lowx	= -50;
	const int rpd_highy	=  50;
	const int rpd_lowy	= -50;
		
	TH2F *rpdCoM = new TH2F("rpdCoM","rpdCoM",200,rpd_lowx,rpd_highx,200,rpd_lowy,rpd_highy);
	
	TH2F *fiberContribtuion = new TH2F("fiber_tile","fiber_tile",4,0,4,4,0,4);
	//END HISTOGRAM CREATION
	

	//begin loop over events
	for (int q=0;q<RPDtree->GetEntries() ; q++)
	{
		Long64_t tentry = RPDtree->LoadTree(q);
		
		//get info for rpd event q
		bX->GetEntry(tentry);
		bY->GetEntry(tentry);
		bZ->GetEntry(tentry);
		bpX->GetEntry(tentry);
		bpY->GetEntry(tentry);
		bpZ->GetEntry(tentry);
		bPID->GetEntry(tentry);
		bEnergy->GetEntry(tentry); //Energy in units of MeV
		bID->GetEntry(tentry);
		bEventNo->GetEntry(tentry);
		bModNb->GetEntry(tentry);
		bRadNb->GetEntry(tentry);
		bRodNb->GetEntry(tentry); //FOR THE RPD RodNb CORRESPONDS TO THE TILE# 
		bNCherenkovs->GetEntry(tentry);

		tentry = Fibertree->LoadTree(q);

		//get info for fiber event q
		b1X->GetEntry(tentry);
		b1Y->GetEntry(tentry);
		b1Z->GetEntry(tentry);
		b1pX->GetEntry(tentry);
		b1pY->GetEntry(tentry);
		b1pZ->GetEntry(tentry);
		b1PID->GetEntry(tentry);
		b1Energy->GetEntry(tentry); //Energy in units of MeV
		b1ID->GetEntry(tentry);
		b1EventNo->GetEntry(tentry);
		b1ModNb->GetEntry(tentry);
		b1RadNb->GetEntry(tentry);
		b1RodNb->GetEntry(tentry);
		b1NCherenkovs->GetEntry(tentry);

		
		for (int k=0; k<X->size(); k++){	//begin loop over rpd hits 
			if(NCherenkovs->at(k) > 0){		//ignore empty vectors
				
				rpdCol[Col(RodNb->at(k)) - 1] += NCherenkovs->at(k); //FOR THE RPD RodNb CORRESPONDS TO THE TILE# 
				rpdRow[Row(RodNb->at(k)) - 1] += NCherenkovs->at(k); //FOR THE RPD RodNb CORRESPONDS TO THE TILE# 

			}
		}									//end loop over rpd hits 
		for (int k=0; k<Xf->size(); k++){	//begin loop over fiber hits
			if(NCherenkovsf->at(k) > 0){	//ignore empty vectors
					
				tile = 4*(RodNbf->at(k)%4) + RodNbf->at(k)/16; //associates the readout fiber with the tile it reads out

				fiber_contribution[tile]+=NCherenkovsf->at(k); //count how many cherenkovs are added to a tile by its readout fibers
			}
		}									//end loop over fiber hits

				for(int i=0;i<4;i++){					//begin loop over rpd col/row
				
					//CALCULATION OF RPD CoM
					totalCol += rpdCol[i];
					totalRow += rpdRow[i];

					value_holder = tile_midx[i] * (double) rpdCol[i];
					rpdColWeighted[i] += value_holder;
					
					value_holder = tile_midy[i] * (double) rpdRow[i];
					rpdRowWeighted[i] += value_holder;

					rpdColSum += rpdColWeighted[i];
					rpdRowSum += rpdRowWeighted[i];

					//Zero Arrays
					rpdCol[i] = 0;
					rpdRow[i] = 0;
					rpdColWeighted[i]= 0;
					rpdRowWeighted[i]= 0;
					}

				//final calculations
				CoMcol = rpdColSum / totalCol;
				CoMrow = rpdRowSum / totalRow;
				
				
				//fill CoM histograms
				rpdCoM->Fill(CoMcol,CoMrow);
				

				//Zero Variables
				totalCol = 0;
				totalRow = 0;
				rpdColSum = 0;
				rpdRowSum = 0;
	}
	//end loop over events

	
	int counter=0;
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			fiberContribtuion->SetBinContent(4-j,4-i, fiber_contribution[counter]);
			counter++;
		}
	}
	
	
	
	//DRAW HISTOGRAMS
	canv1->cd(1);
	rpdCoM->Draw("colz");
	canv1->cd(2);
	fiberContribtuion->Draw("colz text");
	
		
	RPDtree->ResetBranchAddresses();
	Fibertree->ResetBranchAddresses();
}