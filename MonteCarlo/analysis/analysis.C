#include <iostream>




void analysis() {


      char fileName[256];
      sprintf(fileName,"/home/mike/ZDC/100Neutrons75GeV.root");
      TFile *file = TFile::Open(fileName);
      if (!file) {cerr << "file not open " << endl; return 1;}
      TTree *hitTree = new TTree("Hit","Hit");   
      hitTree = (TTree*) file->Get("tree");
  
      Int_t treeSize = (Int_t) hitTree->GetEntries();

      /*      
      char h1Name[256];
      sprintf(h1Name,"maxTB%d",z);
      TH1F *h1 = new TH1F(h1Name,h1Name,27, -0.5, 26.5);
      maxTBHisto.push_back(h1); 
      char h2Name[256];
      sprintf(h2Name,"meanTime%d",z);
      TH1F *h2 = new TH1F(h2Name,h2Name,115,10 , 18);
      meanTimeHisto.push_back(h2);
      */
      int ModNb = 0; int RadNb = 0; int Pid = 0; int EventNb = 0;
      double EDep = 0.; double X = 0.; double Y = 0.; double Z = 0.; double Charge = 0.;
      

      
      TCanvas *c1 = new TCanvas("c1","c1");
      hitTree->Draw("EDep:X");
      c1->SaveAs("EDepVsX.root");
      c1->SaveAs("EDepVsX.png");
      hitTree->Draw("EDep:Y");
      c1->SaveAs("EDepVsY.root");
      c1->SaveAs("EDepVsY.png");
      hitTree->Draw("EDep:Z");
      c1->SaveAs("EDepVsZ.root");
      c1->SaveAs("EDepVsZ.png");

      hitTree->Draw("X","Pid == 0");
      c1->SaveAs("CherenkovVsX.root");
      c1->SaveAs("CherenkovVsX.png");
      hitTree->Draw("Y","Pid == 0");
      c1->SaveAs("CherenkovVsY.root");
      c1->SaveAs("CherenkovVsY.png");
      hitTree->Draw("Z","Pid == 0");
      c1->SaveAs("CherenkovVsZ.root");
      c1->SaveAs("CherenkovVsZ.png");

      hitTree->Draw("X","Pid == 22");
      c1->SaveAs("GammaVsX.root");
      c1->SaveAs("GammaVsX.png");

      hitTree->Draw("X","Pid == 111");
      c1->SaveAs("Pi0VsX.root");
      c1->SaveAs("Pi0VsX.png");

      hitTree->Draw("X","Pid == 211");
      c1->SaveAs("Pi+VsX.root");
      c1->SaveAs("Pi+VsX.png");

      hitTree->Draw("X");
      c1->SaveAs("TotalX.root");
      c1->SaveAs("TotalX.png");
      hitTree->Draw("Y");
      c1->SaveAs("TotalY.root");
      c1->SaveAs("TotalY.png");
      hitTree->Draw("Z");
      c1->SaveAs("TotalZ.root");
      c1->SaveAs("TotalZ.png");

      delete c1;

      hitTree->SetBranchAddress("ModNb",&ModNb);
      hitTree->SetBranchAddress("RadNb",&RadNb);
      hitTree->SetBranchAddress("EDep",&EDep);
      hitTree->SetBranchAddress("Pid",&Pid);     
      hitTree->SetBranchAddress("X",&X);
      hitTree->SetBranchAddress("Y",&Y);
      hitTree->SetBranchAddress("Z",&Z);
      hitTree->SetBranchAddress("EventNo",&EventNb);
      hitTree->SetBranchAddress("Charge",&Charge);
      TH1F *photonHist = new TH1F("cherenkovYield","cherenkovYield",10,1000,100000);
      TF1 *f1 = new TF1("f1","[0]+[1]*x+[2]*x*x",-45,45);
      f1->SetParameter(0,0.162133);
      f1->SetParameter(1,-4.27359e-05);
      f1->SetParameter(2,6.55089e-05);
      float prob = 0;
      float aggregateprob = 0;        
      int prevEventNb = 0;
      for (int i=0; i < treeSize; ++i) {
        hitTree->GetEntry(i);                     
        while (eventNb == prevEventNb) {
          prob = f1->Eval(X);
          aggregateprob += prob;                 
	}
        float totalProb = aggregateProb/nOptPhotons;
        float eventYield = totalProb * nOptPhotons;
	photonHist->Fill(eventYield);
	std::cout << " yield " << eventYield << std::endl;
	prevEventNb = eventNb;
      } 
      
      
      delete hitTree;

      
  }











}
