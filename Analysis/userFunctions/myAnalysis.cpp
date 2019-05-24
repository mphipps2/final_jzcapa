/** @file AnalysisExample.cpp
 *  @brief Example of a simple analysis: read and plot all the w
 *
 *
 *  @author Yakov Kulinich, Riccardo Longo
 *  @bug No known bugs.
 */

#include <stdlib.h>
#include "DataReader.h"
#include "WFAnalysis.h"
#include "ZDCAnalysis.h"
#include "RPDAnalysis.h"
#include "EventTimer.h"
#include "TSystem.h"

using namespace std;


int main(int argc, char *argv[]){

  // can put arguments here, for now I will just do defaults.

  int nCh    = 20;   // 5 DRS4 x 4 ch/board - 16 RPD channels
  int nSamp  = 1024; // Default number of samples?
  int runNum = atoi(argv[1]);

  //Make the file name from the run number
  int start[] = {79,  152, 190, 202, 215, 258, 280, 312, 335, 351, 362, 382, 412 };
  int stop[]  = {112, 171, 200, 213, 231, 277, 296, 322, 350, 359, 381, 391, 413 };
  int scanNum = 0;
  for(int i = 0; i < 13; i++){
      if(start[i]<=runNum && stop[i]>=runNum){
          scanNum = i+1;
      }
  }
  if(scanNum == 0){
      cout << "File does not exist... exiting" << endl;
      return 0;
  }
  string fNameIn = Form("/data/phenix/data/TestBeam2018/new_processed_runs/Scan%d/ZDCBeamTestRun%d.root", scanNum, runNum); // !! Change for your test !!
  gSystem->Exec( Form("mkdir -p /data/phenix/data/TestBeam2018/Post_processing/run%d", runNum) );

  // DataReader is the main class. It reads data and also
  // has analysis classes in it. User should only have to
  // modify the analysis classes and add output in them.
  // User has to add their analysis to DataReader.
  DataReader* r = new DataReader( nCh, nSamp, fNameIn, runNum );
  r->SetVerbosity(1);

  r->SelectDetectorForAnalysis(true,true,true);
  r->AddPreAnalysis( new WFAnalysis() );
  r->AddDetectorAnalysis( new ZDCAnalysis() );
  r->AddDetectorAnalysis( new RPDAnalysis() );
  r->LoadConfigurationFile();
  r->LoadAlignmentFile();
  r->LoadTimingFile();
  r->SetOutputDirectory( Form("/data/phenix/data/TestBeam2018/Post_processing/run%d/", runNum) );
  r->EnablePlotLabel();
  
  EventTimer timer(1000, r, kFALSE);
  timer.TurnOn();
  
  r->Run();
  
  timer.TurnOff();
  std::cout << std::endl << "Finished!" << std::endl;
  
  delete r;
  
  return 0;
}
