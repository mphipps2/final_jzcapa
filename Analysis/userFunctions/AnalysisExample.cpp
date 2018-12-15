/** @file AnalysisExample.cpp
 *  @brief Example of a simple analysis: read and plot all the w
 *
 *
 *  @author Yakov Kulinich, Riccardo Longo
 *  @bug No known bugs.
 */

#include "DataReader.h"

using namespace std;

int main(int argc, char *argv[]){

  // can put arguments here, for now I will just do defaults.

  int nCh    = 20;   // 5 DRS4 x 4 ch/board - 16 RPD channels
  int nSamp  = 1024; // Default number of samples?
  int runNum = 171;   // !! Change for your test !!

  string fNameIn = "TreeZDCBeamTestRun171.root"; // !! Change for your test !!
 
  DataReader* r = new DataReader( nCh, nSamp, fNameIn, runNum );
  r->Initialize();
  r->ProcessEvents();
  r->Finalize();
  
  return 0;
}

