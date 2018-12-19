/** @file SharedData.h
 *  @brief Function prototypes for SharedData.
 *
 *  This contains the prototypes and members 
 *  for SharedData.
 *
 *  @author Yakov Kulinich
 *  @bug No known bugs.
 */

#ifndef YKANALYSIS_SHAREDDATA_H
#define YKANALYSIS_SHAREDDATA_H

#include <TEnv.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <iostream>
#include <string>

class SharedData{
  
public:
  SharedData();
  SharedData( const std::string&, const std::string& );
  ~SharedData();

  // We do not want any copies of this class
  SharedData            ( const SharedData& ) = delete ;
  SharedData& operator= ( const SharedData& ) = delete ;

  // Only want MyRunManager to access these 
public:
  void   Initialize();
  void   EndOfEvent();
  void   Finalize(); 
   
public:
  template<class T> 
  void   AddOutputToTree    ( const std::string&, T*);

  void   AddOutputHistogram ( TH1* );
   
  TEnv*  GetConfig        () { return m_config; }

  bool   DoPrint          ();
  inline TTree* GetTree          () {return m_tree;}
  inline int GetEventNo       () {return m_eventCounter;}
private:
  TEnv*   m_config;
  TFile*  m_fout;
  TTree*  m_tree;
  std::vector< TH1* > m_v_hists;
  
  int     m_eventCounter;

  std::string  m_outputFileName;
  std::string  m_configFileName;
};

/** @brief Function to add an object to the tree
 *
 *  @param1 Name of branch
 *  @param2 Pointer to object
 *
 *  @return void
 */

template<class T>

void SharedData :: AddOutputToTree( const std::string& name, T* pObj ){
  std::cout << "Adding " << name << std::endl;
  m_tree->Branch( name.c_str(), pObj );
}



#endif
