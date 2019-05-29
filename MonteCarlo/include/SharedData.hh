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
#include "XMLSettingsReader.hh"

#include <iostream>
#include <string>

class Survey {

 public :
  
	/** Type of detector - ZDC or RPD **/
    std::string detector;
    /** x_pos for this survey */
    double x_pos;
	/** y_pos for this survey */
    double y_pos;
	/** z_pos for this survey */
    double z_pos;
	/** cos_x for this survey */
    double cos_x;
	/** cos_y for this survey */
    double cos_y;
	/** cos_z for this survey */
    double cos_z;
	
};

class Alignment {

 public:

    /** X position of the Desy Table **/
    double x_table;
    /** Y position of the Desy Table **/
    double y_table;
    /** First detector met by the beam **/
    std::string upstream_Det;
    /** Second detector met by the beam **/
    std::string mid_Det;
    /** Third detector met by the beam **/
    std::string downstream_Det;
    /** GOLIATH magnet status **/
    bool magnet_On;
    /** Target in **/
    bool target_In;
    /** Lead absorber in **/
    bool lead_In;

};


class SharedData{
  
public:
  SharedData();
  SharedData( const std::string&, const std::string& );
  ~SharedData();

  // We do not want any copies of this class
  SharedData            ( const SharedData& ) = delete ;
  SharedData& operator= ( const SharedData& ) = delete ;

  // Only want MyRunManager to access these 

  void   Initialize();
  void   EndOfEvent();
  void   Finalize(); 
   
  template<class T> 
  void   AddOutputToZDCTree    ( const std::string&, T*);
  template<class T> 
  void   AddOutputToRPDTree    ( const std::string&, T*);
  template<class T> 
  void   AddOutputToFiberTree    ( const std::string&, T*);
  
  void   AddOutputHistogram ( TH1* );
   
  TEnv*  GetConfig        () { return m_config; }

  bool   DoPrint          ();
  inline TTree* GetZDCTree          () {return m_treeZDC;}
  inline TTree* GetRPDTree          () {return m_treeRPD;}
  inline TTree* GetFiberTree          () {return m_treeFiber;}
  
  inline int GetEventNo       () {return m_eventCounter;}
  
  Survey* GetSurvey(std::string name);
  Alignment* GetAlignment();
  
  void LoadConfigurationFile(int runNum, std::string _inFile = std::getenv("JZCaPA") + std::string("/Utils/Survey_2018.xml"));
  void LoadAlignmentFile(int runNum, std::string _inFile = std::getenv("JZCaPA") + std::string("/Utils/Alignment_2018.xml"));
   
   
   
private:
  TEnv*   m_config;
  TFile*  m_fout;
  TTree*  m_treeZDC;
  TTree*  m_treeRPD;
  TTree*  m_treeFiber;
  
  std::vector< TH1* > m_v_hists;
  
  int     m_eventCounter;

  std::string  m_outputFileName;
  std::string  m_configFileName;
  
   //XML parser
  XMLSettingsReader *m_XMLparser;

  //Alignment information for the given run
	Alignment* 	m_alignment;
  //Survey information for the given run
	Survey* m_survey;
	std::vector < Survey* > surveyEntries;
	
};

/** @brief Function to add an object to the tree
 *
 *  @param1 Name of branch
 *  @param2 Pointer to object
 *
 *  @return void
 */

template<class T>

void SharedData :: AddOutputToZDCTree( const std::string& name, T* pObj ){
  std::cout << "Adding " << name << std::endl;
   m_treeZDC->Branch( name.c_str(), pObj );
}

template<class T>

void SharedData :: AddOutputToRPDTree( const std::string& name, T* pObj ){
  std::cout << "Adding " << name << std::endl;
   m_treeRPD->Branch( name.c_str(), pObj );
}

template<class T>

void SharedData :: AddOutputToFiberTree( const std::string& name, T* pObj ){
  std::cout << "Adding " << name << std::endl;
   m_treeFiber->Branch( name.c_str(), pObj );
}

#endif
