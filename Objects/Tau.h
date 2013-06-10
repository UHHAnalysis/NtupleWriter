#ifndef Tau_H
#define Tau_H

#include "Particle.h"

/**
 *  @short tau class
 *  @author Thomas Peiffer
 */

class Tau : public Particle{

 public:
  Tau(){
/*     m_leadPFCand_px=0; */
/*     m_leadPFCand_py=0; */
/*     m_leadPFCand_pz=0; */
    m_decayModeFinding=false; 
    //m_byVLooseCombinedIsolationDeltaBetaCorr =false;
    m_byLooseCombinedIsolationDeltaBetaCorr=false; 
    m_byMediumCombinedIsolationDeltaBetaCorr=false;
    m_byTightCombinedIsolationDeltaBetaCorr=false; 
    m_byLooseIsolationMVA=false;
    m_byMediumIsolationMVA=false;   
    m_byTightIsolationMVA=false;
    m_byLooseIsolationMVA2=false;
    m_byMediumIsolationMVA2=false;   
    m_byTightIsolationMVA2=false;
    m_byLooseCombinedIsolationDeltaBetaCorr3Hits=false;
    m_byMediumCombinedIsolationDeltaBetaCorr3Hits=false;
    m_byTightCombinedIsolationDeltaBetaCorr3Hits=false;
    m_againstElectronLooseMVA3 =false;
    m_againstElectronMediumMVA3=false;
    m_againstElectronTightMVA3=false ;
    m_againstElectronVTightMVA3=false ;
    m_againstMuonLoose2=false ;
    m_againstMuonMedium2=false;
    m_againstMuonTight2=false;
    m_byIsolationMVAraw=0;
    m_byIsolationMVA2raw=0;
    m_decayMode=-1;
    m_byCombinedIsolationDeltaBetaCorrRaw=-1;
    m_byCombinedIsolationDeltaBetaCorrRaw3Hits=-1;  

  };

  ~Tau(){
  };

/*   float leadPFCand_px() const{return m_leadPFCand_px;} */
/*   float leadPFCand_py() const{return m_leadPFCand_py;} */
/*   float leadPFCand_pz() const{return m_leadPFCand_pz;} */

  bool decayModeFinding() const{return m_decayModeFinding;} 
  //bool byVLooseCombinedIsolationDeltaBetaCorr () const{return m_byVLooseCombinedIsolationDeltaBetaCorr;}
  bool byLooseCombinedIsolationDeltaBetaCorr() const{return m_byLooseCombinedIsolationDeltaBetaCorr;} 
  bool byMediumCombinedIsolationDeltaBetaCorr() const{return m_byMediumCombinedIsolationDeltaBetaCorr;}
  bool byTightCombinedIsolationDeltaBetaCorr() const{return m_byTightCombinedIsolationDeltaBetaCorr;} 
  bool byLooseIsolationMVA() const{return m_byLooseIsolationMVA;}
  bool byMediumIsolationMVA() const{return m_byMediumIsolationMVA;}   
  bool byTightIsolationMVA() const{return m_byTightIsolationMVA ;}
  bool byLooseIsolationMVA2() const{return m_byLooseIsolationMVA2;}
  bool byMediumIsolationMVA2() const{return m_byMediumIsolationMVA2;}   
  bool byTightIsolationMVA2() const{return m_byTightIsolationMVA2;}
  bool byLooseCombinedIsolationDeltaBetaCorr3Hits() const{return m_byLooseCombinedIsolationDeltaBetaCorr3Hits;}
  bool byMediumCombinedIsolationDeltaBetaCorr3Hits() const{return m_byMediumCombinedIsolationDeltaBetaCorr3Hits;}
  bool byTightCombinedIsolationDeltaBetaCorr3Hits() const{return m_byTightCombinedIsolationDeltaBetaCorr3Hits;}
  bool againstElectronLooseMVA3 () const{return m_againstElectronLooseMVA3;}
  bool againstElectronMediumMVA3() const{return m_againstElectronMediumMVA3;}
  bool againstElectronTightMVA3 () const{return m_againstElectronTightMVA3;}
  bool againstElectronVTightMVA3 () const{return m_againstElectronVTightMVA3;}
  bool againstMuonLoose2() const{return m_againstMuonLoose2;}
  bool againstMuonMedium2() const{return m_againstMuonMedium2;}
  bool againstMuonTight2() const{return m_againstMuonTight2;}
  float byIsolationMVAraw() const{return m_byIsolationMVAraw;}
  float byIsolationMVA2raw() const{return m_byIsolationMVA2raw;}
  float byCombinedIsolationDeltaBetaCorrRaw() const{return m_byCombinedIsolationDeltaBetaCorrRaw;}
  float byCombinedIsolationDeltaBetaCorrRaw3Hits() const{return m_byCombinedIsolationDeltaBetaCorrRaw3Hits;}
    
  int decayMode() const{return m_decayMode;}

/*   void set_leadPFCand_px(float x){m_leadPFCand_px=x;} */
/*   void set_leadPFCand_py(float x){m_leadPFCand_py=x;} */
/*   void set_leadPFCand_pz(float x){m_leadPFCand_pz=x;} */

  void set_decayModeFinding(bool x){m_decayModeFinding=x;} 
  //void set_byVLooseCombinedIsolationDeltaBetaCorr (bool x){m_byVLooseCombinedIsolationDeltaBetaCorr=x;}
  void set_byLooseCombinedIsolationDeltaBetaCorr(bool x){m_byLooseCombinedIsolationDeltaBetaCorr=x;} 
  void set_byMediumCombinedIsolationDeltaBetaCorr(bool x){m_byMediumCombinedIsolationDeltaBetaCorr=x;}
  void set_byTightCombinedIsolationDeltaBetaCorr(bool x){m_byTightCombinedIsolationDeltaBetaCorr=x;} 
  void set_byLooseIsolationMVA(bool x) {m_byLooseIsolationMVA=x;}
  void set_byMediumIsolationMVA(bool x) {m_byMediumIsolationMVA=x;}   
  void set_byTightIsolationMVA(bool x) {m_byTightIsolationMVA=x ;}
  void set_byLooseIsolationMVA2(bool x) {m_byLooseIsolationMVA2=x;}
  void set_byMediumIsolationMVA2(bool x) {m_byMediumIsolationMVA2=x;}   
  void set_byTightIsolationMVA2(bool x) {m_byTightIsolationMVA2=x;}
  void set_byLooseCombinedIsolationDeltaBetaCorr3Hits(bool x) {m_byLooseCombinedIsolationDeltaBetaCorr3Hits=x;}
  void set_byMediumCombinedIsolationDeltaBetaCorr3Hits(bool x) {m_byMediumCombinedIsolationDeltaBetaCorr3Hits=x;}
  void set_byTightCombinedIsolationDeltaBetaCorr3Hits(bool x) {m_byTightCombinedIsolationDeltaBetaCorr3Hits=x;}
  void set_againstElectronLooseMVA3 (bool x){m_againstElectronLooseMVA3=x;}
  void set_againstElectronMediumMVA3(bool x){m_againstElectronMediumMVA3=x;}
  void set_againstElectronTightMVA3 (bool x){m_againstElectronTightMVA3=x;}
  void set_againstElectronVTightMVA3 (bool x){m_againstElectronVTightMVA3=x;}
  void set_againstMuonLoose2 (bool x){m_againstMuonLoose2=x;}
  void set_againstMuonMedium2(bool x){m_againstMuonMedium2=x;}
  void set_againstMuonTight2(bool x){m_againstMuonTight2=x;}
  void set_byIsolationMVAraw(float x){ m_byIsolationMVAraw=x;}
  void set_byIsolationMVA2raw(float x){ m_byIsolationMVA2raw=x;}  
  void set_decayMode(int x){m_decayMode=x;}
  void set_byCombinedIsolationDeltaBetaCorrRaw(float x) {m_byCombinedIsolationDeltaBetaCorrRaw=x;}
  void set_byCombinedIsolationDeltaBetaCorrRaw3Hits(float x) {m_byCombinedIsolationDeltaBetaCorrRaw3Hits=x;} 

 private:
/*   float m_leadPFCand_px; */
/*   float m_leadPFCand_py; */
/*   float m_leadPFCand_pz; */

  bool m_decayModeFinding; 
  //bool m_byVLooseCombinedIsolationDeltaBetaCorr ;
  bool m_byLooseCombinedIsolationDeltaBetaCorr; 
  bool m_byMediumCombinedIsolationDeltaBetaCorr;
  bool m_byTightCombinedIsolationDeltaBetaCorr; 
  bool m_byLooseIsolationMVA;
  bool m_byMediumIsolationMVA;  
  bool m_byTightIsolationMVA;
  bool m_byLooseIsolationMVA2;
  bool m_byMediumIsolationMVA2;  
  bool m_byTightIsolationMVA2;
  bool m_byLooseCombinedIsolationDeltaBetaCorr3Hits;
  bool m_byMediumCombinedIsolationDeltaBetaCorr3Hits;
  bool m_byTightCombinedIsolationDeltaBetaCorr3Hits;
  bool m_againstElectronLooseMVA3 ;
  bool m_againstElectronMediumMVA3;
  bool m_againstElectronTightMVA3 ;
  bool m_againstElectronVTightMVA3 ;
  bool m_againstMuonLoose2 ;
  bool m_againstMuonMedium2;
  bool m_againstMuonTight2;
  float m_byIsolationMVAraw;
  float m_byIsolationMVA2raw;
  float m_byCombinedIsolationDeltaBetaCorrRaw;
  float m_byCombinedIsolationDeltaBetaCorrRaw3Hits;  
  int m_decayMode;

};

#endif
