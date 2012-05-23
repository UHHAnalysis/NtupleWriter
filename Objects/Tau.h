#ifndef Tau_H
#define Tau_H

#include "Particle.h"

class Tau : public Particle{

 public:
  Tau(){
    m_leadPFCand_px=0;
    m_leadPFCand_py=0;
    m_leadPFCand_pz=0;
    m_decayModeFinding=false; 
    m_byVLooseCombinedIsolationDeltaBetaCorr =false;
    m_byLooseCombinedIsolationDeltaBetaCorr=false; 
    m_byMediumCombinedIsolationDeltaBetaCorr=false;
    m_byTightCombinedIsolationDeltaBetaCorr=false; 
    m_againstElectronLoose =false;
    m_againstElectronMedium=false;
    m_againstElectronTight=false ;
    m_againstElectronMVA =false;
    m_againstMuonLoose=false ;
    m_againstMuonMedium=false;
    m_againstMuonTight=false;
  };

  ~Tau(){
  };

  float leadPFCand_px() const{return m_leadPFCand_px;}
  float leadPFCand_py() const{return m_leadPFCand_py;}
  float leadPFCand_pz() const{return m_leadPFCand_pz;}

  bool decayModeFinding() const{return m_decayModeFinding;} 
  bool byVLooseCombinedIsolationDeltaBetaCorr () const{return m_byVLooseCombinedIsolationDeltaBetaCorr;}
  bool byLooseCombinedIsolationDeltaBetaCorr() const{return m_byLooseCombinedIsolationDeltaBetaCorr;} 
  bool byMediumCombinedIsolationDeltaBetaCorr() const{return m_byMediumCombinedIsolationDeltaBetaCorr;}
  bool byTightCombinedIsolationDeltaBetaCorr() const{return m_byTightCombinedIsolationDeltaBetaCorr;} 
  bool againstElectronLoose () const{return m_againstElectronLoose;}
  bool againstElectronMedium() const{return m_againstElectronMedium;}
  bool againstElectronTight () const{return m_againstElectronTight;}
  bool againstElectronMVA () const{return m_againstElectronMVA;}
  bool againstMuonLoose () const{return m_againstMuonLoose;}
  bool againstMuonMedium() const{return m_againstMuonMedium;}
  bool againstMuonTight() const{return m_againstMuonTight;}

  void set_leadPFCand_px(float x){m_leadPFCand_px=x;}
  void set_leadPFCand_py(float x){m_leadPFCand_py=x;}
  void set_leadPFCand_pz(float x){m_leadPFCand_pz=x;}

  void set_decayModeFinding(bool x){m_decayModeFinding=x;} 
  void set_byVLooseCombinedIsolationDeltaBetaCorr (bool x){m_byVLooseCombinedIsolationDeltaBetaCorr=x;}
  void set_byLooseCombinedIsolationDeltaBetaCorr(bool x){m_byLooseCombinedIsolationDeltaBetaCorr=x;} 
  void set_byMediumCombinedIsolationDeltaBetaCorr(bool x){m_byMediumCombinedIsolationDeltaBetaCorr=x;}
  void set_byTightCombinedIsolationDeltaBetaCorr(bool x){m_byTightCombinedIsolationDeltaBetaCorr=x;} 
  void set_againstElectronLoose (bool x){m_againstElectronLoose=x;}
  void set_againstElectronMedium(bool x){m_againstElectronMedium=x;}
  void set_againstElectronTight (bool x){m_againstElectronTight=x;}
  void set_againstElectronMVA (bool x){m_againstElectronMVA=x;}
  void set_againstMuonLoose (bool x){m_againstMuonLoose=x;}
  void set_againstMuonMedium(bool x){m_againstMuonMedium=x;}
  void set_againstMuonTight(bool x){m_againstMuonTight=x;}

 private:
  float m_leadPFCand_px;
  float m_leadPFCand_py;
  float m_leadPFCand_pz;

  bool m_decayModeFinding; 
  bool m_byVLooseCombinedIsolationDeltaBetaCorr ;
  bool m_byLooseCombinedIsolationDeltaBetaCorr; 
  bool m_byMediumCombinedIsolationDeltaBetaCorr;
  bool m_byTightCombinedIsolationDeltaBetaCorr; 
  bool m_againstElectronLoose ;
  bool m_againstElectronMedium;
  bool m_againstElectronTight ;
  bool m_againstElectronMVA ;
  bool m_againstMuonLoose ;
  bool m_againstMuonMedium;
  bool m_againstMuonTight;

};

#endif
