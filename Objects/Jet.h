#ifndef Jet_H
#define Jet_H

#include "Particle.h"


class Jet : public Particle{

 public:

  Jet(){
     m_nTracks=0;
     m_jetArea=0;
     m_pileup=0;
     m_numberOfDaughters=0; 
     m_neutralEmEnergyFraction=0;
     m_neutralHadronEnergyFraction=0;
     m_chargedEmEnergyFraction=0;
     m_chargedHadronEnergyFraction=0;
     m_muonEnergyFraction=0;
     m_photonEnergyFraction=0;
     m_chargedMultiplicity=0;
     m_neutralMultiplicity=0;
     m_muonMultiplicity=0; 
     m_electronMultiplicity=0;
     m_photonMultiplicity=0;
     m_btag_simpleSecondaryVertexHighEff=0;
     m_btag_simpleSecondaryVertexHighPur=0;
     m_btag_combinedSecondaryVertex=0;
     m_btag_combinedSecondaryVertexMVA=0;
     m_btag_jetBProbability=0;
     m_btag_jetProbability=0;
     m_JEC_uncertainty=0;
     m_JEC_factor_raw=0;
     m_genjet_pt=0;
     m_genjet_eta=0;
     m_genjet_phi=0;
     m_genjet_energy=0;
     m_genparticles_indices.clear();
  };

  ~Jet(){
  };

  LorentzVector genjet_v4() const{
    LorentzVector v4;
    v4.SetPt(m_genjet_pt);
    v4.SetEta(m_genjet_eta);
    v4.SetPhi(m_genjet_phi);
    v4.SetE(m_genjet_energy);
    return v4;
  };

  int nTracks() const{return m_nTracks;}
  float jetArea() const{return m_jetArea;}
  float pileup() const{return m_pileup;}
  int numberOfDaughters() const{return m_numberOfDaughters;} 
  float neutralEmEnergyFraction() const{return m_neutralEmEnergyFraction;}
  float neutralHadronEnergyFraction() const{return m_neutralHadronEnergyFraction;}
  float chargedEmEnergyFraction() const{return m_chargedEmEnergyFraction;}
  float chargedHadronEnergyFraction() const{return m_chargedHadronEnergyFraction;}
  float muonEnergyFraction() const{return m_muonEnergyFraction;}
  float photonEnergyFraction() const{return m_photonEnergyFraction;}
  int chargedMultiplicity() const{return m_chargedMultiplicity;}
  int neutralMultiplicity() const{return m_neutralMultiplicity;}
  int muonMultiplicity() const{return m_muonMultiplicity;} 
  int electronMultiplicity() const{return m_electronMultiplicity;}
  int photonMultiplicity() const{return m_photonMultiplicity;}
  float btag_simpleSecondaryVertexHighEff() const{return m_btag_simpleSecondaryVertexHighEff;}
  float btag_simpleSecondaryVertexHighPur() const{return m_btag_simpleSecondaryVertexHighPur;}
  float btag_combinedSecondaryVertex() const{return m_btag_combinedSecondaryVertex;}
  float btag_combinedSecondaryVertexMVA() const{return m_btag_combinedSecondaryVertexMVA;}
  float btag_jetBProbability() const{return m_btag_jetBProbability;}
  float btag_jetProbability() const{return m_btag_jetProbability;}
  float JEC_uncertainty() const{return m_JEC_uncertainty;}
  float JEC_factor_raw() const{return m_JEC_factor_raw;}
  float genjet_pt() const{return m_genjet_pt;}
  float genjet_eta() const{return m_genjet_eta;}
  float genjet_phi() const{return m_genjet_phi;}
  float genjet_energy() const{return m_genjet_energy;}
  std::vector<unsigned int> genparticles_indices() const{return m_genparticles_indices;}

  void set_nTracks(int x){m_nTracks=x;}
  void set_jetArea(float x){m_jetArea=x;}
  void set_pileup(float x){m_pileup=x;}
  void set_numberOfDaughters(int x){m_numberOfDaughters=x;} 
  void set_neutralEmEnergyFraction(float x){m_neutralEmEnergyFraction=x;}
  void set_neutralHadronEnergyFraction(float x){m_neutralHadronEnergyFraction=x;}
  void set_chargedEmEnergyFraction(float x){m_chargedEmEnergyFraction=x;}
  void set_chargedHadronEnergyFraction(float x){m_chargedHadronEnergyFraction=x;}
  void set_muonEnergyFraction(float x){m_muonEnergyFraction=x;}
  void set_photonEnergyFraction(float x){m_photonEnergyFraction=x;}
  void set_chargedMultiplicity(int x){m_chargedMultiplicity=x;}
  void set_neutralMultiplicity(int x){m_neutralMultiplicity=x;}
  void set_muonMultiplicity(int x){m_muonMultiplicity=x;} 
  void set_electronMultiplicity(int x){m_electronMultiplicity=x;}
  void set_photonMultiplicity(int x){m_photonMultiplicity=x;}
  void set_btag_simpleSecondaryVertexHighEff(float x){m_btag_simpleSecondaryVertexHighEff=x;}
  void set_btag_simpleSecondaryVertexHighPur(float x){m_btag_simpleSecondaryVertexHighPur=x;}
  void set_btag_combinedSecondaryVertex(float x){m_btag_combinedSecondaryVertex=x;}
  void set_btag_combinedSecondaryVertexMVA(float x){m_btag_combinedSecondaryVertexMVA=x;}
  void set_btag_jetBProbability(float x){m_btag_jetBProbability=x;}
  void set_btag_jetProbability(float x){m_btag_jetProbability=x;}
  void set_JEC_uncertainty(float x){m_JEC_uncertainty=x;}
  void set_JEC_factor_raw(float x){m_JEC_factor_raw=x;}
  void set_genjet_pt(float x){m_genjet_pt=x;}
  void set_genjet_eta(float x){m_genjet_eta=x;}
  void set_genjet_phi(float x){m_genjet_phi=x;}
  void set_genjet_energy(float x){m_genjet_energy=x;}
  void add_genparticles_index(unsigned int x){m_genparticles_indices.push_back(x);}

 private:
  
  int m_nTracks;
  float m_jetArea;
  float m_pileup;
  int m_numberOfDaughters; 
  float m_neutralEmEnergyFraction;
  float m_neutralHadronEnergyFraction;
  float m_chargedEmEnergyFraction;
  float m_chargedHadronEnergyFraction;
  float m_muonEnergyFraction;
  float m_photonEnergyFraction;
  int m_chargedMultiplicity;
  int m_neutralMultiplicity;
  int m_muonMultiplicity; 
  int m_electronMultiplicity;
  int m_photonMultiplicity;
  float m_btag_simpleSecondaryVertexHighEff;
  float m_btag_simpleSecondaryVertexHighPur;
  float m_btag_combinedSecondaryVertex;
  float m_btag_combinedSecondaryVertexMVA;
  float m_btag_jetBProbability;
  float m_btag_jetProbability;
  float m_JEC_uncertainty;
  float m_JEC_factor_raw;
  float m_genjet_pt;
  float m_genjet_eta;
  float m_genjet_phi;
  float m_genjet_energy;

  std::vector<unsigned int> m_genparticles_indices;

};

#endif
