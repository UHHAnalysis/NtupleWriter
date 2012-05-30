#ifndef Particle_H
#define Particle_H

#include <vector>
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiE4D.h"
#include "TObject.h"


#ifndef PI
#define PI       3.14159265358979323846264338328     
#endif

typedef ROOT::Math::LorentzVector< ROOT::Math::PtEtaPhiE4D< Double32_t > > LorentzVector;

class Particle{
 public:
  Particle(){
    m_charge=0;
    m_pt=0; 
    m_eta=0; 
    m_phi=0; 
    m_energy=0; 
  };

  ~Particle(){
  };

  LorentzVector v4() const{
    LorentzVector v4;
    v4.SetPt(m_pt);
    v4.SetEta(m_eta);
    v4.SetPhi(m_phi);
    v4.SetE(m_energy);
    return v4;
  };

  float charge() const{return m_charge;}
  float pt() const {return m_pt;}
  float eta() const{return m_eta;}
  float phi() const{return m_phi;}
  float energy() const{return m_energy;}

  void set_charge(float charge){m_charge=charge;}
  void set_pt(float pt){m_pt=pt;}  
  void set_eta(float eta){m_eta=eta;}
  void set_phi(float phi){m_phi=phi;}
  void set_energy(float energy){m_energy=energy;}

  void set_v4(LorentzVector v4){
    set_pt(v4.Pt());
    set_eta(v4.Eta());
    set_phi(v4.Phi());
    set_energy(v4.E());
  }

  double deltaPhi(const Particle p2) const{
    double deltaphi = fabs(this->phi() - p2.phi());
    if(deltaphi > PI) deltaphi = 2* PI - deltaphi;
    return deltaphi;
  }
  double deltaR(const Particle p2) const{
    double deltaeta = this->eta() - p2.eta();
    return sqrt(deltaeta*deltaeta+deltaPhi(p2)*deltaPhi(p2));
  }

 private:

  float m_charge;
  float m_pt; 
  float m_eta; 
  float m_phi; 
  float m_energy; 

};

#endif
