#ifndef TopGenJet_H
#define TopGenJet_H

#include "Particle.h"
#include "TLorentzVector.h"


class TopGenJet : public Particle {
 public:
  TopGenJet(){
    m_genparticles_indices.clear();
  };

  ~TopGenJet(){
  };

  std::vector<unsigned int> genparticles_indices() const{return m_genparticles_indices;}
  void add_genparticles_index(int ind){m_genparticles_indices.push_back(ind);}

 private:
std::vector<unsigned int> m_genparticles_indices;


};

#endif
