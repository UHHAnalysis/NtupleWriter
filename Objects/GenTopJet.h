#ifndef GenTopJet_H
#define GenTopJet_H

#include "Particle.h"
#include "TLorentzVector.h"


/**
 *  @short top-jet class with only subjets, used for generator information
 *  @author Roman Kogler
 */

class GenTopJet : public Particle {
 public:
  GenTopJet(){
    m_subjets.clear();
  };

  ~GenTopJet(){
  };


  std::vector<Particle> subjets() const{return m_subjets;}
  void add_subjet(Particle p){m_subjets.push_back(p);}
  std::vector<unsigned int> genparticles_indices() const{return m_genparticles_indices;}
  void add_genparticles_index(int ind){m_genparticles_indices.push_back(ind);}

 private:
  std::vector<Particle> m_subjets;
std::vector<unsigned int> m_genparticles_indices;


};

#endif
