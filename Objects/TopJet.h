#ifndef TopJet_H
#define TopJet_H

#include "Jet.h"

/**
 *  @short top-jet class
 *  @author Thomas Peiffer
 */

class TopJet : public Jet{
 public:
  TopJet(){
    m_subjets.clear();
    m_btagsub_combinedSecondaryVertex.clear();
    m_flavorsub.clear();
  };
  ~TopJet(){
  };

  std::vector<Particle> subjets() const{return m_subjets;}
  std::vector<float> btagsub_combinedSecondaryVertex() const{return m_btagsub_combinedSecondaryVertex;}
  std::vector<int> flavorsub() const{return m_flavorsub;}

  void add_subjet(Particle p){m_subjets.push_back(p);}
  void add_btagsub_combinedSecondaryVertex(float discr){m_btagsub_combinedSecondaryVertex.push_back(discr);}
  void add_flavorsub(int flav){m_flavorsub.push_back(flav);}

 private:
  std::vector<Particle> m_subjets;
  std::vector<float> m_btagsub_combinedSecondaryVertex;
  std::vector<int> m_flavorsub;
};

#endif
