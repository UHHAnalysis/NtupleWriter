#ifndef TopJet_H
#define TopJet_H

#include "Jet.h"

class TopJet : public Jet{
 public:
  TopJet(){
    m_subjets.clear();
  };
  ~TopJet(){
  };

  std::vector<Particle> subjets(){return m_subjets;}

  void add_subjet(Particle p){m_subjets.push_back(p);}

 private:
  std::vector<Particle> m_subjets;

};

#endif
