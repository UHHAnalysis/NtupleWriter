#include "UHHAnalysis/NtupleWriter/Objects/Particle.h"
#include "UHHAnalysis/NtupleWriter/Objects/Jet.h"
#include "UHHAnalysis/NtupleWriter/Objects/Electron.h"
#include "UHHAnalysis/NtupleWriter/Objects/Muon.h"
#include "UHHAnalysis/NtupleWriter/Objects/Tau.h"
#include "UHHAnalysis/NtupleWriter/Objects/Photon.h"
#include "UHHAnalysis/NtupleWriter/Objects/MET.h"
#include "UHHAnalysis/NtupleWriter/Objects/PrimaryVertex.h"
#include "UHHAnalysis/NtupleWriter/Objects/TopJet.h"
#include "UHHAnalysis/NtupleWriter/Objects/GenInfo.h"
#include "UHHAnalysis/NtupleWriter/Objects/GenParticle.h"

#include <vector>

namespace {
  namespace {
    Particle p;
    std::vector<Particle> ps;
    Jet jet;
    std::vector<Jet> jets;
    TopJet topjet;
    std::vector<TopJet> topjets;
    Electron ele; 
    std::vector<Electron> eles; 
    Muon mu; 
    std::vector<Muon> mus; 
    Tau tau;
    std::vector<Tau> taus; 
    Photon ph; 
    std::vector<Photon> phs; 
    MET met;
    PrimaryVertex pv;
    std::vector<PrimaryVertex> pvs; 
    GenInfo genInfo;
    GenParticle genp;
    std::vector<GenParticle> genps;
  }
}
