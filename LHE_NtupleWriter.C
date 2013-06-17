

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include "TString.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include "Objects/Particle.h"
#include "Objects/GenInfo.h"
#include "Objects/GenParticle.h"

#pragma link C++ class Particle+;
#pragma link C++ class std::vector<Particle>+;
#pragma link C++ class GenParticle+;
#pragma link C++ class std::vector<GenParticle>+;

using namespace std;

struct particle{

  int pdg_id;
  int mo_1;
  int mo_2;
  double px;
  double py;
  double pz;
  double E;
  double M;
  int spin;
};



void LHE_NtupleWriter(TString infilename, TString outfilename){

  ifstream file;

  file.open(infilename);  
  

  TFile* outfile = new TFile(outfilename,"RECREATE");
  outfile->cd();
  TTree* tr = new TTree("AnalysisTree","AnalysisTree");

  int run=1;
  int luminosityBlock=1;
  int event=0;
  bool isRealData=false;
  float rho=0;
  std::vector<GenParticle> genps;

  tr->Branch("run",&run);
  tr->Branch("event",&event);
  tr->Branch("luminosityBlock",&luminosityBlock);
  tr->Branch("isRealData",&isRealData);
  tr->Branch("rho",&rho);

  tr->Branch("GenParticles","std::vector<GenParticle>", &genps);

  while(!file.eof()){

    TString word;
    file>>word;
    //std::cout <<word <<std::endl;

    bool valid_event=true;

    if(word=="<event>"){
      event++;
      vector<particle> particles;
      
      //komische Zahlen am Anfang vom Event
      for(int i=0; i<6;++i){
	 file>>word;
      }
      file>>word;
      while(word!="</event>" && word!="#" && word !="#IHPRO"){

	particle p;

	p.pdg_id = atoi((const char*)word);

	file>>word;
	file>>p.mo_1;
	file>>p.mo_2;
	file>>word;
	file>>word;
	file>>p.px;
	file>>p.py;	
	file>>p.pz;
	file>>p.E;
	file>>p.M;
	file>>word;
	file>>word;
	p.spin = atoi((const char*)word);

	particles.push_back(p);

	file>>word;
      }

      if(event%1000==0){
	std::cout << "-------------------- Event " << event << "  ------------------------------------"<<std::endl;
      }

      genps.clear();

      for(size_t j=0; j<particles.size(); ++j){

	//remove events with pdg_id=0
	if(particles[j].pdg_id==0) {
	  std::cout << "Found strange particle with pdgId=0 in event " << event <<" at position " << j << " ---- event will be skipped" << std::endl;
	  valid_event=false;
	  break;
	}

	TLorentzVector v4;
	v4.SetPxPyPzE(particles[j].px,particles[j].py,particles[j].pz,particles[j].E);

	GenParticle genp;
	genp.set_charge(0);
	if(particles[j].pdg_id==2 ||  particles[j].pdg_id==4 || particles[j].pdg_id==6)
	  genp.set_charge(2./3.);
	if(particles[j].pdg_id==-2 ||  particles[j].pdg_id==-4 || particles[j].pdg_id==-6)
	  genp.set_charge(-2./3.);
	if(particles[j].pdg_id==1 ||  particles[j].pdg_id==3 || particles[j].pdg_id==5)
	  genp.set_charge(1./3.);
	if(particles[j].pdg_id==-1 ||  particles[j].pdg_id==-3 || particles[j].pdg_id==-5)
	  genp.set_charge(-1./3.);
	if(particles[j].pdg_id==11 ||  particles[j].pdg_id==13 || particles[j].pdg_id==15)
	  genp.set_charge(-1);
	if(particles[j].pdg_id==-11 ||  particles[j].pdg_id==-13 || particles[j].pdg_id==-15)
	  genp.set_charge(1);
	if(particles[j].pdg_id==24)
	  genp.set_charge(1);
	if(particles[j].pdg_id==-24)
	  genp.set_charge(-1);

	genp.set_pt(v4.Pt());
	if(v4.Pt()>0) genp.set_eta(v4.Eta());
	else if(v4.Pz()>=0)
	  genp.set_eta(10e10);
	else 
	  genp.set_eta(-10e10);

	genp.set_phi(v4.Phi());
	genp.set_energy(v4.E());
	genp.set_index(j+1);
	genp.set_status( 3);
	genp.set_pdgId( particles[j].pdg_id );
	genp.set_spin(particles[j].spin);
	
	genp.set_mother1(particles[j].mo_1);
	genp.set_mother2(particles[j].mo_2);
	genp.set_daughter1(-1);
	genp.set_daughter2(-1);

	for(size_t k=j+1; k<particles.size(); ++k){
	  if(particles[k].mo_1==j+1 && genp.daughter1()<0) {
	    genp.set_daughter1(k+1);
	  }
	  else if(particles[k].mo_1==j+1 && genp.daughter2()<0) {
	    genp.set_daughter2(k+1);
	  }
	  if(particles[k].mo_2==j+1 && particles[k].mo_2!=particles[k].mo_1 && genp.daughter1()<0) {
	    genp.set_daughter1(k+1);
	  } 
	  else if(particles[k].mo_2==j+1 && particles[k].mo_2!=particles[k].mo_1 && genp.daughter2()<0) {
	    genp.set_daughter2(k+1);
	  }
	  if(genp.daughter1()>=0 && genp.daughter2()>=0) break;
	  
	}
	genps.push_back(genp);

      }
      if(valid_event)
	tr->Fill();

    }
  }
  tr->Write();
  outfile->Close();
  file.close();
  
}
