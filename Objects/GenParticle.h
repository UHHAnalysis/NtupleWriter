#ifndef GenParticle_H
#define GenParticle_H

#include "Particle.h"

/**
 *  @short generator particle class
 *  @author Thomas Peiffer
 */

class GenParticle : public Particle{
 public:
  GenParticle(){
    m_pdgId=0;
    m_status=0;
    m_index=0;
    m_mother1=0;
    m_mother2=0;
    m_daughter1=0;
    m_daughter2=0;
    m_spin=0;
  };
  ~GenParticle(){
  };

  int pdgId() const{return m_pdgId;}
  int status() const{return m_status;}
  int index() const{return m_index;}
  int mother1() const{return m_mother1;}
  int mother2() const{return m_mother2;}
  int daughter1() const{return m_daughter1;}
  int daughter2() const{return m_daughter2;}
  int spin() const{return m_spin;}

  //return mother 1 or 2 (ind<=1 or ind>=2)
  GenParticle* mother(std::vector<GenParticle> *gplist, int ind=1){
    for(unsigned int i=0; i< gplist->size(); ++i){
      if(ind<=1){
	if(this->m_mother1 == gplist->at(i).index()){
	  return &(gplist->at(i));
	}
      }
      else{
	if(this->m_mother2 == gplist->at(i).index()){
	  return &(gplist->at(i));
	}	
      }
    }
    //std::cout << "WARNING: Mother " << ind << " not found in list of GenParticles" << std::endl;
    return 0;
  }
  //return daughter 1 or 2 (ind<=1 or ind>=2)
  GenParticle* daughter(std::vector<GenParticle> *gplist, int ind=1){
    for(unsigned int i=0; i< gplist->size(); ++i){
      if(ind<=1){
	if(this->m_daughter1 == gplist->at(i).index()){
	  return &(gplist->at(i));
	}
      }
      else{
	if(this->m_daughter2 == gplist->at(i).index()){
	  return &(gplist->at(i));
	}	
      }
    }
    //std::cout << "WARNING: Daughter " << ind << " not found in list of GenParticles" << std::endl;
    return 0;
  }

  void set_pdgId(int x){  m_pdgId=x;}
  void set_status(int x){  m_status=x;}
  void set_index(int x){  m_index=x;}
  void set_mother1(int x){  m_mother1=x;}
  void set_mother2(int x){  m_mother2=x;}
  void set_daughter1(int x){  m_daughter1=x;}
  void set_daughter2(int x){  m_daughter2=x;}
  void set_spin(int x){  m_spin=x;}

 private:
  int m_pdgId;
  int m_status;
  int m_index;

  int m_mother1;
  int m_mother2;
  int m_daughter1;
  int m_daughter2;
  int m_spin;
 

};

#endif
