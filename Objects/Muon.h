#ifndef Muon_H
#define Muon_H

#include "Particle.h"



class Muon : public Particle{

 public:
  Muon(){
    m_vertex_x=0; 
    m_vertex_y=0; 
    m_vertex_z=0; 
    m_dB=0; 
    // m_particleIso=0; 
    m_neutralHadronIso=0; 
    m_chargedHadronIso=0; 
    m_trackIso=0; 
    m_photonIso=0;
    m_puChargedHadronIso=0;
    m_isGlobalMuon=false;
    m_isStandAloneMuon=false;
    m_isTrackerMuon=false;
    m_numberOfMatchedStations=0;
    m_globalTrack_chi2=0;
    m_globalTrack_ndof=0;
    m_globalTrack_d0=0;
    m_globalTrack_d0Error=0; 
    m_globalTrack_numberOfValidHits=0;  
    m_globalTrack_numberOfLostHits=0;  
    m_innerTrack_chi2=0;
    m_innerTrack_ndof=0;
    m_innerTrack_d0=0;
    m_innerTrack_d0Error=0; 
    m_innerTrack_numberOfValidHits=0;  
    m_innerTrack_numberOfLostHits=0;  
    m_innerTrack_trackerLayersWithMeasurement=0;
    m_innerTrack_numberOfValidPixelHits=0;
    m_outerTrack_chi2=0;
    m_outerTrack_ndof=0;
    m_outerTrack_d0=0;
    m_outerTrack_d0Error=0; 
    m_outerTrack_numberOfValidHits=0;  
    m_outerTrack_numberOfLostHits=0; 
  };

  ~Muon(){
  };

  float vertex_x(){return m_vertex_x;} 
  float vertex_y(){return m_vertex_y;} 
  float vertex_z(){return m_vertex_z;} 
  float dB(){return m_dB;} 
  //float particleIso(){return m_particleIso;} 
  float neutralHadronIso(){return m_neutralHadronIso;} 
  float chargedHadronIso(){return m_chargedHadronIso;} 
  float trackIso(){return m_trackIso;} 
  float photonIso(){return m_photonIso;}
  float puChargedHadronIso(){return m_puChargedHadronIso;}
  bool isGlobalMuon(){return m_isGlobalMuon;}
  bool isStandAloneMuon(){return m_isStandAloneMuon;}
  bool isTrackerMuon(){return m_isTrackerMuon;}
  int numberOfMatchedStations(){return m_numberOfMatchedStations;}
  float globalTrack_chi2(){return m_globalTrack_chi2;}
  float globalTrack_ndof(){return m_globalTrack_ndof;}
  float globalTrack_d0(){return m_globalTrack_d0;}
  float globalTrack_d0Error(){return m_globalTrack_d0Error;} 
  unsigned short globalTrack_numberOfValidHits(){return m_globalTrack_numberOfValidHits;}  
  unsigned short globalTrack_numberOfLostHits(){return m_globalTrack_numberOfLostHits;}  
  float innerTrack_chi2(){return m_innerTrack_chi2;}
  float innerTrack_ndof(){return m_innerTrack_ndof;}
  float innerTrack_d0(){return m_innerTrack_d0;}
  float innerTrack_d0Error(){return m_innerTrack_d0Error;} 
  unsigned short innerTrack_numberOfValidHits(){return m_innerTrack_numberOfValidHits;}  
  unsigned short innerTrack_numberOfLostHits(){return m_innerTrack_numberOfLostHits;}  
  unsigned short innerTrack_trackerLayersWithMeasurement(){return m_innerTrack_trackerLayersWithMeasurement;}
  unsigned short innerTrack_numberOfValidPixelHits(){return m_innerTrack_numberOfValidPixelHits;}
  float outerTrack_chi2(){return m_outerTrack_chi2;}
  float outerTrack_ndof(){return m_outerTrack_ndof;}
  float outerTrack_d0(){return m_outerTrack_d0;}
  float outerTrack_d0Error(){return m_outerTrack_d0Error;} 
  unsigned short outerTrack_numberOfValidHits(){return m_outerTrack_numberOfValidHits;}  
  unsigned short outerTrack_numberOfLostHits(){return m_outerTrack_numberOfLostHits;}  

  void set_vertex_x(float x){m_vertex_x=x;} 
  void set_vertex_y(float x){m_vertex_y=x;} 
  void set_vertex_z(float x){m_vertex_z=x;} 
  void set_dB(float x){m_dB=x;} 
  //void set_particleIso(float x){m_particleIso=x;} 
  void set_neutralHadronIso(float x){m_neutralHadronIso=x;} 
  void set_chargedHadronIso(float x){m_chargedHadronIso=x;} 
  void set_trackIso(float x){m_trackIso=x;} 
  void set_photonIso(float x){m_photonIso=x;}
  void set_puChargedHadronIso(float x){m_puChargedHadronIso=x;}
  void set_isGlobalMuon(bool x){m_isGlobalMuon=x;}
  void set_isStandAloneMuon(bool x){m_isStandAloneMuon=x;}
  void set_isTrackerMuon(bool x){m_isTrackerMuon=x;}
  void set_numberOfMatchedStations(int x){m_numberOfMatchedStations=x;}
  void set_globalTrack_chi2(float x){m_globalTrack_chi2=x;}
  void set_globalTrack_ndof(float x){m_globalTrack_ndof=x;}
  void set_globalTrack_d0(float x){m_globalTrack_d0=x;}
  void set_globalTrack_d0Error(float x){m_globalTrack_d0Error=x;} 
  void set_globalTrack_numberOfValidHits(unsigned short x){m_globalTrack_numberOfValidHits=x;}  
  void set_globalTrack_numberOfLostHits(unsigned short x){m_globalTrack_numberOfLostHits=x;}  
  void set_innerTrack_chi2(float x){m_innerTrack_chi2=x;}
  void set_innerTrack_ndof(float x){m_innerTrack_ndof=x;}
  void set_innerTrack_d0(float x){m_innerTrack_d0=x;}
  void set_innerTrack_d0Error(float x){m_innerTrack_d0Error=x;} 
  void set_innerTrack_numberOfValidHits(unsigned short x){m_innerTrack_numberOfValidHits=x;}  
  void set_innerTrack_numberOfLostHits(unsigned short x){m_innerTrack_numberOfLostHits=x;}  
  void set_innerTrack_trackerLayersWithMeasurement(unsigned short x){m_innerTrack_trackerLayersWithMeasurement=x;}
  void set_innerTrack_numberOfValidPixelHits(unsigned short x){m_innerTrack_numberOfValidPixelHits=x;}
  void set_outerTrack_chi2(float x){m_outerTrack_chi2=x;}
  void set_outerTrack_ndof(float x){m_outerTrack_ndof=x;}
  void set_outerTrack_d0(float x){m_outerTrack_d0=x;}
  void set_outerTrack_d0Error(float x){m_outerTrack_d0Error=x;} 
  void set_outerTrack_numberOfValidHits(unsigned short x){m_outerTrack_numberOfValidHits=x;}  
  void set_outerTrack_numberOfLostHits(unsigned short x){m_outerTrack_numberOfLostHits=x;}  

  float relIso(){
    return ( m_chargedHadronIso + std::max( 0.0, m_neutralHadronIso + m_photonIso - 0.5*m_puChargedHadronIso ) ) / pt();
  }

 private:
  float m_vertex_x; 
  float m_vertex_y; 
  float m_vertex_z; 
  float m_dB; 
  //float m_particleIso; 
  float m_neutralHadronIso; 
  float m_chargedHadronIso; 
  float m_trackIso; 
  float m_photonIso;
  float m_puChargedHadronIso;
  bool m_isGlobalMuon;
  bool m_isStandAloneMuon;
  bool m_isTrackerMuon;
  int m_numberOfMatchedStations;
  float m_globalTrack_chi2;
  float m_globalTrack_ndof;
  float m_globalTrack_d0;
  float m_globalTrack_d0Error; 
  unsigned short m_globalTrack_numberOfValidHits;  
  unsigned short m_globalTrack_numberOfLostHits;  
  float m_innerTrack_chi2;
  float m_innerTrack_ndof;
  float m_innerTrack_d0;
  float m_innerTrack_d0Error; 
  unsigned short m_innerTrack_numberOfValidHits;  
  unsigned short m_innerTrack_numberOfLostHits;  
  unsigned short m_innerTrack_trackerLayersWithMeasurement;
  unsigned short m_innerTrack_numberOfValidPixelHits;
  float m_outerTrack_chi2;
  float m_outerTrack_ndof;
  float m_outerTrack_d0;
  float m_outerTrack_d0Error; 
  unsigned short m_outerTrack_numberOfValidHits;  
  unsigned short m_outerTrack_numberOfLostHits;  
  


};

#endif
