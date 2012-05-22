#ifndef MET_H
#define MET_H

class MET{
   
 public:
  MET(){
    m_pt=0; 
    m_phi=0; 
    m_mEtSig=0;

  };

  ~MET(){
  };


  float pt(){return m_pt;}
  float phi(){return m_phi;}
  float mEtSig(){return m_mEtSig;}

  void set_pt(float pt){m_pt=pt;}  
  void set_phi(float phi){m_phi=phi;}
  void set_mEtSig(float mEtSig){m_mEtSig=mEtSig;}

 private:
  float m_pt; 
  float m_phi; 
  float m_mEtSig; 

};

#endif
