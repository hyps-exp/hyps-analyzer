// -*- C++ -*-

#ifndef CFT_FIBER_CLUSTER_HH
#define CFT_FIBER_CLUSTER_HH

#include <vector>

#include <TString.h>

#include "CFTFiberHit.hh"
#include "HodoCluster.hh"

//_____________________________________________________________________________
class CFTFiberCluster : public HodoCluster
{
public:
  static const TString& ClassName();
  CFTFiberCluster(const HodoHC& cont,
               const index_t& index);
  virtual ~CFTFiberCluster();

private:
  CFTFiberCluster(const CFTFiberCluster&);
  CFTFiberCluster& operator =(const CFTFiberCluster&);

protected:
  Double_t m_total_adc_hi;
  Double_t m_max_adc_hi;  
  Double_t m_total_adc_low;
  Double_t m_max_adc_low;  
  Double_t m_total_mip_hi;
  Double_t m_max_mip_hi;  
  Double_t m_total_mip_low;
  Double_t m_max_mip_low;  
  Double_t m_total_de_hi;
  Double_t m_max_de_hi;  
  Double_t m_total_de_low;
  Double_t m_max_de_low;  

  Double_t m_max_ctime;

  Double_t m_mean_x;
  Double_t m_mean_y;  
  Double_t m_mean_phi;
  Double_t m_mean_r;  
  Double_t m_mean_z0;
  Double_t m_slope;

  Double_t m_z;

  Double_t m_mean_x_cor;
  Double_t m_mean_y_cor;
  Double_t m_mean_phi_cor;

  Double_t m_x_cal;
  Double_t m_y_cal;
  Double_t m_z_cal;
  Double_t m_phi_cal;

  Double_t m_pathlength;
  
public:
  Double_t MaxDeltaE() const { return m_max_de_low; }
  Double_t TotalDeltaE() const { return m_total_de_low; }  
  Double_t MeanX() const  { return m_mean_x; }
  Double_t MeanY() const  { return m_mean_y; }
  Double_t MeanPhi() const  { return m_mean_phi; }
  Double_t MeanR() const  { return m_mean_r; }
  Double_t MeanZ0() const  { return m_mean_z0; }
  Double_t Slope() const  { return m_slope; }          
protected:
  void Calculate();
};

//_____________________________________________________________________________
inline const TString&
CFTFiberCluster::ClassName()
{
  static TString s_name("CFTFiberCluster");
  return s_name;
}

#endif
