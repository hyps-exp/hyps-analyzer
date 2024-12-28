// -*- C++ -*-

#include "CFTFiberCluster.hh"

#include <cmath>
#include <string>
#include <limits>

#include <TMath.h>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "FiberHit.hh"
#include "FuncName.hh"

namespace
{
const Bool_t reject_nan = false;
}

//_____________________________________________________________________________
CFTFiberCluster::CFTFiberCluster(const HodoHC& cont,
                           const index_t& index)
  : HodoCluster(cont, index),
    m_total_adc_hi(),
    m_max_adc_hi(),  
    m_total_adc_low(),
    m_max_adc_low(),  
    m_total_mip_hi(),
    m_max_mip_hi(),  
    m_total_mip_low(),
    m_max_mip_low(),  
    m_total_de_hi(),
    m_max_de_hi(),  
    m_total_de_low(),
    m_max_de_low(),  
    m_max_ctime(),
    m_mean_x(),
    m_mean_y(),  
    m_mean_phi(),
    m_mean_r(),  
    m_mean_z0(),
    m_slope(),
    m_z(),
    m_mean_x_cor(),
    m_mean_y_cor(),
    m_mean_phi_cor(),
    m_x_cal(),
    m_y_cal(),
    m_z_cal(),
    m_phi_cal(),
    m_pathlength()
{
  Calculate();
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
CFTFiberCluster::~CFTFiberCluster()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
CFTFiberCluster::Calculate()
{
  if(!m_is_good)
    return;

  m_total_adc_hi = 0.;
  m_max_adc_hi = 0.;  
  m_total_adc_low = 0.;
  m_max_adc_low = 0.;  
  m_total_mip_hi = 0.;
  m_max_mip_hi = 0.;  
  m_total_mip_low = 0.;
  m_max_mip_low = 0.;  
  m_total_de_hi = 0.;
  m_max_de_hi = 0.;  
  m_total_de_low = 0.;
  m_max_de_low = 0.;

  m_de = 0.;
  
  m_max_ctime = 0.;
  m_mean_x = TMath::QuietNaN();
  m_mean_y = TMath::QuietNaN();
  m_mean_phi = TMath::QuietNaN();
  m_mean_r = TMath::QuietNaN();
  m_mean_z0 = TMath::QuietNaN();
  m_slope = TMath::QuietNaN();
  m_z = 0.;
  m_mean_x_cor = 0.;
  m_mean_y_cor = 0.;
  m_mean_phi_cor = 0.;
  m_x_cal = 0.;
  m_y_cal = 0.;
  m_z_cal = 0.;
  m_phi_cal = 0.;
  m_pathlength = 0.;
  
  m_tot           = 0.;
  m_mean_position = 0.;

  for(Int_t i=0; i<m_cluster_size; ++i){
    const auto& hit = dynamic_cast<CFTFiberHit*>(m_hit_container[i]);
    const auto& index = m_index[i];
    m_tot = TMath::Max(m_tot, hit->MeanTOT(index));
    // m_tot += hit->TOT(index);
    m_mean_position += hit->Position();

    Double_t adccorHi  = hit->GetAdcCorHigh();
    if (adccorHi>0) {
      m_total_adc_hi += adccorHi;
      if (adccorHi > m_max_adc_hi)
	m_max_adc_hi = adccorHi;
    }
    Double_t adccorLow  = hit->GetAdcCorLow();
    if (adccorLow>0) {
      m_total_adc_low += adccorLow;
      if (adccorLow > m_max_adc_low)
	m_max_adc_low = adccorLow;
    }

    Double_t dE_Low  = hit->DeltaELowGain();
    if (dE_Low>0) {
      m_total_de_low += dE_Low;
      m_de += dE_Low;      
      if (dE_Low > m_max_de_low)
	m_max_de_low = dE_Low;
    }
  }
  
  m_mean_position /= Double_t(m_cluster_size);


  Int_t n_true_hit=0;
  for(Int_t i=0; i<m_cluster_size; ++i){
    const auto& hit = dynamic_cast<CFTFiberHit*>(m_hit_container[i]);
    Double_t dE_Low  = hit->DeltaELowGain();

    if (dE_Low > m_max_de_low*0.3) {
      n_true_hit++;
      if (n_true_hit == 1) {
	m_mean_x   = hit->GetX();
	m_mean_y   = hit->GetY();
	m_mean_phi = hit->GetPhi();
	m_mean_r   = hit->GetR();
	m_mean_z0  = hit->GetZ0();
	m_slope    = hit->GetSlope();
      } else {
	m_mean_x   += hit->GetX();
	m_mean_y   += hit->GetY();
	m_mean_phi += hit->GetPhi();
	m_mean_r   += hit->GetR();
	m_mean_z0  += hit->GetZ0();
      }
    }
  }
  if (n_true_hit>0) {
    m_mean_x   /= double(n_true_hit);
    m_mean_y   /= double(n_true_hit);
    m_mean_phi /= double(n_true_hit);
    m_mean_r   /= double(n_true_hit);
    m_mean_z0  /= double(n_true_hit);
  }
  
}
