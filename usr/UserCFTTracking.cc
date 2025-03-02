// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <sstream>

#include <UnpackerManager.hh>

// #include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "CFTFiberHit.hh"
#include "CFTFiberCluster.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "HodoHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "DCAnalyzer.hh"
#include "CFTLocalTrack.hh"
#include "HypsLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "DCGeomMan.hh"
#include "CFTPedCorMan.hh"
#include "CFTPosParamMan.hh"

// #define TimeCut    1 // in cluster analysis
#define FHitBranch 0 // make FiberHit branches (becomes heavy)
#define HodoHitPos 0

namespace
{
using namespace root;
const auto qnan = TMath::QuietNaN();
auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto& gRM       = RMAnalyzer::GetInstance();
auto& gUser     = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t spill;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Int_t cftadc_h[NumOfPlaneCFT][NumOfSegCFT_PHI4];
  Int_t cftadc_l[NumOfPlaneCFT][NumOfSegCFT_PHI4];
  Int_t cfttdc[NumOfPlaneCFT][NumOfSegCFT_PHI4];

  Int_t cftadc_cor_h[NumOfPlaneCFT][NumOfSegCFT_PHI4];
  Int_t cftadc_cor_l[NumOfPlaneCFT][NumOfSegCFT_PHI4];


  Int_t    ntCFT;
  Double_t u0[MaxDepth];
  Double_t v0[MaxDepth];
  Double_t x0[MaxDepth];
  Double_t y0[MaxDepth];
  Double_t dist[MaxDepth];
  Double_t distMeanx[MaxDepth];
  Double_t distMeany[MaxDepth];
  Double_t distMeanz[MaxDepth];
  Double_t cos_16layer[MaxDepth];

  Int_t    ntCFT_16layer;
  Int_t    nhXY_16layer[MaxDepth];
  Int_t    nhZ_16layer[MaxDepth];
  Double_t chisqrXY_16layer[MaxDepth];
  Double_t chisqrZ_16layer[MaxDepth];
  Double_t u0_16layer[MaxDepth];
  Double_t v0_16layer[MaxDepth];
  Double_t x0_16layer[MaxDepth];
  Double_t y0_16layer[MaxDepth];
  Double_t z0_y_16layer[MaxDepth];
  Double_t z0_x_16layer[MaxDepth];
  Double_t dist_16layer[MaxDepth];
  Double_t costheta[MaxDepth];
  Double_t theta_16layer[MaxDepth];
  Double_t cost_16layer[MaxDepth];
  Int_t    nhPhi_16layer[4];
  Double_t hitMeanSegPhi_16layer[4][MaxDepth];
  Double_t phi_16layer[4][MaxDepth];
  Double_t phicor_16layer[4][MaxDepth];
  Double_t phical_16layer[4][MaxDepth];
  Double_t resPhi_16layer[4][MaxDepth];
  Double_t calZPhi_16layer[4][MaxDepth];
  Int_t    nhU_16layer[4];
  Double_t hitMeanSegU_16layer[4][MaxDepth];
  Double_t phiU_16layer[4][MaxDepth];
  Double_t calZU_16layer[4][MaxDepth];
  Double_t resZU_16layer[4][MaxDepth];


  ////////// Normalized

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum     = 0;
  spill     = 0;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for( int it=0; it<NumOfPlaneCFT; ++it ){
    for(int m = 0; m<NumOfSegCFT_PHI4; ++m){
      cftadc_h[it][m] = qnan;
      cftadc_l[it][m] = qnan;
      cfttdc[it][m] = qnan;

      cftadc_cor_h[it][m] = qnan;
      cftadc_cor_l[it][m] = qnan;
    }
  }

  ntCFT = 0;
  ntCFT_16layer = 0;
  for (Int_t j=0; j<4; j++) {
    nhPhi_16layer[j] = 0;
    nhU_16layer[j] = 0;
  }

  for (Int_t j=0; j<MaxDepth; j++) {
    x0[j] = qnan;
    y0[j] = qnan;
    u0[j] = qnan;
    v0[j] = qnan;
    dist[j] = qnan;
    distMeanx[j] = qnan;
    distMeany[j] = qnan;
    distMeanz[j] = qnan;
    cos_16layer[j] = qnan;
    nhXY_16layer[j] = 0;
    nhZ_16layer[j] = 0;
    chisqrXY_16layer[j] = qnan;
    chisqrZ_16layer[j] = qnan;
    x0_16layer[j] = qnan;
    y0_16layer[j] = qnan;
    u0_16layer[j] = qnan;
    v0_16layer[j] = qnan;
    z0_y_16layer[j] = qnan;
    z0_x_16layer[j] = qnan;
    dist_16layer[j] = qnan;
    costheta[j] = qnan;
    theta_16layer[j] = qnan;
    cost_16layer[j] = qnan;
    for (Int_t i=0; i<4; i++) {
      hitMeanSegPhi_16layer[i][j] = qnan;
      phi_16layer[i][j] = qnan;
      phicor_16layer[i][j] = qnan;
      phical_16layer[i][j] = qnan;
      resPhi_16layer[i][j] = qnan;
      calZPhi_16layer[i][j] = qnan;

      hitMeanSegU_16layer[i][j] = qnan;
      phiU_16layer[i][j] = qnan;
      calZU_16layer[i][j] = qnan;
      resZU_16layer[i][j] = qnan;
    }
  }
}

//_____________________________________________________________________________
struct Dst
{
  Int_t evnum;
  Int_t spill;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Int_t    nhBh2;
  Int_t    csBh2[NumOfSegBH2*MaxDepth];
  Double_t Bh2Seg[NumOfSegBH2*MaxDepth];
  Double_t tBh2[NumOfSegBH2*MaxDepth];
  Double_t t0Bh2[NumOfSegBH2*MaxDepth];
  Double_t dtBh2[NumOfSegBH2*MaxDepth];
  Double_t deBh2[NumOfSegBH2*MaxDepth];
  Double_t posBh2[NumOfSegBH2*MaxDepth];

  Int_t    nhT0;
  Int_t    csT0[NumOfSegT0*MaxDepth];
  Double_t T0Seg[NumOfSegT0*MaxDepth];
  Double_t tT0[NumOfSegT0*MaxDepth];

  // for HodoParam
  Double_t t0la[NumOfSegT0];
  Double_t t0lt[NumOfSegT0][MaxDepth];
  Double_t t0ra[NumOfSegT0];
  Double_t t0rt[NumOfSegT0][MaxDepth];
  Double_t ltT0Seg[NumOfSegT0][MaxDepth];
  Double_t rtT0Seg[NumOfSegT0][MaxDepth];

  Int_t    nhSac;
  Int_t    csSac[NumOfSegSAC*MaxDepth];
  Double_t SacSeg[NumOfSegSAC*MaxDepth];
  Double_t tSac[NumOfSegSAC*MaxDepth];
  Double_t deSac[NumOfSegSAC*MaxDepth];

  // Time0
  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  Int_t    nhTof;
  Int_t    csTof[NumOfSegTOF*MaxDepth];
  Double_t TofSeg[NumOfSegTOF*MaxDepth];
  Double_t tTof[NumOfSegTOF*MaxDepth];
  Double_t dtTof[NumOfSegTOF*MaxDepth];
  Double_t deTof[NumOfSegTOF*MaxDepth];

  // for HodoParam
  Double_t tofua[NumOfSegTOF];
  Double_t tofut[NumOfSegTOF][MaxDepth];
  Double_t tofda[NumOfSegTOF];
  Double_t tofdt[NumOfSegTOF][MaxDepth];
  Double_t utTofSeg[NumOfSegTOF][MaxDepth];
  Double_t dtTofSeg[NumOfSegTOF][MaxDepth];
  Double_t udeTofSeg[NumOfSegTOF];
  Double_t ddeTofSeg[NumOfSegTOF];

  void clear();
};

//_____________________________________________________________________________
void
Dst::clear()
{
  nhT0     = 0;
  nhSac    = 0;
  nhBh2    = 0;
  nhTof    = 0;
  evnum    = 0;
  spill    = 0;
  Time0Seg = qnan;
  deTime0  = qnan;
  Time0    = qnan;
  CTime0   = qnan;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<NumOfSegBH2; ++it){
    for(Int_t m=0; m<MaxDepth; ++m){
      csBh2[MaxDepth*it + m]  = 0;
      Bh2Seg[MaxDepth*it + m] = qnan;
      tBh2[MaxDepth*it + m]   = qnan;
      t0Bh2[MaxDepth*it + m]  = qnan;
      dtBh2[MaxDepth*it + m]  = qnan;
      deBh2[MaxDepth*it + m]  = qnan;
      posBh2[MaxDepth*it + m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegT0; it++){
    t0la[it]  = qnan;
    t0ra[it]  = qnan;
    T0Seg[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      t0lt[it][m]        = qnan;
      t0rt[it][m]        = qnan;
      tT0[MaxDepth*it+m] = qnan;
      ltT0Seg[it][m]     = qnan;
      rtT0Seg[it][m]     = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegSAC; it++){
    SacSeg[it] = qnan;
    deSac[it]  = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tSac[MaxDepth*it+m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegTOF; it++){
    tofua[it]     = qnan;
    tofda[it]     = qnan;
    udeTofSeg[it] = qnan;
    ddeTofSeg[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tofut[it][m]            = qnan;
      tofdt[it][m]            = qnan;
      utTofSeg[it][m]         = qnan;
      dtTofSeg[it][m]         = qnan;
      csTof[MaxDepth*it + m]  = 0;
      TofSeg[MaxDepth*it + m] = qnan;
      tTof[MaxDepth*it + m]   = qnan;
      dtTof[MaxDepth*it + m]  = qnan;
      deTof[MaxDepth*it + m]  = qnan;
    }
  }

}


//_____________________________________________________________________________
namespace root
{
Event  event;
Dst    dst;
TH1*   h[MaxHist];
TTree* tree;
// TTree* hodo;
enum eDetHid {
  BH1Hid  = 10000,
  BH2Hid  = 20000,
  T0Hid   = 40000,
  SACHid  = 30000,
  TOFHid  = 60000,
};
}

//_____________________________________________________________________________
Bool_t
ProcessingBegin()
{
  event.clear();
  dst.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingNormal()
{
  static const auto MinTdcCFT = gUser.GetParameter("TdcCFT", 0);
  static const auto MaxTdcCFT = gUser.GetParameter("TdcCFT", 1);
  static const auto MinTimeCFT = gUser.GetParameter("TimeCFT", 0);
  static const auto MaxTimeCFT = gUser.GetParameter("TimeCFT", 1);
  static const auto MinAdcCFT = gUser.GetParameter("AdcCFT", 0);
  static const auto MaxAdcCFT = gUser.GetParameter("AdcCFT", 1);
  //static const auto MinTdcT0  = gUser.GetParameter("TdcT0", 0);
  //static const auto MaxTdcT0  = gUser.GetParameter("TdcT0", 1);
  //static const auto MinTdcSAC = gUser.GetParameter("TdcSAC", 0);
  //static const auto MaxTdcSAC = gUser.GetParameter("TdcSAC", 1);
  //static const auto MinTdcTOF = gUser.GetParameter("TdcTOF", 0);
  //static const auto MaxTdcTOF = gUser.GetParameter("TdcTOF", 1);
#if HodoHitPos
  static const auto PropVelBH2 = gUser.GetParameter("PropagationBH2");
#endif

  RawData rawData;
  HodoAnalyzer hodoAna(rawData);
  DCAnalyzer   DCAna(rawData);

  gRM.Decode();

  // event.evnum = gRM.EventNumber();
  event.evnum = gUnpacker.get_event_number();
  event.spill = gRM.SpillNumber();
  dst.evnum   = gRM.EventNumber();
  dst.spill   = gRM.SpillNumber();

  HF1(1, 0);


  //*****************  RawData  *****************

  // Trigger Flag
  rawData.DecodeHits("TFlag");
  std::bitset<NumOfSegTrig> trigger_flag;
  for(const auto& hit: rawData.GetHodoRawHC("TFlag")){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc();
    if(tdc > 0){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      dst.trigpat[trigger_flag.count()] = seg;
      dst.trigflag[seg] = tdc;
      trigger_flag.set(seg);
    }
  }

  if(trigger_flag[trigger::kSpillOnEnd] || trigger_flag[trigger::kSpillOffEnd])
    return true;

  HF1(1, 1);

  rawData.DecodeHits("CFT");
  // CFT
  {
    const auto& cont = rawData.GetHodoRawHC("CFT");
    const auto& U = HodoRawHit::kUp;
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t plane = hit->PlaneId();
      Int_t seg   = hit->SegmentId();

      Int_t NhitAH = hit->GetSizeAdcHigh();
      Int_t NhitAL = hit->GetSizeAdcLow();
      Int_t NhitT  = hit->GetSizeTdcUp();

      bool flag_tdc = false;
      for(Int_t m = 0; m<NhitT; ++m){
	Int_t tdc_l = hit->GetTdcLeading(U, m);
	Int_t tdc_t = hit->GetTdcTrailing(U, m);
	HF2 (1000*(plane+1)+100, seg, tdc_l);  //TDC Nhits
	HF2 (1000*(plane+1)+101, seg, tdc_t);

	Int_t width = tdc_l - tdc_t;
	HF2 (1000*(plane+1)+103, seg, width);


	if (tdc_l>MinTdcCFT && tdc_l < MaxTdcCFT) {
	  flag_tdc = true;
	  event.cfttdc[plane][seg] = tdc_l;
	}
      }
      if (flag_tdc) {
	HF1 (1000*(plane+1)+102, seg);
      }


      //ADC Hi
      for(Int_t m = 0; m<NhitAH; ++m){
	Int_t adcH = hit->GetAdcHigh();
	HF2 (1000*(plane+1)+104, seg, adcH);
	event.cftadc_h[plane][seg] = adcH;

	if (flag_tdc) {
	  HF2 (1000*(plane+1)+106, seg, adcH);
	}
      }

      //ADC Low
      for(Int_t m = 0; m<NhitAL; ++m){
	Int_t adcL = hit->GetAdcLow();
	HF2 (1000*(plane+1)+105, seg, adcL);
	event.cftadc_l[plane][seg] = adcL;
      }
    }
  }

  hodoAna.DecodeHits<CFTFiberHit>("CFT");
  {
    const auto& U = HodoRawHit::kUp;
    Int_t nh=hodoAna.GetNHits("CFT");
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit<CFTFiberHit>("CFT", i);
      if(!hit) continue;
      Int_t seg   = hit->SegmentId();
      Int_t plane = hit->PlaneId();


      Double_t adccorHi  = hit->GetAdcCorHigh();
      Double_t adccorLow = hit->GetAdcCorLow();
      Double_t mipHi  = hit->GetMipHigh();
      Double_t mipLow = hit->GetMipLow();
      Double_t deLow  = hit->DeltaELowGain();

      HF2 (1000*(plane+1)+204, seg, adccorHi);
      HF2 (1000*(plane+1)+205, seg, adccorLow);
      HF2 (1000*(plane+1)+207, seg, mipHi);
      HF2 (1000*(plane+1)+208, seg, mipLow);
      HF2 (1000*(plane+1)+209, seg, deLow);
      event.cftadc_cor_h[plane][seg]  = (int)adccorHi;
      event.cftadc_cor_l[plane][seg] = (int)adccorLow;

      if (adccorHi>1500) {
	HF2 (1000*(plane+1)+201, seg, event.cfttdc[plane][seg]);
      }

      Int_t NhitT = hit->GetEntries(U);
      bool flagTime = false;
      for(Int_t m = 0; m<NhitT; ++m){
	Double_t time = hit->GetTUp(m);
	HF2 (1000*(plane+1)+200, seg, time);
	if (std::abs(time)<100)
	  flagTime = true;
      }
      if (flagTime)
	HF2 (1000*(plane+1)+206, seg, adccorHi);
    }

    Int_t nc=hodoAna.GetNClusters("CFT");

    Int_t ncluster[NumOfPlaneCFT];
    for (Int_t i=0; i<NumOfPlaneCFT; ++i)
      ncluster[i]=0;

    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster<CFTFiberCluster>("CFT", i);
      if(!cl) continue;
      Int_t plane = cl->PlaneId();
      Int_t cs=cl->ClusterSize();
      Double_t ms = cl->MeanSeg();
      Double_t cmt= cl->CMeanTime();
      Double_t total_de = cl->TotalDeltaE();
      Double_t max_de = cl->MaxDeltaE();

      ncluster[plane]++;
      HF1 (1000*(plane+1)+301, cs);
      HF2 (1000*(plane+1)+302, ms, cmt);
      HF2 (1000*(plane+1)+303, ms, total_de);
      HF2 (1000*(plane+1)+304, ms, max_de);

      HF1 (100*(plane+1) + 0, cmt);
      HF1 (100*(plane+1) + 1, ms);

      if ((cl->PlaneName()).Contains("PHI")) {
	//std::cout << "(x, y) = (" << cl->MeanX() << ", " << cl->MeanY() << ")" << std::endl;
	double meanPhi = cl->MeanPhi();
	HF1(1000*(plane+1)+305, meanPhi);
      } else if ((cl->PlaneName()).Contains("UV")) {
	//std::cout << "(r, z0) = (" << cl->MeanR() << ", " << cl->MeanZ0() << ")" << std::endl;
	double meanZ0 = cl->MeanZ0();
	HF1(1000*(plane+1)+305, meanZ0);
      }

      //std::cout << "Time : " << cmt << ", ";
      for (Int_t m=0; m<cs; ++m) {
	CFTFiberHit* hit = (CFTFiberHit*)cl->GetHit(m);
	Int_t planeId = hit->PlaneId();
	Int_t seg = hit->SegmentId();
	Double_t adccorHi  = hit->GetAdcCorHigh();
	//std::cout << seg << " (" << planeId << ", " << adccorHi << "), ";
	//if (m==cs-1)
	//std::cout << std::endl;
      }

    }
    for (Int_t i=0; i<NumOfPlaneCFT; ++i)
      HF1 (1000*(i+1)+300, ncluster[i]);
  }

  hodoAna.TimeCut("CFT", MinTimeCFT, MaxTimeCFT);
  hodoAna.AdcCut("CFT",  MinAdcCFT,  MaxAdcCFT);
  const auto& CFTClCont = hodoAna.GetClusterContainer("CFT");
  DCAna.DecodeCFTHits(CFTClCont);
  DCAna.TrackSearchCFT();

  Int_t ntCFT=DCAna.GetNtracksCFT();
  HF1(10, ntCFT);
  event.ntCFT = ntCFT;

  //for costheta
  Double_t u0_cos = 0;
  Double_t u1_cos = 0;
  Double_t v0_cos = 0;
  Double_t v1_cos = 0;

  for( Int_t i=0; i<ntCFT; ++i ){
    const auto& tp=DCAna.GetTrackCFT(i);

    Int_t nh   = tp->GetNHit();
    Int_t nhUV = tp->GetNHitUV();
    Double_t chisqrXY=tp->GetChiSquareXY();
    Double_t chisqrZ=tp->GetChiSquareZ();
    Double_t theta_cft =tp->GetThetaCFT();
    Double_t phi_cft   =tp->GetPhiCFT();
    Double_t total_dE  =tp->TotalDELowGain();
    Double_t max_dE    =tp->MaxDELowGain();
    Double_t norm_total_dE  =tp->NormalizedTotalDELowGain();
    Double_t norm_max_dE    =tp->NormalizedMaxDELowGain();
    Double_t total_dE_hi  =tp->TotalDEHiGain();

    Double_t x0 = tp->GetX0(), u0 = tp->GetU0();
    Double_t y0 = tp->GetY0(), v0 = tp->GetV0();

    HF1(11, nh+nhUV);
    HF1(12, nh);
    HF1(13, nhUV);
    HF1(14, chisqrXY);
    HF1(15, chisqrZ);
    HF1(17, theta_cft);
    HF1(18, phi_cft);
    HF2(19, phi_cft, theta_cft);
    HF1(20, total_dE);
    HF1(21, max_dE);
    HF1(22, norm_total_dE);
    HF1(23, norm_max_dE);
    HF1(24, total_dE_hi);
    HF2(25, max_dE, theta_cft);
    HF2(26, norm_max_dE, theta_cft);

    event.x0[0]=x0;
    event.y0[0]=y0;
    event.u0[0]=u0;
    event.v0[0]=v0;
    event.dist[0]=fabs(v0*x0-u0*y0)/sqrt(u0*u0+v0*v0);
    event.distMeanx[0]=v0*(v0*x0-u0*y0)/(2*(u0*u0+v0*v0));
    event.distMeany[0]=u0*(u0*y0-v0*x0)/(2*(u0*u0+v0*v0));
    event.distMeanz[0]=-(u0*x0+v0*y0)/(u0*u0+v0*v0);

    if(i==0){
      u0_cos = u0;
      v0_cos = v0;
    }
    if(i==1){
      u1_cos = u0;
      v1_cos = v0;
    }

    for (Int_t j=0; j<nh; ++j) {
      const auto& cl = tp->GetHit(j);
      Int_t plane = cl->PlaneId();
      HF1(16, plane);

      Double_t cmt= cl->CMeanTime();
      HF1 (100*(plane+1)+10, cmt);
      Double_t ms = cl->MeanSeg();
      HF1 (100*(plane+1)+11, ms);

      Double_t res = cl->GetResidual();
      Double_t res_phi = cl->GetResidualPhi();
      Double_t res_phi_cor = cl->GetResidualPhiCor();
      HF1 (100*(plane+1)+12, res);
      HF1 (100*(plane+1)+13, res_phi);
      HF1 (100*(plane+1)+14, res_phi_cor);

      Double_t z_cal   = cl->GetZcal();
      Double_t phi_cal = cl->GetCalPhi();
      HF2 (100*(plane+1)+15, z_cal, res_phi_cor);
      HF2 (100*(plane+1)+16, phi_cal, res_phi_cor);
    }

    for (Int_t j=0; j<nhUV; ++j) {
      const auto& cl = tp->GetHitUV(j);
      Int_t plane = cl->PlaneId();
      HF1(16, plane);

      Double_t cmt= cl->CMeanTime();
      HF1 (100*(plane+1)+10, cmt);
      Double_t ms = cl->MeanSeg();
      HF1 (100*(plane+1)+11, ms);

      Double_t res_z = cl->GetResidualZ();
      HF1 (100*(plane+1)+12, res_z);
      Double_t z_cal   = cl->GetZcal();
      Double_t phi_cal = cl->GetCalPhi();
      HF2 (100*(plane+1)+15, z_cal, res_z);
      HF2 (100*(plane+1)+16, phi_cal, res_z);
    }

    //std::cout << "nhPhi = " << nh << ", nhUV = " << nhUV
    //<< ", chisqrXY = " << chisqrXY << ", chisqrZ = "  << chisqrZ
    //<< ", theta = " << theta << std::endl;
  }

  // For Cosmic ray data
  if(ntCFT == 2){
    // 16 layer tracking
    DCAna.TrackSearchCFT_16layer();
    Int_t ntCFT_16layer = DCAna.GetNtracksCFT_16layer();
    event.ntCFT_16layer = ntCFT_16layer;

    HF1(50, ntCFT_16layer);
    if (ntCFT_16layer == 1) {
      const auto& tp = DCAna.GetTrackCFT_16layer(0);

      Int_t nhXY = tp->GetNHit();
      Int_t nhUV = tp->GetNHitU();
      Double_t chisqrXY      = tp->GetChiSquareXY();
      Double_t chisqrZ       = tp->GetChiSquareZ();
      Double_t theta_cft     = tp->GetThetaCFT();
      Double_t phi_cft       = tp->GetPhiCFT();
      Double_t total_dE      = tp->TotalDELowGain();
      Double_t max_dE        = tp->MaxDELowGain();
      Double_t norm_total_dE = tp->NormalizedTotalDELowGain();
      Double_t norm_max_dE   = tp->NormalizedMaxDELowGain();
      Double_t total_dE_hi  = tp->TotalDEHiGain();

      Double_t x0 = tp->GetX0(), u0 = tp->GetU0();
      Double_t y0 = tp->GetY0(), v0 = tp->GetV0();

      HF1(51, nhXY+nhUV);
      HF1(52, nhXY);
      HF1(53, nhUV);
      HF1(54, chisqrXY);
      HF1(55, chisqrZ);
      HF1(57, theta_cft);
      HF1(58, phi_cft);
      HF2(59, phi_cft, theta_cft);
      HF1(60, total_dE);
      HF1(61, max_dE);
      HF1(62, norm_total_dE);
      HF1(63, norm_max_dE);
      HF1(64, total_dE_hi);
      HF2(65, max_dE, theta_cft);
      HF2(66, norm_max_dE, theta_cft);


      event.nhXY_16layer[0]=nhXY;
      event.nhZ_16layer[0]=nhUV;
      event.chisqrXY_16layer[0]=chisqrXY;
      event.chisqrZ_16layer[0]=chisqrZ;
      event.x0_16layer[0]=x0;
      event.y0_16layer[0]=y0;
      event.u0_16layer[0]=u0;
      event.v0_16layer[0]=v0;
      event.z0_y_16layer[0]=-y0/v0;
      event.z0_x_16layer[0]=-x0/u0;
      event.dist_16layer[0]=fabs(v0*x0-u0*y0)/sqrt(u0*u0+v0*v0);
      event.costheta[0]= (u0_cos*u1_cos+v0_cos*v1_cos+1)/sqrt((u0_cos*u0_cos+v0_cos*v0_cos+1)*(u1_cos*u1_cos+v1_cos*v1_cos+1));

      for (Int_t ih=0; ih<nhXY; ih++) {
        CFTFiberCluster *cl = tp->GetHit(ih);
	Int_t  plane = cl->PlaneId();
	Double_t ms = cl->MeanSeg();
	Double_t cmt= cl->CMeanTime();
	Double_t res = cl->GetResidual();
	Double_t res_phi = cl->GetResidualPhi();
	Double_t res_phi_cor = cl->GetResidualPhiCor();
        Double_t phi = cl->MeanPhi();
        Double_t phi_cor = cl->MeanPhiCor();
        Double_t phi_cal1 = cl->GetCalPhi();
	Double_t z_cal   = cl->GetZcal();
	Double_t phi_cal = cl->GetCalPhi();
        if (0) {
	  std::cout << "layer : " << plane << ", MeanSeg : " << ms
                    << ", phi : " << phi << ", Z : " << z_cal
                    << ", ResidualPhi : " << res_phi_cor << std::endl;
        }

	HF1(56, plane);
	HF1 (100*(plane+1)+50, cmt);
	HF1 (100*(plane+1)+51, ms);
	HF1 (100*(plane+1)+52, res);
	HF1 (100*(plane+1)+53, res_phi);
	HF1 (100*(plane+1)+54, res_phi_cor);

	HF2 (100*(plane+1)+55, z_cal, res_phi_cor);
	HF2 (100*(plane+1)+56, phi_cal, res_phi_cor);

	/*
        int layer2 = (layer-layerId_PHI1)/2;
        int hid = 20000+layer2*1000;
        if (nhXY >= 7) {
          HF1(hid, res_phi);
          HF2(hid+1, calZ, res_phi);
          HF2(hid+1+int(phi/10.)+1, calZ, res_phi);
        }
	*/
        Int_t layer2 = (Int_t)(plane/2);
        event.hitMeanSegPhi_16layer[layer2][event.nhPhi_16layer[layer2]]  = ms;
        event.resPhi_16layer[layer2][event.nhPhi_16layer[layer2]]  = res_phi_cor;
        event.calZPhi_16layer[layer2][event.nhPhi_16layer[layer2]]  = z_cal;
        event.phi_16layer[layer2][event.nhPhi_16layer[layer2]]  = phi;
        event.phicor_16layer[layer2][event.nhPhi_16layer[layer2]]  = phi_cor;
        event.phical_16layer[layer2][event.nhPhi_16layer[layer2]]  = phi_cal1;
        event.nhPhi_16layer[layer2]++;

        /*
	     std::cout << "layer : " << layer2 << ", phi : " <<  phi
	     << ", int(phi) : " << int(phi)
	     << ", (int)phi" << (int)phi << std::endl;
	*/
      }

      for (Int_t ih=0; ih<nhUV; ih++) {
        CFTFiberCluster *cl = tp->GetHitU(ih);
	Int_t  plane     = cl->PlaneId();
	Double_t ms = cl->MeanSeg();
	Double_t cmt= cl->CMeanTime();
        Double_t phi_cal = cl->GetCalPhi();
        Double_t z_cal   = cl->GetZcal();
        Double_t res_z   = cl->GetResidualZ();

        if (0) {
	  std::cout << "layer : " << plane << ", MeanSeg : " << ms
                    << ", phi : " << phi_cal << ", Z : " << z_cal
                    << ", ResidualZ : " << res_z << std::endl;
        }

	HF1(56, plane);

	HF1 (100*(plane+1)+50, cmt);
	HF1 (100*(plane+1)+51, ms);
	HF1 (100*(plane+1)+52, res_z);
	HF2 (100*(plane+1)+55, z_cal, res_z);
	HF2 (100*(plane+1)+56, phi_cal, res_z);

        Int_t layer2 = (Int_t)(plane/2);
        event.hitMeanSegU_16layer[layer2][event.nhU_16layer[layer2]]  = ms;
        event.resZU_16layer[layer2][event.nhU_16layer[layer2]]  = res_z;
        event.calZU_16layer[layer2][event.nhU_16layer[layer2]]  = z_cal;
        event.phiU_16layer[layer2][event.nhU_16layer[layer2]]  = phi_cal;
        event.nhU_16layer[layer2]++;

        /*
	   std::cout << "layer : " << layer2 << ", phi : " <<  phi
	   << ", int(phi) : " << int(phi)
	   << ", (int)phi" << (int)phi << std::endl;
	*/
      }

      ThreeVector Pos0 = tp->GetPos0();
      ThreeVector Dir = tp->GetDir();
      // double scale = 10.;
      ThreeVector bMom = ThreeVector(0, 0, 1);
      double cost = Dir*bMom/(Dir.Mag()*bMom.Mag());
      double theta = acos(cost)*TMath::RadToDeg();

      event.theta_16layer[0]=theta;
      event.cost_16layer[0]=cost;
    }
  }



  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingEnd()
{
  tree->Fill();
  // hodo->Fill();
  return true;
}

//_____________________________________________________________________________
namespace
{
const Int_t    NbinAdc = 4096;
const Double_t MinAdc  =    0.;
const Double_t MaxAdc  = 4096.;

const Int_t    NbinTdc = 4096;
const Double_t MinTdc  =    0.;
const Double_t MaxTdc  = 4096.;

const Int_t    NbinTdcHr = 1e6/10;
const Double_t MinTdcHr  =  0.;
const Double_t MaxTdcHr  = 1e6;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 20, 0., 20.);

  HB1(10, "#Tracks CFT", 10, 0., 10.);
  HB1(11, "#Hits of Track CFT", 10, 0., 10.);
  HB1(12, "#Hits (PHI) of Track CFT", 5, 0., 5.);
  HB1(13, "#Hits (UV) of Track CFT", 5, 0., 5.);
  HB1(14, "Chisqr(PHI)", 200, 0., 100.);
  HB1(15, "Chisqr(UV)", 200, 0., 100.);
  HB1(16, "LayerId CFT", 10, 0., 10.);
  HB1(17, "Theta of CFT Track", 200, -10., 190.);
  HB1(18, "Phi of CFT Track",   200, -200., 400.);
  HB2(19, "Theta % Phi of CFT Track",   200, -200., 400., 200, -10, 190);
  HB1(20, "Total dE (Low gain) of CFT Track",   200, -2., 18.);
  HB1(21, "Max dE (Low gain) of CFT Track",   200, -2., 18.);
  HB1(22, "Normalized Total dE (Low gain) of CFT Track",   200, 0.0, 2.5);
  HB1(23, "Normalized Max dE (Low gain) of CFT Track",   200, 0.0, 2.5);
  HB1(24, "Total dE (High gain) of CFT Track",   200, 0., 2000.);
  HB2(25, "Theta % Max dE (Low gain) of CFT Track",   200, -2., 18., 200, -10, 190);
  HB2(26, "Theta % Normalized Max dE (Low gain) of CFT Track",   200, 0.0, 2.5, 200, -10, 190);

  HB1(50, "#Tracks CFT", 10, 0., 10.);
  HB1(51, "#Hits of Track CFT", 20, 0., 20.);
  HB1(52, "#Hits (PHI) of Track CFT", 10, 0., 10.);
  HB1(53, "#Hits (UV) of Track CFT", 10, 0., 10.);
  HB1(54, "Chisqr(PHI)", 200, 0., 100.);
  HB1(55, "Chisqr(UV)", 200, 0., 100.);
  HB1(56, "LayerId CFT", 10, 0., 10.);
  HB1(57, "Theta of CFT Track", 200, -10., 190.);
  HB1(58, "Phi of CFT Track",   200, -200., 400.);
  HB2(59, "Theta % Phi of CFT Track",   200, -200., 400., 200, -10, 190);
  HB1(60, "Total dE (Low gain) of CFT Track",   200, -2., 18.);
  HB1(61, "Max dE (Low gain) of CFT Track",   200, -2., 18.);
  HB1(62, "Normalized Total dE (Low gain) of CFT Track",   200, 0.0, 2.5);
  HB1(63, "Normalized Max dE (Low gain) of CFT Track",   200, 0.0, 2.5);
  HB1(64, "Total dE (High gain) of CFT Track",   200, 0., 2000.);
  HB2(65, "Theta % Max dE (Low gain) of CFT Track",   200, -2., 18., 200, -10, 190);
  HB2(66, "Theta % Normalized Max dE (Low gain) of CFT Track",   200, 0.0, 2.5, 200, -10, 190);

  for(Int_t i=0; i<NumOfPlaneCFT; i++){
    TString title100("");
    TString title101("");

    TString title110("");
    TString title111("");
    TString title112("");
    TString title113("");
    TString title114("");
    TString title115("");
    TString title116("");

    TString title150("");
    TString title151("");
    TString title152("");
    TString title153("");
    TString title154("");
    TString title155("");
    TString title156("");
    if(i%2 == 0){// spiral layer
      Int_t layer = (Int_t)i/2 +1;
      title100  = Form("CFT UV %d : Time", layer);
      title101  = Form("CFT UV %d : MeanSeg", layer);

      title110  = Form("CFT UV %d : Time [Track]", layer);
      title111  = Form("CFT UV %d : MeanSeg [Track]", layer);
      title112  = Form("CFT UV %d : Residual (Z) [Track]", layer);
      title113  = Form("CFT UV %d : Not Used", layer);
      title114  = Form("CFT UV %d : Not Used", layer);
      title115  = Form("CFT UV %d : Residual (Z) % Z [Track]", layer);
      title116  = Form("CFT UV %d : Residual (Z) % Phi [Track]", layer);

      title150  = Form("CFT UV %d : Time [16layer Track]", layer);
      title151  = Form("CFT UV %d : MeanSeg [16layer Track]", layer);
      title152  = Form("CFT UV %d : Residual (Z) [16layer Track]", layer);
      title153  = Form("CFT UV %d : Not Used", layer);
      title154  = Form("CFT UV %d : Not Used", layer);
      title155  = Form("CFT UV %d : Residual (Z) % Z [16layer Track]", layer);
      title156  = Form("CFT UV %d : Residual (Z) % Phi [16layer Track]", layer);
    }else if(i%2 == 1){// straight layer
      Int_t layer = (Int_t)i/2 +1;
      title100  = Form("CFT Phi %d : Time", layer);
      title101  = Form("CFT Phi %d : MeanSeg", layer);

      title110  = Form("CFT Phi %d : Time [Track]", layer);
      title111  = Form("CFT Phi %d : MeanSeg [Track]", layer);
      title112  = Form("CFT Phi %d : Residual (X or Y) [Track]", layer);
      title113  = Form("CFT Phi %d : Residual (phi) [Track]", layer);
      title114  = Form("CFT Phi %d : Residual (phi_cor) [Track]", layer);
      title115  = Form("CFT Phi %d : Residual (phi_cor) % Z [Track]", layer);
      title116  = Form("CFT Phi %d : Residual (phi_cor) % Phi [Track]", layer);

      title150  = Form("CFT Phi %d : Time [16layer Track]", layer);
      title151  = Form("CFT Phi %d : MeanSeg [16layer Track]", layer);
      title152  = Form("CFT Phi %d : Residual (X or Y) [16layer Track]", layer);
      title153  = Form("CFT Phi %d : Residual (phi) [16layer Track]", layer);
      title154  = Form("CFT Phi %d : Residual (phi_cor) [16layer Track]", layer);
      title155  = Form("CFT Phi %d : Residual (phi_cor) % Z [16layer Track]", layer);
      title156  = Form("CFT Phi %d : Residual (phi_cor) % Phi [16layer Track]", layer);
    }

    HB1( 100*(i+1) +  0, title100, 1000,-50, 150);
    HB1( 100*(i+1) +  1, title101, NumOfSegCFT[i], 0, NumOfSegCFT[i]);

    HB1( 100*(i+1) + 10, title110, 1000,-50, 150);
    HB1( 100*(i+1) + 11, title111, NumOfSegCFT[i], 0, NumOfSegCFT[i]);
    HB1( 100*(i+1) + 12, title112, 100, -5., 5.);
    HB1( 100*(i+1) + 13, title113, 100, -5., 5.);
    HB1( 100*(i+1) + 14, title114, 100, -5., 5.);
    HB2( 100*(i+1) + 15, title115, 500, -200, 300, 100, -5., 5.);
    HB2( 100*(i+1) + 16, title116, 400, -20, 380, 100, -5., 5.);

    HB1( 100*(i+1) + 50, title110, 1000,-50, 150);
    HB1( 100*(i+1) + 51, title111, NumOfSegCFT[i], 0, NumOfSegCFT[i]);
    HB1( 100*(i+1) + 52, title112, 100, -5., 5.);
    HB1( 100*(i+1) + 53, title113, 100, -5., 5.);
    HB1( 100*(i+1) + 54, title114, 100, -5., 5.);
    HB2( 100*(i+1) + 55, title115, 500, -200, 300, 100, -5., 5.);
    HB2( 100*(i+1) + 56, title116, 400, -20, 380, 100, -5., 5.);
  }



  for(Int_t i=0; i<NumOfPlaneCFT; i++){
    //ADC
    TString title100("");
    TString title101("");
    TString title102("");
    TString title103("");
    TString title104("");
    TString title105("");
    TString title106("");
    TString title200("");
    TString title201("");
    TString title204("");
    TString title205("");
    TString title206("");
    TString title207("");
    TString title208("");
    TString title209("");
    TString title300("");
    TString title301("");
    TString title302("");
    TString title303("");
    TString title304("");
    TString title305("");
    if(i%2 == 0){// spiral layer
      Int_t layer = (Int_t)i/2 +1;
      title100  = Form("CFT UV %d : Tdc(Leading) vs seg", layer);
      title101  = Form("CFT UV %d : Tdc(Trailing) vs seg", layer);
      title102  = Form("CFT UV %d : Hit pattern", layer);
      title103  = Form("CFT UV %d : TOT vs seg", layer);
      title104  = Form("CFT UV %d : Adc(High) vs seg", layer);
      title105  = Form("CFT UV %d : Adc(Low) vs seg", layer);
      title106  = Form("CFT UV %d : Adc(High) vs seg (w/ TDC)", layer);
      title200 = Form("CFT UV %d : Time vs seg", layer);
      title201 = Form("CFT UV %d : TDC (w/ Large ADC) vs seg", layer);
      title204 = Form("CFT UV %d : AdcCor(High) vs seg", layer);
      title205 = Form("CFT UV %d : AdcCor(Low) vs seg", layer);
      title206 = Form("CFT UV %d : AdcCor(High) vs seg (w/ TDC)", layer);
      title207 = Form("CFT UV %d : Mip Calib(High) vs seg", layer);
      title208 = Form("CFT UV %d : Mip Calib(Low) vs seg", layer);
      title209 = Form("CFT UV %d : dE (Low) vs seg", layer);
      title300 = Form("CFT UV %d : nCluster", layer);
      title301 = Form("CFT UV %d : Cluster Size", layer);
      title302 = Form("CFT UV %d (Cluster): Time vs seg", layer);
      title303 = Form("CFT UV %d (Cluster): total dE vs seg", layer);
      title304 = Form("CFT UV %d (Cluster): max dE vs seg", layer);
      title305 = Form("CFT UV %d (Cluster): MeanZ0", layer);
    }else if(i%2 == 1){// straight layer
      Int_t layer = (Int_t)i/2 +1;
      title100  = Form("CFT Phi %d : Tdc(Leading) vs seg", layer);
      title101  = Form("CFT Phi %d : Tdc(Trailing) vs seg", layer);
      title102  = Form("CFT Phi %d : Hit pattern", layer);
      title103  = Form("CFT Phi %d : TOT vs seg", layer);
      title104  = Form("CFT Phi %d : Adc(High) vs seg", layer);
      title105  = Form("CFT Phi %d : Adc(Low) vs seg", layer);
      title106  = Form("CFT Phi %d : Adc(High) vs seg (w/ TDC)", layer);
      title200 = Form("CFT Phi %d : Time vs seg", layer);
      title201 = Form("CFT Phi %d : TDC (w/ Large ADC) vs seg", layer);
      title204 = Form("CFT Phi %d : AdcCor(High) vs seg", layer);
      title205 = Form("CFT Phi %d : AdcCor(Low) vs seg", layer);
      title206 = Form("CFT Phi %d : AdcCor(High) vs seg (w/ TDC)", layer);
      title207 = Form("CFT Phi %d : Mip Calib(High) vs seg", layer);
      title208 = Form("CFT Phi %d : Mip Calib(Low) vs seg", layer);
      title209 = Form("CFT Phi %d : dE (Low) vs seg", layer);
      title300 = Form("CFT Phi %d : nCluster", layer);
      title301 = Form("CFT Phi %d : Cluster Size", layer);
      title302 = Form("CFT Phi %d (Cluster): Time vs seg", layer);
      title303 = Form("CFT Phi %d (Cluster): total dE vs seg", layer);
      title304 = Form("CFT Phi %d (Cluster): max dE vs seg", layer);
      title305 = Form("CFT Phi %d (Cluster): MeanPhi", layer);
    }
    HB2( 1000*(i+1)+100, title100, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1024,0,1024);
    HB2( 1000*(i+1)+101, title101, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1024,0,1024);
    HB1( 1000*(i+1)+102, title102, NumOfSegCFT[i], 0, NumOfSegCFT[i]);
    HB2( 1000*(i+1)+103, title103, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1024,0,1024);
    HB2( 1000*(i+1)+104, title104, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,0,4000);
    HB2( 1000*(i+1)+105, title105, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,0,4000);
    HB2( 1000*(i+1)+106, title106, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,0,4000);
    HB2( 1000*(i+1)+200, title200, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-500,500);
    HB2( 1000*(i+1)+201, title201, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1024,0,1024);
    HB2( 1000*(i+1)+204, title204, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-500,3500);
    HB2( 1000*(i+1)+205, title205, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-500,3500);
    HB2( 1000*(i+1)+206, title206, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-500,3500);
    HB2( 1000*(i+1)+207, title207, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-10,90);
    HB2( 1000*(i+1)+208, title208, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-1,9);
    HB2( 1000*(i+1)+209, title209, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-1,9);
    HB1( 1000*(i+1)+300, title300, 20, 0, 20);
    HB1( 1000*(i+1)+301, title301, 20, 0, 20);
    HB2( 1000*(i+1)+302, title302, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-500,500);
    HB2( 1000*(i+1)+303, title303, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-1,9);
    HB2( 1000*(i+1)+304, title304, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-1,9);
    HB1( 1000*(i+1)+305, title305, 500, -50, 450);
  }



  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spill",     &event.spill,     "spill/I");
  //Trig
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  tree->Branch("cftadc_h",   event.cftadc_h,  Form("cftadc_h[%d][%d]/I", NumOfPlaneCFT, NumOfSegCFT_PHI4));
  tree->Branch("cftadc_l",   event.cftadc_l,  Form("cftadc_l[%d][%d]/I", NumOfPlaneCFT, NumOfSegCFT_PHI4));
  tree->Branch("cfttdc",     event.cfttdc,    Form("cfttdc[%d][%d]/I", NumOfPlaneCFT, NumOfSegCFT_PHI4));

  tree->Branch("cftadc_cor_h",   event.cftadc_cor_h,  Form("cftadc_cor_h[%d][%d]/I", NumOfPlaneCFT, NumOfSegCFT_PHI4));
  tree->Branch("cftadc_cor_l",   event.cftadc_cor_l,  Form("cftadc_cor_l[%d][%d]/I", NumOfPlaneCFT, NumOfSegCFT_PHI4));


  // 16layer tracking
  tree->Branch("ntCFT",   &event.ntCFT,  "ntCFT/I");
  tree->Branch("x0",   event.x0,  "x0[ntCFT]/D");
  tree->Branch("y0",   event.y0,  "y0[ntCFT]/D");
  tree->Branch("u0",   event.u0,  "u0[ntCFT]/D");
  tree->Branch("v0",   event.v0,  "v0[ntCFT]/D");
  tree->Branch("dist",   event.dist,  "dist[ntCFT]/D");
  tree->Branch("distMeanx",   event.distMeanx,  "distMeanx[ntCFT]/D");
  tree->Branch("distMeany",   event.distMeany,  "distMeany[ntCFT]/D");
  tree->Branch("distMeanz",   event.distMeanz,  "distMeanz[ntCFT]/D");
  tree->Branch("cos_16layer", event.cos_16layer, "cos_16layer[ntCFT]/D");
  tree->Branch("ntCFT_16layer",   &event.ntCFT_16layer,  "ntCFT_16layer/I");
  tree->Branch("nhXY_16layer",   event.nhXY_16layer,  "nhXY_16layer[ntCFT_16layer]/I");
  tree->Branch("nhZ_16layer",   event.nhZ_16layer,  "nhZ_16layer[ntCFT_16layer]/I");
  tree->Branch("chisqrXY_16layer",   event.chisqrXY_16layer,  "chisqrXY_16layer[ntCFT_16layer]/D");
  tree->Branch("chisqrZ_16layer",   event.chisqrZ_16layer,  "chisqrZ_16layer[ntCFT_16layer]/D");
  tree->Branch("x0_16layer",   event.x0_16layer,  "x0_16layer[ntCFT_16layer]/D");
  tree->Branch("y0_16layer",   event.y0_16layer,  "y0_16layer[ntCFT_16layer]/D");
  tree->Branch("u0_16layer",   event.u0_16layer,  "u0_16layer[ntCFT_16layer]/D");
  tree->Branch("v0_16layer",   event.v0_16layer,  "v0_16layer[ntCFT_16layer]/D");
  tree->Branch("dist_16layer",   event.dist_16layer,  "dist_16layer[ntCFT_16layer]/D");
  tree->Branch("costheta",   event.costheta,  "costheta[ntCFT_16layer]/D");
  tree->Branch("theta_16layer",   event.theta_16layer,  "theta_16layer[ntCFT_16layer]/D");
  tree->Branch("cost_16layer",   event.cost_16layer,  "cost_16layer[ntCFT_16layer]/D");
  tree->Branch("nhPhi1_16layer", &event.nhPhi_16layer[0], "nhPhi1_16layer/I");
  tree->Branch("nhPhi2_16layer", &event.nhPhi_16layer[1], "nhPhi2_16layer/I");
  tree->Branch("nhPhi3_16layer", &event.nhPhi_16layer[2], "nhPhi3_16layer/I");
  tree->Branch("nhPhi4_16layer", &event.nhPhi_16layer[3], "nhPhi4_16layer/I");

  tree->Branch("hitMeanSegPhi1_16layer",   event.hitMeanSegPhi_16layer[0],  "hitMeanSegPhi1_16layer[nhPhi1_16layer]/D");
  tree->Branch("hitMeanSegPhi2_16layer",   event.hitMeanSegPhi_16layer[1],  "hitMeanSegPhi2_16layer[nhPhi2_16layer]/D");
  tree->Branch("hitMeanSegPhi3_16layer",   event.hitMeanSegPhi_16layer[2],  "hitMeanSegPhi3_16layer[nhPhi3_16layer]/D");
  tree->Branch("hitMeanSegPhi4_16layer",   event.hitMeanSegPhi_16layer[3],  "hitMeanSegPhi4_16layer[nhPhi4_16layer]/D");

  tree->Branch("phi1_16layer",   event.phi_16layer[0],  "phi1_16layer[nhPhi1_16layer]/D");
  tree->Branch("phi2_16layer",   event.phi_16layer[1],  "phi2_16layer[nhPhi2_16layer]/D");
  tree->Branch("phi3_16layer",   event.phi_16layer[2],  "phi3_16layer[nhPhi3_16layer]/D");
  tree->Branch("phi4_16layer",   event.phi_16layer[3],  "phi4_16layer[nhPhi4_16layer]/D");

  tree->Branch("phi1cor_16layer",   event.phicor_16layer[0],  "phi1cor_16layer[nhPhi1_16layer]/D");
  tree->Branch("phi2cor_16layer",   event.phicor_16layer[1],  "phi2cor_16layer[nhPhi2_16layer]/D");
  tree->Branch("phi3cor_16layer",   event.phicor_16layer[2],  "phi3cor_16layer[nhPhi3_16layer]/D");
  tree->Branch("phi4cor_16layer",   event.phicor_16layer[3],  "phi4cor_16layer[nhPhi4_16layer]/D");

    tree->Branch("phi1cal_16layer",   event.phical_16layer[0],  "phi1cal_16layer[nhPhi1_16layer]/D");
  tree->Branch("phi2cal_16layer",   event.phical_16layer[1],  "phi2cal_16layer[nhPhi2_16layer]/D");
  tree->Branch("phi3cal_16layer",   event.phical_16layer[2],  "phi3cal_16layer[nhPhi3_16layer]/D");
  tree->Branch("phi4cal_16layer",   event.phical_16layer[3],  "phi4cal_16layer[nhPhi4_16layer]/D");

  tree->Branch("resPhi1_16layer",   event.resPhi_16layer[0],  "resPhi1_16layer[nhPhi1_16layer]/D");
  tree->Branch("resPhi2_16layer",   event.resPhi_16layer[1],  "resPhi2_16layer[nhPhi2_16layer]/D");
  tree->Branch("resPhi3_16layer",   event.resPhi_16layer[2],  "resPhi3_16layer[nhPhi3_16layer]/D");
  tree->Branch("resPhi4_16layer",   event.resPhi_16layer[3],  "resPhi4_16layer[nhPhi4_16layer]/D");

  tree->Branch("calZPhi1_16layer",   event.calZPhi_16layer[0],  "calZPhi1_16layer[nhPhi1_16layer]/D");
  tree->Branch("calZPhi2_16layer",   event.calZPhi_16layer[1],  "calZPhi2_16layer[nhPhi2_16layer]/D");
  tree->Branch("calZPhi3_16layer",   event.calZPhi_16layer[2],  "calZPhi3_16layer[nhPhi3_16layer]/D");
  tree->Branch("calZPhi4_16layer",   event.calZPhi_16layer[3],  "calZPhi4_16layer[nhPhi4_16layer]/D");

  tree->Branch("nhU1_16layer", &event.nhU_16layer[0], "nhU1_16layer/I");
  tree->Branch("nhU2_16layer", &event.nhU_16layer[1], "nhU2_16layer/I");
  tree->Branch("nhU3_16layer", &event.nhU_16layer[2], "nhU3_16layer/I");
  tree->Branch("nhU4_16layer", &event.nhU_16layer[3], "nhU4_16layer/I");

  tree->Branch("hitMeanSegU1_16layer",   event.hitMeanSegU_16layer[0],  "hitMeanSegU1_16layer[nhU1_16layer]/D");
  tree->Branch("hitMeanSegU2_16layer",   event.hitMeanSegU_16layer[1],  "hitMeanSegU2_16layer[nhU2_16layer]/D");
  tree->Branch("hitMeanSegU3_16layer",   event.hitMeanSegU_16layer[2],  "hitMeanSegU3_16layer[nhU3_16layer]/D");
  tree->Branch("hitMeanSegU4_16layer",   event.hitMeanSegU_16layer[3],  "hitMeanSegU4_16layer[nhU4_16layer]/D");

  tree->Branch("phiU1_16layer",   event.phiU_16layer[0],  "phiU1_16layer[nhU1_16layer]/D");
  tree->Branch("phiU2_16layer",   event.phiU_16layer[1],  "phiU2_16layer[nhU2_16layer]/D");
  tree->Branch("phiU3_16layer",   event.phiU_16layer[2],  "phiU3_16layer[nhU3_16layer]/D");
  tree->Branch("phiU4_16layer",   event.phiU_16layer[3],  "phiU4_16layer[nhU4_16layer]/D");

  tree->Branch("resZU1_16layer",   event.resZU_16layer[0],  "resZU1_16layer[nhU1_16layer]/D");
  tree->Branch("resZU2_16layer",   event.resZU_16layer[1],  "resZU2_16layer[nhU2_16layer]/D");
  tree->Branch("resZU3_16layer",   event.resZU_16layer[2],  "resZU3_16layer[nhU3_16layer]/D");
  tree->Branch("resZU4_16layer",   event.resZU_16layer[3],  "resZU4_16layer[nhU4_16layer]/D");

  tree->Branch("calZU1_16layer",   event.calZU_16layer[0],  "calZU1_16layer[nhU1_16layer]/D");
  tree->Branch("calZU2_16layer",   event.calZU_16layer[1],  "calZU2_16layer[nhU2_16layer]/D");
  tree->Branch("calZU3_16layer",   event.calZU_16layer[2],  "calZU3_16layer[nhU3_16layer]/D");
  tree->Branch("calZU4_16layer",   event.calZU_16layer[3],  "calZU4_16layer[nhU4_16layer]/D");


  // HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
     InitializeParameter<HodoParamMan>("HDPRM") &&
     InitializeParameter<HodoPHCMan>("HDPHC")   &&
     InitializeParameter<UserParamMan>("USER")   &&
     InitializeParameter<CFTPedCorMan>("CFTPED") &&
     InitializeParameter<CFTPosParamMan>("CFTPOS")
     );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
