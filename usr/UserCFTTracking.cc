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
#include "Kinematics.hh"

// #define TimeCut    1 // in cluster analysis
#define FHitBranch 0 // make FiberHit branches (becomes heavy)
#define HodoHitPos 0
#define ForPHC 0

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

  Double_t track_theta[MaxDepth];
  Double_t track_phi[MaxDepth];

  Int_t    nVer_product;
  Double_t dist_product[MaxDepth];
  Double_t xtar_product[MaxDepth];
  Double_t ytar_product[MaxDepth];
  Double_t ztar_product[MaxDepth];

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
  nVer_product = 0;

  for (Int_t j=0; j<MaxDepth; j++) {
    x0[j] = qnan;
    y0[j] = qnan;
    u0[j] = qnan;
    v0[j] = qnan;
    dist[j] = qnan;
    distMeanx[j] = qnan;
    distMeany[j] = qnan;
    distMeanz[j] = qnan;
    track_theta[j] = qnan;
    track_phi[j] = qnan;
    dist_product[j] = qnan;
    xtar_product[j] = qnan;
    ytar_product[j] = qnan;
    ztar_product[j] = qnan;
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
#if HodoHitPos
  static const auto PropVelBH2 = gUser.GetParameter("PropagationBH2");
#endif

  RawData rawData;
  HodoAnalyzer hodoAna(rawData);
  DCAnalyzer   DCAna(rawData);

  gRM.Decode();

  event.evnum = gUnpacker.get_event_number();
  event.spill = gRM.SpillNumber();
  dst.evnum   = gRM.EventNumber();
  dst.spill   = gRM.SpillNumber();

  HF1(1, 0);

  // ForPHC
  Int_t maxplane;
  Int_t maxseg;
  Double_t maxadc;
  Double_t maxcmt;

  //***************** RawData  *****************

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

  // if(trigger_flag[kTriggerFlag::SpillOnEnd] || trigger_flag[kTriggerFlag::SpillOffEnd])
  //   return true;

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
        // event.cftadc_h[plane][seg] = adcH;

        if (flag_tdc) {
          HF2 (1000*(plane+1)+106, seg, adcH);
        }
      }

      //ADC Low
      for(Int_t m = 0; m<NhitAL; ++m){
        Int_t adcL = hit->GetAdcLow();
        HF2 (1000*(plane+1)+105, seg, adcL);
        // event.cftadc_l[plane][seg] = adcL;
      }
    }
  }

  hodoAna.DecodeHits<CFTFiberHit>("CFT");//ここ
  {
    const auto& U = HodoRawHit::kUp;
    Int_t nh=hodoAna.GetNHits("CFT");
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit<CFTFiberHit>("CFT", i);
      if(!hit) continue;
      Int_t seg   = hit->SegmentId();
      Int_t plane = hit->PlaneId();//レイヤー


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
        double meanPhi = cl->MeanPhi();
        HF1(1000*(plane+1)+305, meanPhi);
      } else if ((cl->PlaneName()).Contains("UV")) {
        double meanZ0 = cl->MeanZ0();
        HF1(1000*(plane+1)+305, meanZ0);
      }

      for (Int_t m=0; m<cs; ++m) {
        CFTFiberHit* hit = (CFTFiberHit*)cl->GetHit(m);
	Int_t planeId = hit->PlaneId();
	Int_t seg = hit->SegmentId();
	Double_t adccorHi  = hit->GetAdcCorHigh();
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


  if (ntCFT > 0) {

    const Int_t n_all_hits = hodoAna.GetNHits("CFT");
    for (Int_t i = 0; i < n_all_hits; ++i) {
      const auto* hit = hodoAna.GetHit<CFTFiberHit>("CFT", i);
      if (!hit) continue;

      const Int_t plane = hit->PlaneId();
      const Int_t seg = hit->SegmentId();
      const Double_t adccorHi = hit->GetAdcCorHigh();


      HF2(1000 * (plane + 1) + 307, seg, adccorHi);
    }
  }


  ThreeVector track_Dir[MaxDepth];
  ThreeVector track_Pos[MaxDepth];

  for( Int_t i=0; i<ntCFT; ++i ){
    const auto& tp = DCAna.GetTrackCFT(i);

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

    event.x0[i]=x0;
    event.y0[i]=y0;
    event.u0[i]=u0;
    event.v0[i]=v0;
    event.dist[i]=fabs(v0*x0-u0*y0)/sqrt(u0*u0+v0*v0);
    event.distMeanx[i]=v0*(v0*x0-u0*y0)/(2*(u0*u0+v0*v0));
    event.distMeany[i]=u0*(u0*y0-v0*x0)/(2*(u0*u0+v0*v0));
    event.distMeanz[i]=-(u0*x0+v0*y0)/(u0*u0+v0*v0);

    track_Dir[i] = tp->GetDir();
    track_Pos[i] = tp->GetPos0();
    event.track_theta[i] = track_Dir[i].Theta()*TMath::RadToDeg();
    event.track_phi[i] = track_Dir[i].Phi()*TMath::RadToDeg();

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
      HF2 (100*(plane+1)+17, z_cal, cmt);

      Int_t cs = cl->ClusterSize();
      for (Int_t k = 0; k < cs; ++k) {
        const auto* hit = dynamic_cast<const CFTFiberHit*>(cl->GetHit(k));
        if (!hit) continue;
        Int_t seg = hit->SegmentId();
        Double_t adccorHi = hit->GetAdcCorHigh();
        HF2(1000*(plane+1)+306, seg, adccorHi);
      }

#if ForPHC
      maxplane = plane;
      maxseg   = cl->MaxSegment();
      maxadc   = cl->MaxAdcHigh();
      maxcmt   = cmt;

      if(std::abs(maxcmt)<100 && maxseg>-1 && ntCFT>0){
	Int_t histId = ((maxplane+1)*1000+maxseg)*10;
	if (maxadc>1000){
	  HF1( histId+5, maxcmt);
	}
	HF2( histId+6, 1./sqrt(maxadc), maxcmt);
      }
#endif

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
      HF2 (100*(plane+1)+17, z_cal, cmt);

      Int_t cs = cl->ClusterSize();
      for (Int_t k = 0; k < cs; ++k) {
        const auto* hit = dynamic_cast<const CFTFiberHit*>(cl->GetHit(k));
        if (!hit) continue;
        Int_t seg = hit->SegmentId();
        Double_t adccorHi = hit->GetAdcCorHigh();
        HF2(1000*(plane+1)+306, seg, adccorHi);
      }

#if ForPHC
      maxplane = plane;
      maxseg   = cl->MaxSegment();
      maxadc   = cl->MaxAdcHigh();
      maxcmt   = cmt;

      if(std::abs(maxcmt)<100 && maxseg>-1 && ntCFT>0){
	Int_t histId = ((maxplane+1)*1000+maxseg)*10;
	if (maxadc>1000){
	  HF1( histId+5, maxcmt);
	}
	HF2( histId+6, 1./sqrt(maxadc), maxcmt);
      }
#endif

    }
  }

  //For production data
  if(ntCFT>1){
    Int_t nVer = 0;
    TVector3 target;
    Double_t dist_product;
    for(Int_t i=0; i<ntCFT-1; i++){
      for(Int_t j=i+1; j<ntCFT; j++){
        target = Kinematics::VertexPoint3D(track_Pos[i], track_Pos[j], track_Dir[i], track_Dir[j], dist_product);
        event.xtar_product[nVer] = target.X();
        event.ytar_product[nVer] = target.Y();
        event.ztar_product[nVer] = target.Z();
        event.dist_product[nVer] = dist_product;

        HF2 (30, target.X(), target.Y());
        HF2 (31, target.Z(), target.X());
        HF2 (32, target.Z(), target.Y());

        nVer++;
      }
    }
    event.nVer_product = nVer;
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingEnd()
{
#if ForPHC
#else
  tree->Fill();
  // hodo->Fill();
#endif
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

  HB2(30, "x vs y of CFT Track",   100, -50, 50, 100, -50, 50);
  HB2(31, "z vs x of CFT Track",   500, -200, 300, 100, -50, 50);
  HB2(32, "z vs y of CFT Track",   500, -200, 300, 100, -50, 50);

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
    TString title117("");

    if(i%2 == 0){// spiral layer
      Int_t layer = (Int_t)i/2 +1;
      title100  = Form("CFT UV %d : Time", layer);
      title101  = Form("CFT UV %d : MeanSeg", layer);

      title110  = Form("CFT UV %d : Time [Track]", layer);
      title111  = Form("CFT UV %d : MeanSeg [Track]", layer);
      title112  = Form("CFT UV %d : Residual (Z) [Track]", layer);
      title113  = Form("CFT UV %d : Not Used", layer);
      title114  = Form("CFT UV %d : Not Used", layer);
      title115  = Form("CFT UV %d : Residual (Z) vs Z [Track]", layer);
      title116  = Form("CFT UV %d : Residual (Z) vs Phi [Track]", layer);
      title117  = Form("CFT UV %d : Z vs Time [Track]", layer);

    }else if(i%2 == 1){// straight layer
      Int_t layer = (Int_t)i/2 +1;
      title100  = Form("CFT Phi %d : Time", layer);
      title101  = Form("CFT Phi %d : MeanSeg", layer);

      title110  = Form("CFT Phi %d : Time [Track]", layer);
      title111  = Form("CFT Phi %d : MeanSeg [Track]", layer);
      title112  = Form("CFT Phi %d : Residual (X or Y) [Track]", layer);
      title113  = Form("CFT Phi %d : Residual (phi) [Track]", layer);
      title114  = Form("CFT Phi %d : Residual (phi_cor) [Track]", layer);
      title115  = Form("CFT Phi %d : Residual (phi_cor) vs Z [Track]", layer);
      title116  = Form("CFT Phi %d : Residual (phi_cor) vs Phi [Track]", layer);
      title117  = Form("CFT Phi %d : Z vs Time [Track]", layer);

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
    HB2( 100*(i+1) + 17, title117, 1000, -200, 300, 700, -20., 50.);

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
    TString title306("");
    TString title307("");
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
      title306 = Form("CFT UV %d : AdcCor(High) vs seg (Track_only)", layer);
      title307 = Form("CFT UV %d : AdcCor(High) vs seg (Track)", layer);
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
      title306 = Form("CFT Phi %d : AdcCor(High) vs seg (Track_only)", layer);
      title307 = Form("CFT Phi %d : AdcCor(High) vs seg (Track)", layer);

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
    HB2( 1000*(i+1)+306, title306, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000, 0, 4000 );
    HB2( 1000*(i+1)+307, title307, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000, 0, 4000 );
  }

#if ForPHC
  for (int l=0; l<NumOfPlaneCFT; l++) {
    Int_t NumOfSeg = NumOfSegCFT[l];
    for (int seg=0; seg<NumOfSeg; seg++ ) {
      Int_t hid = ((l+1)*1000+seg)*10;
      TString title5 = Form("CTime (ADCcor>1000) %d-%d", l, seg);
      TString title6 = Form("CTime : 1/sqrt(MaxClADCcor) %d-%d", l, seg);
      HB1(hid+5, title5, 100, -50, 50);
      HB2(hid+6, title6, 100, 0, 0.15, 100, -50, 50);
    }
  }
#else

  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spill",     &event.spill,     "spill/I");
  //Trig
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  // tree->Branch("cftadc_h",   event.cftadc_h,  Form("cftadc_h[%d][%d]/I", NumOfPlaneCFT, NumOfSegCFT_PHI4));
  // tree->Branch("cftadc_l",   event.cftadc_l,  Form("cftadc_l[%d][%d]/I", NumOfPlaneCFT, NumOfSegCFT_PHI4));
  tree->Branch("cfttdc",     event.cfttdc,    Form("cfttdc[%d][%d]/I", NumOfPlaneCFT, NumOfSegCFT_PHI4));

  tree->Branch("cftadc_cor_h",   event.cftadc_cor_h,  Form("cftadc_cor_h[%d][%d]/I", NumOfPlaneCFT, NumOfSegCFT_PHI4));
  tree->Branch("cftadc_cor_l",   event.cftadc_cor_l,  Form("cftadc_cor_l[%d][%d]/I", NumOfPlaneCFT, NumOfSegCFT_PHI4));


  tree->Branch("ntCFT",   &event.ntCFT,  "ntCFT/I");
  tree->Branch("x0",   event.x0,  "x0[ntCFT]/D");
  tree->Branch("y0",   event.y0,  "y0[ntCFT]/D");
  tree->Branch("u0",   event.u0,  "u0[ntCFT]/D");
  tree->Branch("v0",   event.v0,  "v0[ntCFT]/D");
  tree->Branch("dist",   event.dist,  "dist[ntCFT]/D");
  tree->Branch("distMeanx",   event.distMeanx,  "distMeanx[ntCFT]/D");
  tree->Branch("distMeany",   event.distMeany,  "distMeany[ntCFT]/D");
  tree->Branch("distMeanz",   event.distMeanz,  "distMeanz[ntCFT]/D");
  tree->Branch("track_theta", event.track_theta,"track_theta[ntCFT]/D");
  tree->Branch("track_phi",   event.track_phi,  "track_phi[ntCFT]/D");

  tree->Branch("nVer_product",   &event.nVer_product, "nVer_product/I");
  tree->Branch("dist_product",   event.dist_product,  "dist_product[nVer_product]/D");
  tree->Branch("xtar_product",   event.xtar_product,  "xtar_product[nVer_product]/D");
  tree->Branch("ytar_product",   event.ytar_product,  "ytar_product[nVer_product]/D");
  tree->Branch("ztar_product",   event.ztar_product,  "ztar_product[nVer_product]/D");

#endif // For PHC

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
