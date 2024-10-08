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
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "HodoHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "S2sLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "DCGeomMan.hh"

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

  Int_t bh2nhits;
  Int_t bh2hitpat[MaxHits];
  Double_t bh2ua[NumOfSegBH2];
  Double_t bh2ut[NumOfSegBH2][MaxDepth];
  Double_t bh2da[NumOfSegBH2];
  Double_t bh2dt[NumOfSegBH2][MaxDepth];

  Int_t bacnhits;
  Int_t bachitpat[MaxHits];
  Double_t baca[NumOfSegBAC];
  Double_t bact[NumOfSegBAC][MaxDepth];

  Int_t bac1nhits;
  Int_t bac2nhits;


  Int_t tofnhits;
  Int_t tofhitpat[MaxHits];
  Double_t tofua[NumOfSegTOF];
  Double_t tofut[NumOfSegTOF][MaxDepth];
  Double_t tofda[NumOfSegTOF];
  Double_t tofdt[NumOfSegTOF][MaxDepth];

  ////////// Normalized
  Double_t bh2mt[NumOfSegBH2][MaxDepth];
  Double_t bh2cmt[NumOfSegBH2][MaxDepth];
  Double_t bh2utime[NumOfSegBH2][MaxDepth];
  Double_t bh2uctime[NumOfSegBH2][MaxDepth];
  Double_t bh2dtime[NumOfSegBH2][MaxDepth];
  Double_t bh2dctime[NumOfSegBH2][MaxDepth];
  Double_t bh2hitpos[NumOfSegBH2][MaxDepth];
  Double_t bh2de[NumOfSegBH2];
  Double_t bh2ude[NumOfSegBH2];
  Double_t bh2dde[NumOfSegBH2];

  Double_t bacmt[NumOfSegBAC][MaxDepth];
  Double_t bacde[NumOfSegBAC];

  Double_t t0[NumOfSegBH2][MaxDepth];
  Double_t ct0[NumOfSegBH2][MaxDepth];

  Double_t tofmt[NumOfSegTOF][MaxDepth];
  Double_t tofde[NumOfSegTOF];
  Double_t tofude[NumOfSegTOF];
  Double_t tofdde[NumOfSegTOF];
  Double_t tofctu[NumOfSegTOF][MaxDepth];
  Double_t tofctd[NumOfSegTOF][MaxDepth];
  Double_t tofcmt[NumOfSegTOF][MaxDepth];


  // Time0
  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum     = 0;
  spill     = 0;
  bacnhits  = 0;
  bac1nhits = 0;
  bac2nhits = 0;
  bh2nhits  = 0;
  tofnhits  = 0;
  Time0Seg  = qnan;
  deTime0   = qnan;
  Time0     = qnan;
  CTime0    = qnan;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<MaxHits; ++it){
    bh2hitpat[it]   = -1;
    bachitpat[it]   = -1;
    tofhitpat[it]   = -1;
  }

  for(Int_t it=0; it<NumOfSegBH2; ++it){
    bh2ua[it]  = qnan;
    bh2da[it]  = qnan;
    bh2de[it]  = qnan;
    bh2dde[it] = qnan;
    bh2ude[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      bh2ut[it][m]     = qnan;
      bh2dt[it][m]     = qnan;
      bh2mt[it][m]     = qnan;
      t0[it][m]        = qnan;
      ct0[it][m]       = qnan;
      bh2utime[it][m]  = qnan;
      bh2dtime[it][m]  = qnan;
      bh2uctime[it][m] = qnan;
      bh2dctime[it][m] = qnan;
      bh2hitpos[it][m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegBAC; it++){
    baca[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      bact[it][m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegTOF; it++){
    tofua[it]  = qnan;
    tofda[it]  = qnan;
    tofde[it]  = qnan;
    tofude[it] = qnan;
    tofdde[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tofut[it][m]  = qnan;
      tofdt[it][m]  = qnan;
      tofmt[it][m]  = qnan;
      tofctu[it][m] = qnan;
      tofctd[it][m] = qnan;
      tofcmt[it][m] = qnan;
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

  Int_t    nhBac;
  Int_t    csBac[NumOfSegBAC*MaxDepth];
  Double_t BacSeg[NumOfSegBAC*MaxDepth];
  Double_t tBac[NumOfSegBAC*MaxDepth];
  Double_t deBac[NumOfSegBAC*MaxDepth];

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
  nhBac    = 0;
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

  for(Int_t it=0; it<NumOfSegBAC; it++){
    BacSeg[it] = qnan;
    deBac[it]  = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tBac[MaxDepth*it+m] = qnan;
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
  BACHid  = 30000,
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
  // static const auto MinTdcBH2 = gUser.GetParameter("TdcBH2", 0);
  // static const auto MaxTdcBH2 = gUser.GetParameter("TdcBH2", 1);
  // static const auto MinTdcBAC = gUser.GetParameter("TdcBAC", 0);
  // static const auto MaxTdcBAC = gUser.GetParameter("TdcBAC", 1);
  static const auto MinTdcTOF = gUser.GetParameter("TdcTOF", 0);
  static const auto MaxTdcTOF = gUser.GetParameter("TdcTOF", 1);
#if HodoHitPos
  static const auto PropVelBH2 = gUser.GetParameter("PropagationBH2");
#endif

  RawData rawData;
  HodoAnalyzer hodoAna(rawData);


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
      HF1(10, seg);
      HF1(10+seg, tdc);
    }
  }

  if(trigger_flag[trigger::kSpillOnEnd] || trigger_flag[trigger::kSpillOffEnd])
    return true;

  HF1(1, 1);

#if 0 // BH2, BAC
  ///// BH2
  rawData.DecodeHits("BH2");
  {
    Int_t bh2_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("BH2");
    Int_t nh = cont.size();
    HF1(BH2Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(BH2Hid +1, seg-0.5);
      // Up
      Int_t Au = hit->GetAdcUp();
      HF1(BH2Hid +100*seg +1, Au);
      event.bh2ua[seg-1] = Au;
      Bool_t is_hit_u = false;
      Int_t m_u = 0;
      for(const auto& T: hit->GetArrayTdcUp()){
        HF1(BH2Hid +100*seg +3, T);
        if(m_u < MaxDepth) event.bh2ut[seg-1][m_u++] = T;
        if(MinTdcBH2 < T && T < MaxTdcBH2) is_hit_u = true;
      }
      if(is_hit_u) HF1(BH2Hid +100*seg +5, Au);
      else         HF1(BH2Hid +100*seg +7, Au);
      // Down
      Int_t Ad = hit->GetAdcDown();
      HF1(BH2Hid +100*seg +2, Double_t(Ad));
      event.bh2da[seg-1] = Ad;
      Bool_t is_hit_d = false;
      Int_t m_d = 0;
      for(const auto& T: hit->GetArrayTdcDown()){
        HF1(BH2Hid +100*seg +4, Double_t(T));
        if(m_d < MaxDepth) event.bh2dt[seg-1][m_d++] = T;
        if(MinTdcBH2 < T && T < MaxTdcBH2) is_hit_d = true;
      }
      if(is_hit_d) HF1(BH2Hid +100*seg +6, Ad);
      else         HF1(BH2Hid +100*seg +8, Ad);
      // HitPat
      if(is_hit_u || is_hit_d){
        ++nh1; HF1(BH2Hid +3, seg-0.5);
      }
      if(is_hit_u && is_hit_d){
        event.bh2hitpat[bh2_nhits++] = seg;
        ++nh2; HF1(BH2Hid +5, seg-0.5);
      }
    }
    HF1(BH2Hid +2, nh1); HF1(BH2Hid +4, nh2);
    event.bh2nhits = bh2_nhits;
  }

  ///// BAC
  rawData.DecodeHits("BAC");
  {
    Int_t bac_nhits = 0;
    Int_t bac1_nhits = 0;
    Int_t bac2_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("BAC");
    Int_t nh = cont.size();
    HF1(BACHid, nh);
    Int_t nh1 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(BACHid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      HF1(BACHid+100*seg+1, A);
      event.baca[seg-1] = A;
      Bool_t is_hit = false;
      Int_t m = 0;
      for(const auto& T: hit->GetArrayTdcLeading()){
        HF1(BACHid+100*seg+3, T);
        if(m < MaxDepth) event.bact[seg-1][m++] = T;
        if(MinTdcBAC < T && T < MaxTdcBAC){
	  is_hit = true;
	  if(seg==1)++bac1_nhits;
	  if(seg==2)++bac2_nhits;
	}
      }
      if(is_hit) HF1(BACHid+100*seg+5, A);
      else       HF1(BACHid+100*seg+7, A);
      // Hitpat
      if(is_hit){
        event.bachitpat[bac_nhits++] = seg;
        ++nh1; HF1(BACHid+3, seg-0.5);
      }
    }
    HF1(BACHid+2, nh1);
    event.bacnhits = bac_nhits;
    event.bac1nhits = bac1_nhits;
    event.bac2nhits = bac2_nhits;
  }
#endif // BH2, BAC

  ///// TOF
  rawData.DecodeHits("TOF");
  {
    Int_t tof_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("TOF");
    Int_t nh = cont.size();
    HF1(TOFHid, Double_t(nh));
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId();
      HF1(TOFHid+1, seg+0.5);

      // Up
      Int_t Au = hit->GetAdcUp();
      HF1(TOFHid+100*(seg+1)+1, Au);
      event.tofua[seg] = Au;
      dst.tofua[seg]   = Au;

      Bool_t is_hit_u = false;
      Int_t m_u = 0;
      for(const auto& T: hit->GetArrayTdcUp()){
        HF1(TOFHid +100*(seg+1) +3, T);
        if(m_u < MaxDepth){
          event.tofut[seg][m_u] = T;
          dst.tofut[seg][m_u]   = T;
          ++m_u;
        }
        if(MinTdcTOF < T && T < MaxTdcTOF) is_hit_u = true;
      }
      if(is_hit_u) HF1(TOFHid+100*(seg+1)+5, Au);
      else         HF1(TOFHid+100*(seg+1)+7, Au);

      // Down
      Int_t Ad = hit->GetAdcDown();
      HF1(TOFHid+100*(seg+1)+2, Ad);
      event.tofda[seg] = Ad;
      dst.tofda[seg]   = Ad;

      Bool_t is_hit_d = false;
      Int_t m_d = 0;
      for(const auto& T: hit->GetArrayTdcDown()){
        HF1(TOFHid +100*(seg+1) +4, T);
        if(m_d < MaxDepth){
          event.tofdt[seg][m_d] = T;
          dst.tofdt[seg][m_d]   = T;
          ++m_d;
        }
        if(MinTdcTOF < T && T < MaxTdcTOF)  is_hit_d = true;
      }
      if(is_hit_d) HF1(TOFHid+100*(seg+1)+6, Ad);
      else         HF1(TOFHid+100*(seg+1)+8, Ad);

      // Hitpat
      if(is_hit_u || is_hit_d){
        ++nh1; HF1(TOFHid+3, seg+0.5);
      }
      if(is_hit_u && is_hit_d){
        event.tofhitpat[tof_nhits++] = seg;
        ++nh2; HF1(TOFHid+5, seg+0.5);
      }
    }
    HF1(TOFHid+2, nh1); HF1(TOFHid+4, nh2);
    event.tofnhits = tof_nhits;
  }


  //*****************  Normalzed Data  *****************

#if 0 // Normalized

#if 0 // BH2, BAC
  // BH2
  hodoAna.DecodeHits<BH2Hit>("BH2");
  // hodoAna.TimeCut("BH2", -2, 2);
  {
    Int_t nh = hodoAna.GetNHits("BH2");
    HF1(BH2Hid+10, Double_t(nh));
    Int_t nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit<BH2Hit>("BH2", i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId()+1;

      Int_t n_mhit = hit->GetEntries();
      for(Int_t m=0; m<n_mhit; ++m){
        HF1(BH2Hid+11, seg-0.5);
        Double_t au  = hit->GetAUp(),   ad  = hit->GetADown();
        Double_t tu  = hit->GetTUp(m),  td  = hit->GetTDown(m);
        Double_t ctu = hit->GetCTUp(m), ctd = hit->GetCTDown(m);
        Double_t mt  = hit->MeanTime(m),cmt = hit->CMeanTime(m);
        Double_t de  = hit->DeltaE();
        Double_t ude  = hit->UDeltaE();
        Double_t dde  = hit->DDeltaE();
        Double_t ut0  = hit->UTime0(m), dt0  = hit->DTime0(m);
        Double_t uct0 = hit->UCTime0(m),dct0 = hit->DCTime0(m);
        Double_t t0   = hit->Time0(m),  ct0  = hit->CTime0(m);
#if HodoHitPos
	Double_t ctdiff = hit->TimeDiff(m);
	event.bh2hitpos[seg-1][m] = 0.5*PropVelBH2*ctdiff;
#endif
        event.bh2mt[seg-1][m] = mt;
        event.bh2utime[seg-1][m] = tu;
        event.bh2dtime[seg-1][m] = td;
        event.bh2uctime[seg-1][m] = ctu;
        event.bh2dctime[seg-1][m] = ctd;
        event.bh2de[seg-1]    = de;
        event.bh2ude[seg-1]    = ude;
        event.bh2dde[seg-1]    = dde;
        event.t0[seg-1][m]    = t0;
        event.ct0[seg-1][m]   = ct0;

	HF1(BH2Hid+100*seg+11, tu);      HF1(BH2Hid+100*seg+12, td);
        HF1(BH2Hid+100*seg+13, mt);
        HF1(BH2Hid+100*seg+17, ctu);     HF1(BH2Hid+100*seg+18, ctd);
        HF1(BH2Hid+100*seg+19, cmt);     HF1(BH2Hid+100*seg+20, ctu-ctd);
        HF1(BH2Hid+100*seg+21, ut0);     HF1(BH2Hid+100*seg+22, dt0);
        HF1(BH2Hid+100*seg+23, uct0);    HF1(BH2Hid+100*seg+24, dct0);
        HF1(BH2Hid+100*seg+25, t0);      HF1(BH2Hid+100*seg+26, ct0);
        HF2(BH2Hid+100*seg+27, tu, au);  HF2(BH2Hid+100*seg+28, td, ad);
        HF2(BH2Hid+100*seg+29, ctu, au); HF2(BH2Hid+100*seg+30, ctd, ad);

        // HF1(BH2Hid+100*seg+11, tu); HF1(BH2Hid+100*seg+13, mt);
        // HF1(BH2Hid+100*seg+14, au); HF1(BH2Hid+100*seg+16, de);
        // HF1(BH2Hid+100*seg+17, ctu); HF1(BH2Hid+100*seg+19, cmt);
        // HF2(BH2Hid+100*seg+21, tu, au); HF2(BH2Hid+100*seg+23, ctu, au);
        // HF1(BH2Hid+12, cmt);

        if(m == 0){
          HF1(BH2Hid+100*seg+14, au);	   HF1(BH2Hid+100*seg+15, ad);
          HF1(BH2Hid+100*seg+16, de);    HF1(BH2Hid+13, de);
        }

        if(de>0.5){
          ++nh2; HF1(BH2Hid+15, seg-0.5); HF1(BH2Hid+16, cmt);
        }
      }
    }//for(i)
    HF1(BH2Hid+14, Double_t(nh2));
    for(Int_t i1=0; i1<nh; ++i1){
      const auto& hit1 = hodoAna.GetHit<BH2Hit>("BH2", i1);
      if(!hit1 || hit1->DeltaE()<=0.5) continue;
      Int_t seg1 = hit1->SegmentId()+1;
      for(Int_t i2=0; i2<nh; ++i2){
        if(i1==i2) continue;
        const auto& hit2 = hodoAna.GetHit<BH2Hit>("BH2", i2);
        if(!hit2 || hit2->DeltaE()<=0.5) continue;
        Int_t seg2 = hit2->SegmentId()+1;

        if(1 == hit1->GetEntries() && 1 == hit2->GetEntries()){
          Double_t ct1 = hit1->CMeanTime(), ct2 = hit2->CMeanTime();
          HF2(BH2Hid+21, seg1-0.5, seg2-0.5);
          HF2(BH2Hid+22, ct1, ct2);
          HF1(BH2Hid+23, ct2-ct1);
          if(std::abs(ct2-ct1)<2.0){
            HF2(BH2Hid+24, seg1-0.5, seg2-0.5);
          }
        }
      }//for(i2)
    }//for(i1)

    Int_t nc=hodoAna.GetNClusters("BH2");
    HF1(BH2Hid+30, Double_t(nc));
    Int_t nc2=0;

    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster("BH2", i);
      if(!cl) continue;
      Int_t cs=cl->ClusterSize();
      Double_t ms = cl->MeanSeg()+1;
      Double_t cmt= cl->CMeanTime();
      Double_t de = cl->DeltaE();
      // Double_t mt = cl->MeanTime();
      HF1(BH2Hid+31, Double_t(cs));
      HF1(BH2Hid+32, ms-0.5);
      HF1(BH2Hid+33, cmt); HF1(BH2Hid+34, de);
      if(de>0.5){
        ++nc2; HF1(BH2Hid+36, cmt);
      }

      for(Int_t i2=0; i2<nc; ++i2){
        if(i2==i) continue;
        const auto& cl2 = hodoAna.GetCluster("BH2", i2);
        if(!cl2) continue;
        Double_t ms2=cl2->MeanSeg()+1, cmt2=cl2->CMeanTime(),
          de2=cl2->DeltaE();
        if(de<=0.5 || de2<=0.5) continue;
        HF2(BH2Hid+41, ms-0.5, ms2-0.5);
        HF2(BH2Hid+42, cmt, cmt2);
        HF1(BH2Hid+43, cmt2-cmt);
        if(std::abs(cmt2-cmt)<2.0){
          HF2(BH2Hid+44, ms-0.5, ms2-0.5);
        }
      }//for(i2)
    }//for(i)
    HF1(BH2Hid+35, Double_t(nc2));

    const auto& cl_time0 = hodoAna.GetTime0BH2Cluster();
    if(cl_time0){
      event.Time0Seg = cl_time0->MeanSeg()+1;
      event.deTime0 = cl_time0->DeltaE();
      event.Time0 = cl_time0->Time0();
      event.CTime0 = cl_time0->CTime0();
      dst.Time0Seg = cl_time0->MeanSeg()+1;
      dst.deTime0 = cl_time0->DeltaE();
      dst.Time0 = cl_time0->Time0();
      dst.CTime0 = cl_time0->CTime0();
      HF1(100, cl_time0->MeanTime());
      HF1(101, cl_time0->CTime0());
    }

  }

#if 0
  // BH1-BH2 PHC
  {
    Int_t nh1 = hodoAna.GetNHits("BH1");
    Int_t nh2 = hodoAna.GetNHits("BH2");
    for(Int_t i2=0; i2<nh2; ++i2){
      const auto& hit2 = hodoAna.GetHit<BH2Hit>("BH2", i2);
      Int_t     seg2 = hit2->SegmentId()+1;
      Double_t  au2  = hit2->GetAUp(),  ad2  = hit2->GetADown();
      Int_t n_mhit2  = hit2->GetEntries();
      for(Int_t m2 = 0; m2<n_mhit2; ++m2){
        Double_t  tu2  = hit2->GetTUp(m2),  td2  = hit2->GetTDown(m2);
        Double_t  ctu2 = hit2->GetCTUp(m2), ctd2 = hit2->GetCTDown(m2);
        // Double_t  t0   = hit2->Time0();
        Double_t  ct0  = hit2->CTime0(m2);
        Double_t  tofs = ct0-(ctu2+ctd2)/2.;
        for(Int_t i1=0; i1<nh1; ++i1){
          const auto& hit1 = hodoAna.GetHit("BH1", i1);
          Int_t       seg1 = hit1->SegmentId()+1;
          Double_t    au1  = hit1->GetAUp(),  ad1  = hit1->GetADown();
          Int_t    n_mhit1 = hit1->GetEntries();
          for(Int_t m1 = 0; m1<n_mhit1; ++m1){
            Double_t    tu1  = hit1->GetTUp(m1),  td1  = hit1->GetTDown(m1);
            Double_t    ctu1 = hit1->GetCTUp(m1), ctd1 = hit1->GetCTDown(m1);
            Double_t    cmt1 = hit1->CMeanTime(m1);
            // if(trigger_flag[trigger::kBeamA] == 0) continue;
            // if(event.bacnhits > 0) continue;
            HF2(100*seg1+BH1Hid+81, au1, ct0-0.5*(ctu1+ctd1));
            HF2(100*seg1+BH1Hid+82, ad1, ct0-0.5*(ctu1+ctd1));
            HF2(100*seg1+BH1Hid+83, au1, ct0-tu1);
            HF2(100*seg1+BH1Hid+84, ad1, ct0-td1);
            HF2(100*seg2+BH2Hid+81, au2, (cmt1-tofs)-0.5*(ctu2+ctd2));
            HF2(100*seg2+BH2Hid+82, ad2, (cmt1-tofs)-0.5*(ctu2+ctd2));
            HF2(100*seg2+BH2Hid+83, au2, (cmt1-tofs)-tu2);
            HF2(100*seg2+BH2Hid+84, ad2, (cmt1-tofs)-td2);
          }// for(m1)
        }// for(bh1:seg)
      }// for(m2)
    }// for(bh2:seg)
  }
#endif

  // BAC
  hodoAna.DecodeHits("BAC");
  {
    Int_t nh=hodoAna.GetNHits("BAC");
    dst.nhBac = nh;
    HF1(BACHid+10, Double_t(nh));
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit("BAC", i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId()+1;
      HF1(BACHid+11, seg-0.5);
      Double_t a = hit->GetAUp();
      event.bacde[seg-1]  = a;
      Int_t n_mhit = hit->GetEntries();
      for(Int_t m=0; m<n_mhit; ++m){
        Double_t t = hit->GetTUp(m);
        Double_t ct = hit->GetCTUp(m);
        event.bacmt[i][m] = ct;
        HF1(BACHid+100*seg+11, t);
        HF1(BACHid+100*seg+12, a);
        HF1(BACHid+100*seg+13, ct);
      }
    }
  }
#endif // BH2, BAC

  // TOF
  hodoAna.DecodeHits("TOF");
  {
    Int_t nh = hodoAna.GetNHits("TOF");
    HF1(TOFHid+10, Double_t(nh));
    Int_t nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit("TOF", i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId();
      Int_t n_mhit = hit->GetEntries();
      for(Int_t m=0; m<n_mhit; ++m){
        HF1(TOFHid+11, seg+0.5);
        Double_t au   = hit->GetAUp(),   ad  = hit->GetADown();
        Double_t tu   = hit->GetTUp(),   td  = hit->GetTDown();
        Double_t ctu  = hit->GetCTUp(),  ctd = hit->GetCTDown();
        Double_t mt   = hit->MeanTime(), cmt = hit->CMeanTime();
        Double_t de   = hit->DeltaE();
        Double_t ude  = hit->UDeltaE();
        Double_t dde  = hit->DDeltaE();
        event.tofmt[seg][m]  = mt;
        event.tofde[seg]     = de;
        event.tofude[seg]    = ude;
        event.tofdde[seg]    = dde;
        event.tofctu[seg][m] = ctu;
        event.tofctd[seg][m] = ctd;
        event.tofcmt[seg][m] = cmt;
        HF1(TOFHid+100*seg+11, tu);      HF1(TOFHid+100*seg+12, td);
        HF1(TOFHid+100*seg+13, mt);
        HF1(TOFHid+100*seg+17, ctu);     HF1(TOFHid+100*seg+18, ctd);
        HF1(TOFHid+100*seg+19, cmt);     HF1(TOFHid+100*seg+20, ctu-ctd);
        //HF2(TOFHid+100*seg+21, tu, au);  HF2(TOFHid+100*seg+22, td, ad);
        //HF2(TOFHid+100*seg+23, ctu, au); HF2(TOFHid+100*seg+24, ctd, ad);
        HF1(TOFHid+12, cmt);

        dst.utTofSeg[seg][m] = tu;
        dst.dtTofSeg[seg][m] = td;
        dst.udeTofSeg[seg]   = au;
        dst.ddeTofSeg[seg]   = ad;

        if(m == 0){
          HF1(TOFHid+100*seg+14, au);    HF1(TOFHid+100*seg+15, ad);
          HF1(TOFHid+100*seg+16, de);    HF1(TOFHid+13, de);
        }

        if(de>0.5){
          HF1(TOFHid+15, seg+0.5);
          ++nh2;
        }
      }
    }

    HF1(TOFHid+14, Double_t(nh2));
    for(Int_t i1=0; i1<nh; ++i1){
      const auto& hit1 = hodoAna.GetHit("TOF", i1);
      if(!hit1 || hit1->DeltaE()<=0.5) continue;
      Int_t seg1 = hit1->SegmentId();
      for(Int_t i2=0; i2<nh; ++i2){
        if(i1==i2) continue;
        const auto& hit2 = hodoAna.GetHit("TOF", i2);
        if(!hit2 || hit2->DeltaE()<=0.5) continue;
        Int_t seg2 = hit2->SegmentId();

        if(1 == hit1->GetEntries() && 1 == hit2->GetEntries()){
          Double_t ct1 = hit1->CMeanTime(), ct2 = hit2->CMeanTime();
          HF2(TOFHid+21, seg1+0.5, seg2+0.5);
          HF2(TOFHid+22, ct1, ct2);
          HF1(TOFHid+23, ct2-ct1);
          if(std::abs(ct2-ct1)<3.0){
            HF2(TOFHid+24, seg1+0.5, seg2+0.5);
          }
        }//for(i2)
      }//for(i1)
    }

    Int_t nc = hodoAna.GetNClusters("TOF");
    HF1(TOFHid+30, Double_t(nc));
    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster("TOF", i);
      if(!cl) continue;
      Int_t cs = cl->ClusterSize();
      Double_t ms  = cl->MeanSeg();
      Double_t cmt = cl->CMeanTime();
      Double_t de  = cl->DeltaE();
      HF1(TOFHid+31, Double_t(cs));
      HF1(TOFHid+32, ms+0.5);
      HF1(TOFHid+33, cmt); HF1(TOFHid+34, de);
    }
  }
#endif // Normalized

#if 0 // Dst
  ////////// Dst
  {
    Int_t nc = hodoAna.GetNClusters("BH2");
    dst.nhBh2 = nc;
    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster<BH2Cluster>("BH2", i);
      if(!cl) continue;
      dst.csBh2[i]  = cl->ClusterSize();
      dst.Bh2Seg[i] = cl->MeanSeg()+1;
      dst.tBh2[i]   = cl->CMeanTime();
      dst.t0Bh2[i]  = cl->CTime0();
      dst.dtBh2[i]  = cl->TimeDiff();
      dst.deBh2[i]  = cl->DeltaE();
#if HodoHitPos
      dst.posBh2[i]  = 0.5*PropVelBH2*cl->TimeDiff();
#endif
    }
  }

  {
    Int_t nc = hodoAna.GetNClusters("BAC");
    dst.nhBac = nc;
    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster("BAC", i);
      if(!cl) continue;
      dst.csBac[i]  = cl->ClusterSize();
      dst.BacSeg[i] = cl->MeanSeg()+1;
      dst.tBac[i]   = cl->CMeanTime();
      dst.deBac[i]  = cl->DeltaE();
    }
  }

  {
    Int_t nc = hodoAna.GetNClusters("TOF");
    dst.nhTof = nc;
    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster("TOF", i);
      if(!cl) continue;
      dst.csTof[i]  = cl->ClusterSize();
      dst.TofSeg[i] = cl->MeanSeg();
      dst.tTof[i]   = cl->CMeanTime();
      dst.dtTof[i]  = cl->TimeDiff();
      dst.deTof[i]  = cl->DeltaE();
    }
  }
#endif // Dst

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
  HB1(10, "Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig));
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000);
  }

#if 0 // BH2, BAC
  // BH2
  HB1(BH2Hid +0, "#Hits BH2",        NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +1, "Hitpat BH2",       NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +2, "#Hits BH2(Tor)",   NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +3, "Hitpat BH2(Tor)",  NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +4, "#Hits BH2(Tand)",  NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +5, "Hitpat BH2(Tand)", NumOfSegBH2,   0., Double_t(NumOfSegBH2));

  for(Int_t i=1; i<=NumOfSegBH2; ++i){
    TString title1 = Form("BH2-%d UpAdc", i);
    TString title2 = Form("BH2-%d DownAdc", i);
    TString title3 = Form("BH2-%d UpTdc", i);
    TString title4 = Form("BH2-%d DownTdc", i);
    TString title5 = Form("BH2-%d UpAdc(w Tdc)", i);
    TString title6 = Form("BH2-%d DownAdc(w Tdc)", i);
    TString title7 = Form("BH2-%d UpAdc(w/o Tdc)", i);
    TString title8 = Form("BH2-%d DownAdc(w/o Tdc)", i);
    HB1(BH2Hid +100*i +1, title1, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +2, title2, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(BH2Hid +100*i +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(BH2Hid +100*i +5, title5, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +6, title6, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +7, title7, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +8, title8, NbinAdc,   MinAdc,   MaxAdc);
  }
  HB1(BH2Hid +10, "#Hits BH2[Hodo]",  NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +11, "Hitpat BH2[Hodo]", NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +12, "CMeanTime BH2", 200, -10., 10.);
  HB1(BH2Hid +13, "dE BH2", 200, -0.5, 4.5);
  HB1(BH2Hid +14, "#Hits BH2[HodoGood]",  NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +15, "Hitpat BH2[HodoGood]", NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +16, "CMeanTime BH2[HodoGood]", 200, -10., 10.);

  for(Int_t i=1; i<=NumOfSegBH2; ++i){
    TString title11 = Form("BH2-%d Up Time", i);
    TString title12 = Form("BH2-%d Down Time", i);
    TString title13 = Form("BH2-%d MeanTime", i);
    TString title14 = Form("BH2-%d Up dE", i);
    TString title15 = Form("BH2-%d Down dE", i);
    TString title16 = Form("BH2-%d dE", i);
    TString title17 = Form("BH2-%d Up CTime", i);
    TString title18 = Form("BH2-%d Down CTime", i);
    TString title19 = Form("BH2-%d CMeanTime", i);
    TString title20 = Form("BH2-%d Tup-Tdown", i);
    TString title21 = Form("BH2-%d Up Time0", i);
    TString title22 = Form("BH2-%d Down Time0", i);
    TString title23 = Form("BH2-%d Up CTime", i);
    TString title24 = Form("BH2-%d Down CTime", i);
    TString title25 = Form("BH2-%d MeanTime0", i);
    TString title26 = Form("BH2-%d CMeanTime0", i);
    TString title27 = Form("BH2-%d Up dE%%Time", i);
    TString title28 = Form("BH2-%d Down dE%%Time", i);
    TString title29 = Form("BH2-%d Up dE%%CTime", i);
    TString title30 = Form("BH2-%d Down dE%%CTime", i);
    HB1(BH2Hid +100*i +11, title11, 200, -10., 10.);
    HB1(BH2Hid +100*i +12, title12, 200, -10., 10.);
    HB1(BH2Hid +100*i +13, title13, 200, -10., 10.);
    HB1(BH2Hid +100*i +14, title14, 200, -0.5, 4.5);
    HB1(BH2Hid +100*i +15, title15, 200, -0.5, 4.5);
    HB1(BH2Hid +100*i +16, title16, 200, -0.5, 4.5);
    HB1(BH2Hid +100*i +17, title17, 200, -10., 10.);
    HB1(BH2Hid +100*i +18, title18, 200, -10., 10.);
    HB1(BH2Hid +100*i +19, title19, 200, -10., 10.);
    HB1(BH2Hid +100*i +20, title20, 200, -5.0, 5.0);
    HB1(BH2Hid +100*i +21, title21, 200, -10., 10.);
    HB1(BH2Hid +100*i +22, title22, 200, -10., 10.);
    HB1(BH2Hid +100*i +23, title23, 200, -10., 10.);
    HB1(BH2Hid +100*i +24, title24, 200, -10., 10.);
    HB1(BH2Hid +100*i +25, title25, 200, -10., 10.);
    HB1(BH2Hid +100*i +26, title26, 200, -10., 10.);
    HB2(BH2Hid +100*i +27, title27, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH2Hid +100*i +28, title28, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH2Hid +100*i +29, title29, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH2Hid +100*i +30, title30, 100, -10., 10., 100, -0.5, 4.5);
  }

  HB2(BH2Hid +21, "BH2HitPat%BH2HitPat[HodoGood]", NumOfSegBH2,   0., Double_t(NumOfSegBH2),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB2(BH2Hid +22, "CMeanTimeBH2%CMeanTimeBH2[HodoGood]",
      100, -2.5, 2.5, 100, -2.5, 2.5);
  HB1(BH2Hid +23, "TDiff BH2[HodoGood]", 200, -10., 10.);
  HB2(BH2Hid +24, "BH2HitPat%BH2HitPat[HodoGood2]", NumOfSegBH2,   0., Double_t(NumOfSegBH2),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));

  HB1(BH2Hid +30, "#Clusters BH2", NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +31, "ClusterSize BH2", 5, 0., 5.);
  HB1(BH2Hid +32, "HitPat Cluster BH2", 2*NumOfSegBH2, 0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +33, "CMeamTime Cluster BH2", 200, -10., 10.);
  HB1(BH2Hid +34, "DeltaE Cluster BH2", 100, -0.5, 4.5);
  HB1(BH2Hid +35, "#Clusters BH2(ADCGood)", NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +36, "CMeamTime Cluster BH2(ADCGood)", 200, -10., 10.);

  HB2(BH2Hid +41, "BH2ClP%BH2ClP", NumOfSegBH2,   0., Double_t(NumOfSegBH2),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB2(BH2Hid +42, "CMeanTimeBH2%CMeanTimeBH2[Cluster]",
      100, -2.5, 2.5, 100, -2.5, 2.5);
  HB1(BH2Hid +43, "TDiff BH2[Cluster]", 200, -10., 10.);
  HB2(BH2Hid +44, "BH2ClP%BH2ClP(ADCGood)", NumOfSegBH2,   0., Double_t(NumOfSegBH2),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));

#if 0
  // BH1-BH2 PHC
  for(Int_t i=1; i<=NumOfSegBH1; ++i){
    TString title1 = Form("BH1-%dU CT-TOF%%dE", i);
    TString title2 = Form("BH1-%dD CT-TOF%%dE", i);
    TString title3 = Form("BH1-%dU  T-TOF%%dE", i);
    TString title4 = Form("BH1-%dD  T-TOF%%dE", i);
    HB2(BH1Hid +100*i +81, title1, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH1Hid +100*i +82, title2, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH1Hid +100*i +83, title3, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH1Hid +100*i +84, title4, 200, -0.5, 4.5, 200, -10., 10.);
  }
  for(Int_t i=1; i<=NumOfSegBH2; ++i){
    TString title1 = Form("BH2-%dU CT-TOF%%dE", i);
    TString title2 = Form("BH2-%dD CT-TOF%%dE", i);
    TString title3 = Form("BH2-%dU  T-TOF%%dE", i);
    TString title4 = Form("BH2-%dD  T-TOF%%dE", i);
    HB2(BH2Hid +100*i +81, title1, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH2Hid +100*i +82, title2, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH2Hid +100*i +83, title3, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH2Hid +100*i +84, title4, 200, -0.5, 4.5, 200, -10., 10.);
  }
#endif

  // BTOF
  HB1(100, "BH2 MeanTime0", 400, -4, 4);
  HB1(101, "CTime0", 400, -4, 4);

  // BAC
  HB1(BACHid +0, "#Hits BAC",        NumOfSegBAC+1, 0., Double_t(NumOfSegBAC+1));
  HB1(BACHid +1, "Hitpat BAC",       NumOfSegBAC,   0., Double_t(NumOfSegBAC));
  HB1(BACHid +2, "#Hits BAC(Tor)",   NumOfSegBAC+1, 0., Double_t(NumOfSegBAC+1));
  HB1(BACHid +3, "Hitpat BAC(Tor)",  NumOfSegBAC,   0., Double_t(NumOfSegBAC));
  HB1(BACHid +4, "#Hits BAC(Tand)",  NumOfSegBAC+1, 0., Double_t(NumOfSegBAC+1));
  HB1(BACHid +5, "Hitpat BAC(Tand)", NumOfSegBAC,   0., Double_t(NumOfSegBAC));
  for(Int_t i=1; i<=NumOfSegBAC; ++i){
    TString title1 = Form("BAC-%d UpAdc", i);
    TString title3 = Form("BAC-%d UpTdc", i);
    TString title5 = Form("BAC-%d UpAdc(w Tdc)", i);
    TString title7 = Form("BAC-%d UpAdc(w/o Tdc)", i);
    HB1(BACHid +100*i +1, title1, NbinTdc, MinTdc, MaxTdc);
    HB1(BACHid +100*i +3, title3, NbinTdc, MinTdc, MaxTdc);
    HB1(BACHid +100*i +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(BACHid +100*i +7, title7, NbinAdc, MinAdc, MaxAdc);
  }
  HB1(BACHid +10, "#Hits BAC[Hodo]",     NumOfSegBAC+1, 0., Double_t(NumOfSegBAC+1));
  HB1(BACHid +11, "Hitpat BAC[Hodo]",    NumOfSegBAC,   0., Double_t(NumOfSegBAC));
  for(Int_t i=1; i<=NumOfSegBAC; ++i){
    TString title1 = Form("BAC-%d Time", i);
    TString title3 = Form("BAC-%d dE", i);
    TString title5 = Form("BAC-%d CTime", i);
    HB1(BACHid +100*i +11, title1, 500, -5., 45.);
    HB1(BACHid +100*i +12, title3, 200, -0.5, 4.5);
    HB1(BACHid +100*i +13, title5, 500, -5., 45.);
  }
#endif // BH2, BAC

  // TOF
  HB1(TOFHid +0, "#Hits TOF",        NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +1, "Hitpat TOF",       NumOfSegTOF,   0., Double_t(NumOfSegTOF));
  HB1(TOFHid +2, "#Hits TOF(Tor)",   NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +3, "Hitpat TOF(Tor)",  NumOfSegTOF,   0., Double_t(NumOfSegTOF));
  HB1(TOFHid +4, "#Hits TOF(Tand)",  NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +5, "Hitpat TOF(Tand)", NumOfSegTOF,   0., Double_t(NumOfSegTOF));

  for(Int_t i=0; i<NumOfSegTOF; ++i){
    TString title1 = Form("TOF-%d UpAdc", i);
    TString title2 = Form("TOF-%d DownAdc", i);
    TString title3 = Form("TOF-%d UpTdc", i);
    TString title4 = Form("TOF-%d DownTdc", i);
    TString title5 = Form("TOF-%d UpAdc(w Tdc)", i);
    TString title6 = Form("TOF-%d DownAdc(w Tdc)", i);
    TString title7 = Form("TOF-%d UpAdc(w/o Tdc)", i);
    TString title8 = Form("TOF-%d DownAdc(w/o Tdc)", i);
    HB1(TOFHid +100*(i+1) +1, title1, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*(i+1) +2, title2, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*(i+1) +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(TOFHid +100*(i+1) +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(TOFHid +100*(i+1) +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*(i+1) +6, title6, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*(i+1) +7, title7, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*(i+1) +8, title8, NbinAdc, MinAdc, MaxAdc);
  }

#if 0 // TOF Normalized
  HB1(TOFHid +10, "#Hits Tof[Hodo]",  NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +11, "Hitpat Tof[Hodo]", NumOfSegTOF,   0., Double_t(NumOfSegTOF));
  HB1(TOFHid +12, "CMeanTime Tof", 500, -5., 45.);
  HB1(TOFHid +13, "dE Tof", 200, -0.5, 4.5);
  HB1(TOFHid +14, "#Hits Tof[HodoGood]",  NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +15, "Hitpat Tof[HodoGood]", NumOfSegTOF,   0., Double_t(NumOfSegTOF));

  for(Int_t i=0; i<NumOfSegTOF; ++i){
    TString title11 = Form("TOF-%d Up Time", i);
    TString title12 = Form("TOF-%d Down Time", i);
    TString title13 = Form("TOF-%d MeanTime", i);
    TString title14 = Form("TOF-%d Up dE", i);
    TString title15 = Form("TOF-%d Down dE", i);
    TString title16 = Form("TOF-%d dE", i);
    TString title17 = Form("TOF-%d Up CTime", i);
    TString title18 = Form("TOF-%d Down CTime", i);
    TString title19 = Form("TOF-%d CMeanTime", i);
    TString title20 = Form("TOF-%d Tup-Tdown", i);
    TString title21 = Form("TOF-%d dE (w/ TOF-HT)", i);
    HB1(TOFHid +100*(i+1) +11, title11, 500, -5., 45.);
    HB1(TOFHid +100*(i+1) +12, title12, 500, -5., 45.);
    HB1(TOFHid +100*(i+1) +13, title13, 500, -5., 45.);
    HB1(TOFHid +100*(i+1) +14, title14, 200, -0.5, 4.5);
    HB1(TOFHid +100*(i+1) +15, title15, 200, -0.5, 4.5);
    HB1(TOFHid +100*(i+1) +16, title16, 200, -0.5, 4.5);
    HB1(TOFHid +100*(i+1) +17, title17, 500, -5., 45.);
    HB1(TOFHid +100*(i+1) +18, title18, 500, -5., 45.);
    HB1(TOFHid +100*(i+1) +19, title19, 500, -5., 45.);
    HB1(TOFHid +100*(i+1) +20, title20, 200, -10.0, 10.0);
    HB1(TOFHid +100*(i+1) +21, title21, 200, -0.5, 4.5);
  }

  HB2(TOFHid +21, "TofHitPat%TofHitPat[HodoGood]", NumOfSegTOF,   0., Double_t(NumOfSegTOF),
      NumOfSegTOF,   0., Double_t(NumOfSegTOF));
  HB2(TOFHid +22, "CMeanTimeTof%CMeanTimeTof[HodoGood]",
      120, 10., 40., 120, 10., 40.);
  HB1(TOFHid +23, "TDiff Tof[HodoGood]", 200, -10., 10.);
  HB2(TOFHid +24, "TofHitPat%TofHitPat[HodoGood2]", NumOfSegTOF,   0., Double_t(NumOfSegTOF),
      NumOfSegTOF,   0., Double_t(NumOfSegTOF));

  HB1(TOFHid +30, "#Clusters Tof", NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +31, "ClusterSize Tof", 5, 0., 5.);
  HB1(TOFHid +32, "HitPat Cluster Tof", 2*NumOfSegTOF, 0., Double_t(NumOfSegTOF));
  HB1(TOFHid +33, "CMeamTime Cluster Tof", 500, -5., 45.);
  HB1(TOFHid +34, "DeltaE Cluster Tof", 100, -0.5, 4.5);
#endif // TOF Normalized

  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spill",     &event.spill,     "spill/I");
  //Trig
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  // //BH2
  // tree->Branch("bh2nhits",   &event.bh2nhits,    "bh2nhits/I");
  // tree->Branch("bh2hitpat",   event.bh2hitpat,   Form("bh2hitpat[%d]/I", NumOfSegBH2));
  // tree->Branch("bh2ua",       event.bh2ua,       Form("bh2ua[%d]/D", NumOfSegBH2));
  // tree->Branch("bh2ut",       event.bh2ut,       Form("bh2ut[%d][%d]/D", NumOfSegBH2, MaxDepth));
  // tree->Branch("bh2da",       event.bh2da,       Form("bh2da[%d]/D", NumOfSegBH2));
  // tree->Branch("bh2dt",       event.bh2dt,       Form("bh2dt[%d][%d]/D", NumOfSegBH2, MaxDepth));
  // //BAC
  // tree->Branch("bacnhits",   &event.bacnhits,   "bacnhits/I");
  // tree->Branch("bachitpat",   event.bachitpat,  Form("bachitpat[%d]/I", NumOfSegBAC));
  // tree->Branch("baca",        event.baca,       Form("baca[%d]/D", NumOfSegBAC));
  // tree->Branch("bact",        event.bact,       Form("bact[%d][%d]/D", NumOfSegBAC, MaxDepth));
  // tree->Branch("bac1nhits",   &event.bac1nhits,   "bac1nhits/I");
  // tree->Branch("bac2nhits",   &event.bac2nhits,   "bac2nhits/I");
  //TOF
  tree->Branch("tofnhits",   &event.tofnhits,   "tofnhits/I");
  tree->Branch("tofhitpat",   event.tofhitpat,  Form("tofhitpat[%d]/I", NumOfSegTOF));
  tree->Branch("tofua",       event.tofua,      Form("tofua[%d]/D", NumOfSegTOF));
  tree->Branch("tofut",       event.tofut,      Form("tofut[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofda",       event.tofda,      Form("tofda[%d]/D", NumOfSegTOF));
  tree->Branch("tofdt",       event.tofdt,      Form("tofdt[%d][%d]/D", NumOfSegTOF, MaxDepth));

#if 0 // Normalized, Dst
  //Normalized data
  tree->Branch("bh2mt",     event.bh2mt,     Form("bh2mt[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2utime",     event.bh2utime,     Form("bh2utime[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2dtime",     event.bh2dtime,     Form("bh2dtime[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2uctime",     event.bh2uctime,     Form("bh2uctime[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2dctime",     event.bh2dctime,     Form("bh2dctime[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2de",     event.bh2de,     Form("bh2de[%d]/D", NumOfSegBH2));
  tree->Branch("bh2ude",     event.bh2ude,     Form("bh2ude[%d]/D", NumOfSegBH2));
  tree->Branch("bh2dde",     event.bh2dde,     Form("bh2dde[%d]/D", NumOfSegBH2));
#if HodoHitPos
  tree->Branch("bh2hitpos",     event.bh2hitpos,     Form("bh2hitpos[%d][%d]/D", NumOfSegBH2, MaxDepth));
#endif
  tree->Branch("bacmt",     event.bacmt,     Form("bacmt[%d][%d]/D", NumOfSegBAC, MaxDepth));
  tree->Branch("bacde",     event.bacde,     Form("bacde[%d]/D", NumOfSegBAC));
  tree->Branch("tofmt",     event.tofmt,     Form("tofmt[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofde",     event.tofde,     Form("tofde[%d]/D", NumOfSegTOF));
  tree->Branch("tofude",     event.tofude,     Form("tofude[%d]/D", NumOfSegTOF));
  tree->Branch("tofdde",     event.tofdde,     Form("tofdde[%d]/D", NumOfSegTOF));
  tree->Branch("tofctu",     event.tofctu,     Form("tofctu[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofctd",     event.tofctd,     Form("tofctd[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofcmt",     event.tofcmt,     Form("tofcmt[%d][%d]/D", NumOfSegTOF, MaxDepth));

  tree->Branch("t0",        event.t0,        Form("t0[%d][%d]/D",  NumOfSegBH2, MaxDepth));
  tree->Branch("ct0",       event.ct0,       Form("ct0[%d][%d]/D", NumOfSegBH2, MaxDepth));

  tree->Branch("Time0Seg", &event.Time0Seg,  "Time0Seg/D");
  tree->Branch("deTime0",  &event.deTime0,   "deTime0/D");
  tree->Branch("Time0",    &event.Time0,     "Time0/D");
  tree->Branch("CTime0",   &event.CTime0,    "CTime0/D");

  ////////////////////////////////////////////
  //Dst
  hodo = new TTree("hodo","Data Summary Table of Hodoscope");
  hodo->Branch("evnum",     &dst.evnum,     "evnum/I");
  hodo->Branch("spill",     &dst.spill,     "spill/I");
  hodo->Branch("trigpat",    dst.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  hodo->Branch("trigflag",   dst.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  hodo->Branch("nhBh2",     &dst.nhBh2,     "nhBh2/I");
  hodo->Branch("csBh2",      dst.csBh2,     "csBh2[nhBh2]/I");
  hodo->Branch("Bh2Seg",     dst.Bh2Seg,    "Bh2Seg[nhBh2]/D");
  hodo->Branch("tBh2",       dst.tBh2,      "tBh2[nhBh2]/D");
  hodo->Branch("t0Bh2",      dst.t0Bh2,     "t0Bh2[nhBh2]/D");
  hodo->Branch("dtBh2",      dst.dtBh2,     "dtBh2[nhBh2]/D");
  hodo->Branch("deBh2",      dst.deBh2,     "deBh2[nhBh2]/D");
#if HodoHitPos
  hodo->Branch("posBh2",     dst.posBh2,    "posBh2[nhBh2]/D");
#endif
  hodo->Branch("Time0Seg",  &dst.Time0Seg,  "Time0Seg/D");
  hodo->Branch("deTime0",   &dst.deTime0,   "deTime0/D");
  hodo->Branch("Time0",     &dst.Time0,     "Time0/D");
  hodo->Branch("CTime0",    &dst.CTime0,    "CTime0/D");

  hodo->Branch("nhBac",     &dst.nhBac,     "nhBac/I");
  hodo->Branch("BacSeg",     dst.BacSeg,    "BacSeg[nhBac]/D");
  hodo->Branch("tBac",       dst.tBac,      "tBac[nhBac]/D");
  hodo->Branch("deBac",      dst.deBac,     "deBac[nhBac]/D");

  hodo->Branch("nhTof",     &dst.nhTof,     "nhTof/I");
  hodo->Branch("csTof",      dst.csTof,     "csTof[nhTof]/I");
  hodo->Branch("TofSeg",     dst.TofSeg,    "TofSeg[nhTof]/D");
  hodo->Branch("tTof",       dst.tTof,      "tTof[nhTof]/D");
  hodo->Branch("dtTof",      dst.dtTof,     "dtTof[nhTof]/D");
  hodo->Branch("deTof",      dst.deTof,     "deTof[nhTof]/D");

  hodo->Branch("tofua", event.tofua, Form("tofua[%d]/D", NumOfSegTOF));
  hodo->Branch("tofut", event.tofut, Form("tofut[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("tofda", event.tofda, Form("tofda[%d]/D", NumOfSegTOF));
  hodo->Branch("tofdt", event.tofdt, Form("tofdt[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("utTofSeg", dst.utTofSeg,
               Form("utTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("dtTofSeg", dst.dtTofSeg,
               Form("dtTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("udeTofSeg",  dst.udeTofSeg,
               Form("udeTofSeg[%d]/D", NumOfSegTOF));
  hodo->Branch("ddeTofSeg",  dst.ddeTofSeg,
               Form("ddeTofSeg[%d]/D", NumOfSegTOF));
#endif // Normalized, Dst

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
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
