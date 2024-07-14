// -*- C++ -*-

// File   : UserHypsRaw.cc
// Author : Rintaro Kurata

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

  Int_t    tofnhits;
  Int_t    tofhitpat[MaxHits];
  Double_t tofua[NumOfSegTOF];
  Double_t tofut[NumOfSegTOF][MaxDepth];
  Double_t tofda[NumOfSegTOF];
  Double_t tofdt[NumOfSegTOF][MaxDepth];

  ////////// Normalized
  Double_t tofmt[NumOfSegTOF][MaxDepth];
  Double_t tofde[NumOfSegTOF];
  Double_t tofude[NumOfSegTOF];
  Double_t tofdde[NumOfSegTOF];
  Double_t tofctu[NumOfSegTOF][MaxDepth];
  Double_t tofctd[NumOfSegTOF][MaxDepth];
  Double_t tofcmt[NumOfSegTOF][MaxDepth];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum      = 0;
  spill      = 0;
  tofnhits   = 0;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<MaxHits; ++it){
    tofhitpat[it] = -1;
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

#if 0 // Dst
//_____________________________________________________________________________
struct Dst
{
  Int_t evnum;
  Int_t spill;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

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
  evnum    = 0;
  spill    = 0;
  nhTof    = 0;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
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
#endif // Dst

//_____________________________________________________________________________
namespace root
{
Event  event;
// Dst    dst;
TH1*   h[MaxHist];
TTree* tree;
// TTree* hodo;
enum eDetHid {
  TOFHid    = 60000,
};
}

//_____________________________________________________________________________
Bool_t
ProcessingBegin()
{
  event.clear();
  // dst.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingNormal()
{
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
  // dst.evnum   = gRM.EventNumber();
  // dst.spill   = gRM.SpillNumber();

  HF1(1, 0);
  //**************************************************************************
  //****************** RawData

  // Trigger Flag
  rawData.DecodeHits("TFlag");
  std::bitset<NumOfSegTrig> trigger_flag;
  for(const auto& hit: rawData.GetHodoRawHC("TFlag")){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc();
    if(tdc > 0){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      // dst.trigpat[trigger_flag.count()] = seg;
      // dst.trigflag[seg] = tdc;
      trigger_flag.set(seg);
      HF1(10, seg);
      HF1(10+seg, tdc);
    }
  }

  if(trigger_flag[trigger::kSpillOnEnd] || trigger_flag[trigger::kSpillOffEnd])
    return true;

  HF1(1, 1);

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
      HF1(TOFHid+1, seg-0.5);
      // Up
      Int_t Au = hit->GetAdcUp();
      HF1(TOFHid+100*(seg+1)+1, Au);
      event.tofua[seg] = Au;
      // dst.tofua[seg] = Au;
      Bool_t is_hit_u = false;
      Int_t m_u = 0;
      for(const auto& T: hit->GetArrayTdcUp()){
        HF1(TOFHid +100*(seg+1) +3, T);
        if(m_u < MaxDepth){
          event.tofut[seg][m_u] = T;
          // dst.tofut[seg][m_u] = T;
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
      // dst.tofda[seg] = Ad;
      Bool_t is_hit_d = false;
      Int_t m_d = 0;
      for(const auto& T: hit->GetArrayTdcDown()){
        HF1(TOFHid +100*(seg+1) +4, Double_t(T));
        if(m_d < MaxDepth){
          event.tofdt[seg][m_d] = T;
          // dst.tofdt[seg][m_d] = T;
          ++m_d;
        }
        if(MinTdcTOF < T && T < MaxTdcTOF)  is_hit_d = true;
      }
      if(is_hit_d) HF1(TOFHid+100*(seg+1)+6, Ad);
      else         HF1(TOFHid+100*(seg+1)+8, Ad);
      // Hitpat
      if(is_hit_u || is_hit_d){
        ++nh1; HF1(TOFHid+3, seg-0.5);
      }
      if(is_hit_u && is_hit_d){
        event.tofhitpat[tof_nhits++] = seg;
        ++nh2; HF1(TOFHid+5, seg-0.5);
      }
    }
    HF1(TOFHid+2, nh1); HF1(TOFHid+4, nh2);
    event.tofnhits = tof_nhits;
  }

#if 0 // Normalized, Dst
  //**************************************************************************
  //****************** NormalizedData

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
        HF1(TOFHid+11, seg-0.5);
        Double_t au  = hit->GetAUp(),   ad  = hit->GetADown();
        Double_t tu  = hit->GetTUp(),   td  = hit->GetTDown();
        Double_t ctu = hit->GetCTUp(),  ctd = hit->GetCTDown();
        Double_t mt  = hit->MeanTime(), cmt = hit->CMeanTime();
        Double_t de  = hit->DeltaE();
        Double_t ude  = hit->UDeltaE();
        Double_t dde  = hit->DDeltaE();
        event.tofmt[seg][m]  = mt;
        event.tofde[seg]     = de;
        event.tofude[seg]    = ude;
        event.tofdde[seg]    = dde;
        event.tofctu[seg][m] = ctu;
        event.tofctd[seg][m] = ctd;
        event.tofcmt[seg][m] = cmt;
        HF1(TOFHid+100*(seg+1)+11, tu);
	HF1(TOFHid+100*(seg+1)+12, td);
        HF1(TOFHid+100*(seg+1)+13, mt);
        HF1(TOFHid+100*(seg+1)+17, ctu);
	HF1(TOFHid+100*(seg+1)+18, ctd);
        HF1(TOFHid+100*(seg+1)+19, cmt);
	HF1(TOFHid+100*(seg+1)+20, ctu-ctd);
        // HF2(TOFHid+100*(seg+1)+21, tu, au);
	// HF2(TOFHid+100*(seg+1)+22, td, ad);
        // HF2(TOFHid+100*(seg+1)+23, ctu, au);
	// HF2(TOFHid+100*(seg+1)+24, ctd, ad);
        HF1(TOFHid+12, cmt);

        dst.utTofSeg[seg][m] = tu;
        dst.dtTofSeg[seg][m] = td;
        dst.udeTofSeg[seg]   = au;
        dst.ddeTofSeg[seg]   = ad;

        if(m == 0){
          HF1(TOFHid+100*(seg+1)+14, au);
	  HF1(TOFHid+100*(seg+1)+15, ad);
          HF1(TOFHid+100*(seg+1)+16, de);
	  HF1(TOFHid+13, de);
        }

        if(de>0.5){
          HF1(TOFHid+15, seg-0.5);
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
          HF2(TOFHid+21, seg1-0.5, seg2-0.5);
          HF2(TOFHid+22, ct1, ct2);
          HF1(TOFHid+23, ct2-ct1);
          if(std::abs(ct2-ct1)<3.0){
            HF2(TOFHid+24, seg1-0.5, seg2-0.5);
          }
        }//for(i2)
      }//for(i1)
    }

    Int_t nc = hodoAna.GetNClusters("TOF");
    HF1(TOFHid+30, Double_t(nc));
    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster("TOF", i);
      if(!cl) continue;
      Int_t cs     = cl->ClusterSize();
      Double_t ms  = cl->MeanSeg();
      Double_t cmt = cl->CMeanTime();
      Double_t de  = cl->DeltaE();
      HF1(TOFHid+31, Double_t(cs));
      HF1(TOFHid+32, ms-0.5);
      HF1(TOFHid+33, cmt); HF1(TOFHid+34, de);
    }
  }

  ////////// Dst
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
#endif // Normalized, Dst

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
    HB1(TOFHid +100*(i+1) +1, title1, NbinAdc,   MinAdc,   MaxAdc);
    HB1(TOFHid +100*(i+1) +2, title2, NbinAdc,   MinAdc,   MaxAdc);
    HB1(TOFHid +100*(i+1) +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(TOFHid +100*(i+1) +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(TOFHid +100*(i+1) +5, title5, NbinAdc,   MinAdc,   MaxAdc);
    HB1(TOFHid +100*(i+1) +6, title6, NbinAdc,   MinAdc,   MaxAdc);
    HB1(TOFHid +100*(i+1) +7, title7, NbinAdc,   MinAdc,   MaxAdc);
    HB1(TOFHid +100*(i+1) +8, title8, NbinAdc,   MinAdc,   MaxAdc);
  }

#if 0 //Normalized
  HB1(TOFHid +10, "#Hits Tof[Hodo]",      NumOfSegTOF+1, 0.,   Double_t(NumOfSegTOF+1));
  HB1(TOFHid +11, "Hitpat Tof[Hodo]",     NumOfSegTOF,   0.,   Double_t(NumOfSegTOF));
  HB1(TOFHid +12, "CMeanTime Tof",        500,           -5.,  45.);
  HB1(TOFHid +13, "dE Tof",               200,           -0.5, 4.5);
  HB1(TOFHid +14, "#Hits Tof[HodoGood]",  NumOfSegTOF+1, 0.,   Double_t(NumOfSegTOF+1));
  HB1(TOFHid +15, "Hitpat Tof[HodoGood]", NumOfSegTOF,   0.,   Double_t(NumOfSegTOF));

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
    HB1(TOFHid +100*(i+1) +11, title11, 500, -5.,   45.);
    HB1(TOFHid +100*(i+1) +12, title12, 500, -5.,   45.);
    HB1(TOFHid +100*(i+1) +13, title13, 500, -5.,   45.);
    HB1(TOFHid +100*(i+1) +14, title14, 200, -0.5,  4.5);
    HB1(TOFHid +100*(i+1) +15, title15, 200, -0.5,  4.5);
    HB1(TOFHid +100*(i+1) +16, title16, 200, -0.5,  4.5);
    HB1(TOFHid +100*(i+1) +17, title17, 500, -5.,   45.);
    HB1(TOFHid +100*(i+1) +18, title18, 500, -5.,   45.);
    HB1(TOFHid +100*(i+1) +19, title19, 500, -5.,   45.);
    HB1(TOFHid +100*(i+1) +20, title20, 200, -10.0, 10.0);
    HB1(TOFHid +100*(i+1) +21, title21, 200, -0.5,  4.5);
  }

  HB2(TOFHid +21, "TofHitPat%TofHitPat[HodoGood]", NumOfSegTOF,   0., Double_t(NumOfSegTOF),
      NumOfSegTOF,   0., Double_t(NumOfSegTOF));
  HB2(TOFHid +22, "CMeanTimeTof%CMeanTimeTof[HodoGood]",
      120, 10., 40., 120, 10., 40.);
  HB1(TOFHid +23, "TDiff Tof[HodoGood]", 200, -10., 10.);
  HB2(TOFHid +24, "TofHitPat%TofHitPat[HodoGood2]", NumOfSegTOF,   0., Double_t(NumOfSegTOF),
      NumOfSegTOF,   0., Double_t(NumOfSegTOF));

  HB1(TOFHid +30, "#Clusters Tof",         NumOfSegTOF+1, 0.,   Double_t(NumOfSegTOF+1));
  HB1(TOFHid +31, "ClusterSize Tof",       5,             0.,   5.);
  HB1(TOFHid +32, "HitPat Cluster Tof",    2*NumOfSegTOF, 0.,   Double_t(NumOfSegTOF));
  HB1(TOFHid +33, "CMeamTime Cluster Tof", 500,           -5.,  45.);
  HB1(TOFHid +34, "DeltaE Cluster Tof",    100,           -0.5, 4.5);
#endif // Normalized

  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,      "evnum/I");
  tree->Branch("spill",     &event.spill,      "spill/I");
  //Trig
  tree->Branch("trigpat",    event.trigpat,    Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,   Form("trigflag[%d]/I", NumOfSegTrig));

  //TOF
  tree->Branch("tofnhits",  &event.tofnhits,   "tofnhits/I");
  tree->Branch("tofhitpat",  event.tofhitpat,  Form("tofhitpat[%d]/I", NumOfSegTOF));
  tree->Branch("tofua",      event.tofua,      Form("tofua[%d]/D", NumOfSegTOF));
  tree->Branch("tofut",      event.tofut,      Form("tofut[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofda",      event.tofda,      Form("tofda[%d]/D", NumOfSegTOF));
  tree->Branch("tofdt",      event.tofdt,      Form("tofdt[%d][%d]/D", NumOfSegTOF, MaxDepth));

#if 0 // Normalized, Dst
  //Normalized data
  tree->Branch("tofmt",      event.tofmt,      Form("tofmt[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofde",      event.tofde,      Form("tofde[%d]/D", NumOfSegTOF));
  tree->Branch("tofude",     event.tofude,     Form("tofude[%d]/D", NumOfSegTOF));
  tree->Branch("tofdde",     event.tofdde,     Form("tofdde[%d]/D", NumOfSegTOF));
  tree->Branch("tofctu",     event.tofctu,     Form("tofctu[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofctd",     event.tofctd,     Form("tofctd[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofcmt",     event.tofcmt,     Form("tofcmt[%d][%d]/D", NumOfSegTOF, MaxDepth));

  ////////////////////////////////////////////
  //Dst
  hodo = new TTree("hodo","Data Summary Table of Hodoscope");
  hodo->Branch("evnum",     &dst.evnum,        "evnum/I");
  hodo->Branch("spill",     &dst.spill,        "spill/I");
  hodo->Branch("trigpat",    dst.trigpat,      Form("trigpat[%d]/I", NumOfSegTrig));
  hodo->Branch("trigflag",   dst.trigflag,     Form("trigflag[%d]/I", NumOfSegTrig));

  hodo->Branch("nhTof",     &dst.nhTof,        "nhTof/I");
  hodo->Branch("csTof",      dst.csTof,        "csTof[nhTof]/I");
  hodo->Branch("TofSeg",     dst.TofSeg,       "TofSeg[nhTof]/D");
  hodo->Branch("tTof",       dst.tTof,         "tTof[nhTof]/D");
  hodo->Branch("dtTof",      dst.dtTof,        "dtTof[nhTof]/D");
  hodo->Branch("deTof",      dst.deTof,        "deTof[nhTof]/D");

  hodo->Branch("tofua",      event.tofua,      Form("tofua[%d]/D", NumOfSegTOF));
  hodo->Branch("tofut",      event.tofut,      Form("tofut[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("tofda",      event.tofda,      Form("tofda[%d]/D", NumOfSegTOF));
  hodo->Branch("tofdt",      event.tofdt,      Form("tofdt[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("utTofSeg",   dst.utTofSeg,     Form("utTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("dtTofSeg",   dst.dtTofSeg,     Form("dtTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("udeTofSeg",  dst.udeTofSeg,    Form("udeTofSeg[%d]/D", NumOfSegTOF));
  hodo->Branch("ddeTofSeg",  dst.ddeTofSeg,    Form("ddeTofSeg[%d]/D", NumOfSegTOF));
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
