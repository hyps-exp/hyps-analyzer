// -*- C++ -*-
//-----------------------------------------------------------------
// other : Akao
// date: 2026-10-04
// note: SP8 beam test
//-----------------------------------------------------------------
#include "VEvent.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

#include <UnpackerManager.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "DCGeomMan.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"

namespace
{
using namespace root;
const auto qnan       = TMath::QuietNaN();
auto&      gUnpacker  = hddaq::unpacker::GUnpacker::get_instance();
auto&      gRM        = RMAnalyzer::GetInstance();

// CFT planes
const Int_t PlanePhi3  = 5; // kCFTPlane::PHI3
const Int_t PlanePhi4  = 7; // kCFTPlane::PHI4
const Int_t SegPhi3Min = 480, SegPhi3Max = 511;
const Int_t SegPhi4Min = 544, SegPhi4Max = 575;
const Int_t NumOfSegRayrawPerPlane = 32;

//_____________________________________________________________________________
// histogram numbering (phi3 -> 1-32, phi4 -> 101-132)
Bool_t
RayrawGroupId(Int_t plane, Int_t seg, Int_t& gid)
{
  if(plane==PlanePhi3 && seg>=SegPhi3Min && seg<=SegPhi3Max){
    gid = (seg - SegPhi3Min) + 1;
    return true;
  }
  if(plane==PlanePhi4 && seg>=SegPhi4Min && seg<=SegPhi4Max){
    gid = (seg - SegPhi4Min) + 101;
    return true;
  }
  return false;
}
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  std::vector<Int_t>                 plane;
  std::vector<Int_t>                 seg;
  std::vector<std::vector<Double_t>> fadc;         // waveform sample value   (10 bit)
  std::vector<std::vector<Double_t>> crs_cnt;      // waveform sample time tag ( 13.33 ns/count)
  std::vector<std::vector<Double_t>> tdc_leading;  // TDC leading  ( 0.8333 ns/count)
  std::vector<std::vector<Double_t>> tdc_trailing; // TDC trailing ( 0.8333 ns/count)
  std::vector<Double_t>              adc_max;      // max(fadc), one value per hit (for PHC)
  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum = 0;
  plane.clear();
  seg.clear();
  fadc.clear();
  crs_cnt.clear();
  tdc_leading.clear();
  tdc_trailing.clear();
  adc_max.clear();
}

//_____________________________________________________________________________
namespace root
{
Event  event;
TH1   *h[MaxHist];
TTree *tree;
enum eDetHid
  {
    RayrawCFTHid = 150000,
  };
}

//_____________________________________________________________________________
Bool_t
ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingNormal()
{
  RawData rawData;
  rawData.DecodeHits("CFT");

  event.evnum = gUnpacker.get_event_number();
  HF1(1, 0);

  const auto& cont = rawData.GetHodoRawHC("CFT");
  for(const auto& hit : cont){
    Int_t plane = hit->PlaneId();
    Int_t seg   = hit->SegmentId();

    Int_t gid = 0;
    if(!RayrawGroupId(plane, seg, gid)) continue;

    const auto& fadc     = hit->GetArrayAdcHigh();     // waveform samples
    const auto& crs_cnt  = hit->GetArrayAdcLow();     
    const auto& leading  = hit->GetArrayTdcLeading();
    const auto& trailing = hit->GetArrayTdcTrailing();

    event.plane.push_back(plane);
    event.seg.push_back(seg);
    event.fadc.push_back(std::vector<Double_t>(fadc.begin(), fadc.end()));
    event.crs_cnt.push_back(std::vector<Double_t>(crs_cnt.begin(), crs_cnt.end()));
    event.tdc_leading.push_back(std::vector<Double_t>(leading.begin(), leading.end()));
    event.tdc_trailing.push_back(std::vector<Double_t>(trailing.begin(), trailing.end()));
    event.adc_max.push_back(fadc.empty() ? qnan : *std::max_element(fadc.begin(), fadc.end()));

    Int_t hid_wf      = RayrawCFTHid + gid*10 + 0; // waveform
    Int_t hid_tdc_l   = RayrawCFTHid + gid*10 + 1;
    Int_t hid_tdc_t   = RayrawCFTHid + gid*10 + 2;
    Int_t hid_nsample = RayrawCFTHid + gid*10 + 3;

    Int_t nsample = std::min(fadc.size(), crs_cnt.size());
    HF1(hid_nsample, nsample);
    if(!leading.empty()){ //  fill waveform w/TDC
      for(Int_t s=0; s<nsample; ++s){
        HF2(hid_wf, s, fadc[s]);
      }
    }
    for(const auto& tdc_l : leading)  HF1(hid_tdc_l, tdc_l);
    for(const auto& tdc_t : trailing) HF1(hid_tdc_t, tdc_t);
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  const Double_t MinSample   = 0.;
  const Double_t MaxSample   = 57.;
  const Int_t    NbinSample  = (Int_t)(MaxSample - MinSample);
  const Double_t MinFadc     = 500.;
  const Double_t MaxFadc     = 1000.;
  const Int_t    NbinFadc    = (Int_t)(MaxFadc - MinFadc);
  const Double_t MinTdc      = 2000.;
  const Double_t MaxTdc      = 4000.;
  const Int_t    NbinTdc     = (Int_t)(MaxTdc - MinTdc);
  const Int_t    NbinNsample = 600;
  const Double_t MinNsample  = 0.;
  const Double_t MaxNsample  = 600.;

  HB1( 1, "Status",  20,   0., 20.);

  for(Int_t i=0; i<NumOfSegRayrawPerPlane; ++i){
    // phi3 (gid = i+1), phi4 (gid = i+101)
    for(const auto& gid : {i+1, i+101}){
      Bool_t isPhi4 = (gid>=101);
      Int_t  plane  = isPhi4 ? PlanePhi4 : PlanePhi3;
      Int_t  seg    = isPhi4 ? (SegPhi4Min+i) : (SegPhi3Min+i);

      TString title_wf = Form("RAYRAW CFT plane%d seg%d - Waveform w/TDC", plane, seg);
      TString title_tl = Form("RAYRAW CFT plane%d seg%d - TDC Leading",  plane, seg);
      TString title_tt = Form("RAYRAW CFT plane%d seg%d - TDC Trailing", plane, seg);
      TString title_ns = Form("RAYRAW CFT plane%d seg%d - N fadc sample", plane, seg);

      HB2( RayrawCFTHid + gid*10 + 0, title_wf, NbinSample, MinSample, MaxSample, NbinFadc, MinFadc, MaxFadc );
      HB1( RayrawCFTHid + gid*10 + 1, title_tl, NbinTdc,    MinTdc,   MaxTdc );
      HB1( RayrawCFTHid + gid*10 + 2, title_tt, NbinTdc,    MinTdc,   MaxTdc );
      HB1( RayrawCFTHid + gid*10 + 3, title_ns, NbinNsample,MinNsample,MaxNsample );
    }
  }

  //Tree
  HBTree( "tree", "tree" );
  tree->Branch("evnum",        &event.evnum,  "evnum/I");
  tree->Branch("plane",        &event.plane);
  tree->Branch("seg",          &event.seg);
  tree->Branch("fadc",         &event.fadc);
  tree->Branch("crs_cnt",      &event.crs_cnt);
  tree->Branch("tdc_leading",  &event.tdc_leading);
  tree->Branch("tdc_trailing", &event.tdc_trailing);
  tree->Branch("adc_max",      &event.adc_max);

  HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")    &&
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
