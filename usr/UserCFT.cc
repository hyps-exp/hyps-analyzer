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
#include "CFTPedCorMan.hh"

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

  hodoAna.DecodeHits<FiberHit>("CFT");
  {
    Int_t nh=hodoAna.GetNHits("CFT");

    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit<FiberHit>("CFT", i);
      if(!hit) continue;
      Int_t seg   = hit->SegmentId();
      Int_t plane = hit->PlaneId();


      Double_t adccorHi  = hit->GetAdcCorHigh();
      Double_t adccorLow = hit->GetAdcCorLow();      

      HF2 (1000*(plane+1)+204, seg, adccorHi);
      HF2 (1000*(plane+1)+205, seg, adccorLow);            
      event.cftadc_cor_h[plane][seg]  = (int)adccorHi;
      event.cftadc_cor_l[plane][seg] = (int)adccorLow;
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
  HB1(10, "Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig));
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000);
  }

  for(Int_t i=0; i<NumOfPlaneCFT; i++){
    //ADC

    TString title1("");
    TString title2("");
    TString title3("");    
    TString title4("");    
    TString title5("");
    TString title6("");
    TString title7("");    
    TString title15("");
    TString title16("");
    TString title17("");    
    if(i%2 == 0){// spiral layer
      Int_t layer = (Int_t)i/2 +1;
      title1  = Form("CFT UV %d : Tdc(Leading) vs seg", layer);
      title2  = Form("CFT UV %d : Tdc(Trailing) vs seg", layer);
      title3  = Form("CFT UV %d : Hit pattern", layer);      
      title4  = Form("CFT UV %d : TOT vs seg", layer);            
      title5  = Form("CFT UV %d : Adc(High) vs seg", layer);
      title6  = Form("CFT UV %d : Adc(Low) vs seg", layer);
      title7  = Form("CFT UV %d : Adc(High) vs seg (w/ TDC)", layer);      
      title15 = Form("CFT UV %d : AdcCor(High) vs seg", layer);
      title16 = Form("CFT UV %d : AdcCor(Low) vs seg", layer);
      title17 = Form("CFT UV %d : AdcCor(High) vs seg (w/ TDC)", layer);      
    }else if(i%2 == 1){// straight layer
      Int_t layer = (Int_t)i/2 +1;
      title1  = Form("CFT Phi %d : Tdc(Leading) vs seg", layer);
      title2  = Form("CFT Phi %d : Tdc(Trailing) vs seg", layer);
      title3  = Form("CFT Phi %d : Hit pattern", layer);      
      title4  = Form("CFT Phi %d : TOT vs seg", layer);            
      title5  = Form("CFT Phi %d : Adc(High) vs seg", layer);
      title6  = Form("CFT Phi %d : Adc(Low) vs seg", layer);
      title7  = Form("CFT Phi %d : Adc(High) vs seg (w/ TDC)", layer);
      title15 = Form("CFT Phi %d : AdcCor(High) vs seg", layer);
      title16 = Form("CFT Phi %d : AdcCor(Low) vs seg", layer);
      title17 = Form("CFT Phi %d : AdcCor(High) vs seg (w/ TDC)", layer);
    }
    HB2( 1000*(i+1)+100, title1, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1024,0,1024);
    HB2( 1000*(i+1)+101, title2, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1024,0,1024);
    HB1( 1000*(i+1)+102, title3, NumOfSegCFT[i], 0, NumOfSegCFT[i]);   
    HB2( 1000*(i+1)+103, title4, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1024,0,1024);        
    HB2( 1000*(i+1)+104, title5, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,0,4000);
    HB2( 1000*(i+1)+105, title6, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,0,4000);
    HB2( 1000*(i+1)+106, title7, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,0,4000);    
    HB2( 1000*(i+1)+204, title15, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-500,3500);
    HB2( 1000*(i+1)+205, title16, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-500,3500);
    HB2( 1000*(i+1)+206, title17, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-500,3500);    
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
     InitializeParameter<CFTPedCorMan>("CFTPED") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
