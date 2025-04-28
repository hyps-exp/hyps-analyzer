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
#include "HypsLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "DCGeomMan.hh"
#include "CFTPedCorMan.hh"

// #define TimeCut    1 // in cluster analysis
#define FHitBranch 0 // make FiberHit branches (becomes heavy)
#define HodoHitPos 0
#define ForPHC     0

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

  Double_t cftadc_h[NumOfPlaneCFT][NumOfSegCFT_PHI4];
  Double_t cftadc_l[NumOfPlaneCFT][NumOfSegCFT_PHI4];
  Double_t cfttdc[NumOfPlaneCFT][NumOfSegCFT_PHI4];

  Double_t cftadc_cor_h[NumOfPlaneCFT][NumOfSegCFT_PHI4];
  Double_t cftadc_cor_l[NumOfPlaneCFT][NumOfSegCFT_PHI4];



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

	Double_t width = tdc_l - tdc_t;
	HF2 (1000*(plane+1)+103, seg, width);


	if (tdc_l>MinTdcCFT && tdc_l < MaxTdcCFT) {
	  flag_tdc = true;
	  event.cfttdc[plane][seg] = (Double_t)tdc_l;
	}
      }
      if (flag_tdc) {
	HF1 (1000*(plane+1)+102, seg);
      }


      //ADC Hi
      for(Int_t m = 0; m<NhitAH; ++m){
	Int_t adcH = hit->GetAdcHigh();
	HF2 (1000*(plane+1)+104, seg, adcH);
	event.cftadc_h[plane][seg] = (Double_t)adcH;

	if (flag_tdc) {
	  HF2 (1000*(plane+1)+106, seg, adcH);
	}
      }

      //ADC Low
      for(Int_t m = 0; m<NhitAL; ++m){
	Int_t adcL = hit->GetAdcLow();
	HF2 (1000*(plane+1)+105, seg, adcL);
	event.cftadc_l[plane][seg] = (Double_t)adcL;
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
      event.cftadc_cor_h[plane][seg]  = adccorHi;
      event.cftadc_cor_l[plane][seg] = adccorLow;

      if (adccorHi>1500) {
	HF2 (1000*(plane+1)+201, seg, event.cfttdc[plane][seg]);
      }

      Int_t NhitT = hit->GetEntries(U);
      Bool_t flagTime = false;

      double time0 = -999;
      double ctime0 = -999;
      
      for(Int_t m = 0; m<NhitT; ++m){
	Double_t time = hit->GetTUp(m);
	Double_t ctime = hit->GetCTUp(m);		
	HF2 (1000*(plane+1)+200, seg, time);
	HF2 (1000*(plane+1)+210, seg, ctime);            			
	if (std::abs(time)<100)
	  flagTime = true;

	if (std::abs(ctime0) > std::abs(ctime)) {
	  time0 = time;
	  ctime0 = ctime;
	}
      }
      if (flagTime) {
	HF2 (1000*(plane+1)+206, seg, adccorHi);
#ifdef ForPHC
	int histId = ((plane+1)*1000+seg)*10+1;
	if (adccorHi>1000)
	  HF1( histId, time0);

	histId = ((plane+1)*1000+seg)*10+2;
	HF2( histId, 1./sqrt(adccorHi), time0);
#endif	
      }
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

      if ((cl->PlaneName()).Contains("PHI")) {
	//std::cout << "(x, y) = (" << cl->MeanX() << ", " << cl->MeanY() << ")" << std::endl;
	Double_t meanPhi = cl->MeanPhi();
	HF1(1000*(plane+1)+305, meanPhi);
      } else if ((cl->PlaneName()).Contains("UV")) {
	//std::cout << "(r, z0) = (" << cl->MeanR() << ", " << cl->MeanZ0() << ")" << std::endl;
	Double_t meanZ0 = cl->MeanZ0();
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

  rawData.DecodeHits("PiID");
  // PiID
  {
    const auto& cont = rawData.GetHodoRawHC("PiID");
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
	HF2 (9000+100, seg, tdc_l);  //TDC Nhits 
	HF2 (9000+101, seg, tdc_t);
	
	Int_t width = tdc_l - tdc_t;	  
	HF2 (9000+103, seg, width);


	if (tdc_l>MinTdcCFT && tdc_l < MaxTdcCFT) {
	  flag_tdc = true;
	  //event.cfttdc[plane][seg] = tdc_l;	
	}
      }
      if (flag_tdc) {
	HF1 (9000+102, seg);
      }

      
      //ADC Hi
      for(Int_t m = 0; m<NhitAH; ++m){
	Int_t adcH = hit->GetAdcHigh();	
	HF2 (9000+104, seg, adcH);
	//event.cftadc_h[plane][seg] = adcH;

	if (flag_tdc) {
	  HF2 (9000+106, seg, adcH);
	}
      }
	
      //ADC Low
      for(Int_t m = 0; m<NhitAL; ++m){	    
	Int_t adcL = hit->GetAdcLow();	
	HF2 (9000+105, seg, adcL);
	//event.cftadc_l[plane][seg] = adcL;	
      }
    }
  }

  hodoAna.DecodeHits<CFTFiberHit>("PiID");
  {
    const auto& U = HodoRawHit::kUp;
    Int_t nh=hodoAna.GetNHits("PiID");
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit<CFTFiberHit>("PiID", i);
      if(!hit) continue;
      Int_t seg   = hit->SegmentId();
      Int_t plane = hit->PlaneId();

      Double_t adccorHi  = hit->GetAdcCorHigh();
      Double_t adccorLow = hit->GetAdcCorLow();      
      Double_t mipHi  = hit->GetMipHigh();
      Double_t mipLow = hit->GetMipLow();      
      Double_t deLow  = hit->DeltaELowGain();
      

      HF2 (9000+204, seg, adccorHi);
      HF2 (9000+205, seg, adccorLow);            
      HF2 (9000+207, seg, mipHi);
      HF2 (9000+208, seg, mipLow);
      /*
      HF2 (1000*(plane+1)+209, seg, deLow);                  
      event.cftadc_cor_h[plane][seg]  = (int)adccorHi;
      event.cftadc_cor_l[plane][seg] = (int)adccorLow;

      if (adccorHi>1500) {
	HF2 (1000*(plane+1)+201, seg, event.cfttdc[plane][seg]);
      }
      */
      
      Int_t NhitT = hit->GetEntries(U);
      bool flagTime = false;

      double time0 = -999;
      double ctime0 = -999;

      for(Int_t m = 0; m<NhitT; ++m){
	Double_t time = hit->GetTUp(m);
	Double_t ctime = hit->GetCTUp(m);	

	HF2 (9000+200, seg, time);            	
	if (std::abs(time)<100)
	  flagTime = true;

	if (std::abs(ctime0) > std::abs(ctime)) {
	  time0 = time;
	  ctime0 = ctime;
	}
      }

      if (flagTime) {
	HF2 (9000+206, seg, adccorHi);
      }

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
    TString title210("");        
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
      title210 = Form("CFT UV %d : CTime vs seg", layer);            
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
      title210 = Form("CFT Phi %d : CTime vs seg", layer);            
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
    HB2( 1000*(i+1)+210, title200, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-500,500);    
    HB1( 1000*(i+1)+300, title300, 20, 0, 20);
    HB1( 1000*(i+1)+301, title301, 20, 0, 20);
    HB2( 1000*(i+1)+302, title302, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-500,500);
    HB2( 1000*(i+1)+303, title303, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-1,9);
    HB2( 1000*(i+1)+304, title304, NumOfSegCFT[i], 0, NumOfSegCFT[i], 1000,-1,9);
    HB1( 1000*(i+1)+305, title305, 500, -50, 450);
  }


#ifdef ForPHC
  // Consume too many memories 
  // use only when the each channnel study is necessary 
  for (int l=0; l<NumOfPlaneCFT; l++) {
    int NumOfSeg = NumOfSegCFT[l];
    for (int seg=0; seg<NumOfSeg; seg++ ) {
      int hid = ((l+1)*1000+seg)*10;

      TString title1 = Form("Time (ADCcor>1000) %d-%d", l, seg);
      HB1(hid+1, title1, 100, -50, 50);
      TString title2 = Form("Time : 1/sqrt(ADCcor) %d-%d", l, seg);
      HB2(hid+2, title2, 100, 0, 0.15, 100, -50, 50);
    }
  }
#endif


  // PiID
  {
    //ADC

    TString title100("PiID : Tdc(Leading) vs seg");
    TString title101("PiID : Tdc(Trailing) vs seg");
    TString title102("PiID : Hit pattern");    
    TString title103("PiID : TOT vs seg");    
    TString title104("PiID : Adc(High) vs seg");
    TString title105("PiID : Adc(Low) vs seg");
    TString title106("PiID : Adc(High) vs seg (w/ TDC)");      
    TString title200("PiID : Time vs seg");
    TString title201("PiID : TDC (w/ Large ADC) vs seg");      
    TString title204("PiID : AdcCor(High) vs seg");
    TString title205("PiID : AdcCor(Low) vs seg");
    TString title206("PiID : AdcCor(High) vs seg (w/ TDC)");
    TString title207("PiID : Mip Calib(High) vs seg");
    TString title208("PiID : Mip Calib(Low) vs seg");
    TString title209("PiID : dE (Low) vs seg");

    HB2( 9000+100, title100, NumOfSegPiID, 0, NumOfSegPiID, 1024,0,1024);
    HB2( 9000+101, title101, NumOfSegPiID, 0, NumOfSegPiID, 1024,0,1024);
    HB1( 9000+102, title102, NumOfSegPiID, 0, NumOfSegPiID);   
    HB2( 9000+103, title103, NumOfSegPiID, 0, NumOfSegPiID, 1024,0,1024);        
    HB2( 9000+104, title104, NumOfSegPiID, 0, NumOfSegPiID, 1000,0,4000);
    HB2( 9000+105, title105, NumOfSegPiID, 0, NumOfSegPiID, 1000,0,4000);
    HB2( 9000+106, title106, NumOfSegPiID, 0, NumOfSegPiID, 1000,0,4000);    

    HB2( 9000+200, title200, NumOfSegPiID, 0, NumOfSegPiID, 1000,-500,500);
    HB2( 9000+201, title201, NumOfSegPiID, 0, NumOfSegPiID, 1024,0,1024);
    HB2( 9000+204, title204, NumOfSegPiID, 0, NumOfSegPiID, 1000,-500,3500);
    HB2( 9000+205, title205, NumOfSegPiID, 0, NumOfSegPiID, 1000,-500,3500);
    HB2( 9000+206, title206, NumOfSegPiID, 0, NumOfSegPiID, 1000,-500,3500);
    HB2( 9000+207, title207, NumOfSegPiID, 0, NumOfSegPiID, 1000,-10,90);
    HB2( 9000+208, title208, NumOfSegPiID, 0, NumOfSegPiID, 1000,-1,9);
    HB2( 9000+209, title209, NumOfSegPiID, 0, NumOfSegPiID, 1000,-1,9);
  }



  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spill",     &event.spill,     "spill/I");
  //Trig
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  tree->Branch("cftadc_h",   event.cftadc_h,  Form("cftadc_h[%d][%d]/D", NumOfPlaneCFT, NumOfSegCFT_PHI4));
  tree->Branch("cftadc_l",   event.cftadc_l,  Form("cftadc_l[%d][%d]/D", NumOfPlaneCFT, NumOfSegCFT_PHI4));
  tree->Branch("cfttdc",     event.cfttdc,    Form("cfttdc[%d][%d]/D", NumOfPlaneCFT, NumOfSegCFT_PHI4));

  tree->Branch("cftadc_cor_h",   event.cftadc_cor_h,  Form("cftadc_cor_h[%d][%d]/D", NumOfPlaneCFT, NumOfSegCFT_PHI4));
  tree->Branch("cftadc_cor_l",   event.cftadc_cor_l,  Form("cftadc_cor_l[%d][%d]/D", NumOfPlaneCFT, NumOfSegCFT_PHI4));

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
