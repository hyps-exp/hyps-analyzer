// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <TSystem.h>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DatabasePDG.hh"
#include "DetectorID.hh"
#include "EventDisplay.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "CFTFiberHit.hh"
#include "CFTFiberCluster.hh"
#include "FuncName.hh"
#include "HodoAnalyzer.hh"
#include "HodoHit.hh"
#include "HodoWaveformHit.hh"
#include "HypsLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UnpackerManager.hh"
#include "BH2Filter.hh"
#include "TemplateFitMan.hh"
#include "BGOCalibMan.hh"
#include "CFTPedCorMan.hh"
#include "CFTPosParamMan.hh"
#include "CFTLocalTrack.hh"
#include "CFTParticle.hh"
#include "CATCHPidMan.hh"
#include "TAGPLMatch.hh"
#include "EventSelectMan.hh"

#define SAVEPDF 0
#define DEBUG 0
#define UseTOF 1
#define TAG_PL_FADC 1

namespace
{
using namespace root;
const auto& gGeom   = DCGeomMan::GetInstance();
auto&       gEvDisp = EventDisplay::GetInstance();
const auto& gUser   = UserParamMan::GetInstance();
auto&       gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto&       gTAGPLMth = TAGPLMatch::GetInstance();    
const auto& gEventMan = EventSelectMan::GetInstance();
const Double_t PionMass   = pdg::PionMass();
const Double_t KaonMass   = pdg::KaonMass();
const Double_t ProtonMass = pdg::ProtonMass();
}

//_____________________________________________________________________________
Bool_t
ProcessingBegin()
{
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingNormal()
{
  static const auto MaxMultiHitBcOut = gUser.GetParameter("MaxMultiHitBcOut");
  static const auto MaxMultiHitSdcIn = gUser.GetParameter("MaxMultiHitSdcIn");
  static const auto MaxMultiHitSdcOut = gUser.GetParameter("MaxMultiHitSdcOut");

  // static const auto MinTimeAFT = gUser.GetParameter("TimeAFT", 0);
  // static const auto MaxTimeAFT = gUser.GetParameter("TimeAFT", 1);
  static const auto MinTdcBH2 = gUser.GetParameter("TdcBH2", 0);
  static const auto MaxTdcBH2 = gUser.GetParameter("TdcBH2", 1);
  static const auto MinTdcBH1 = gUser.GetParameter("TdcBH1", 0);
  static const auto MaxTdcBH1 = gUser.GetParameter("TdcBH1", 1);
  static const auto MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const auto MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
  static const auto MinTdcTOF  = gUser.GetParameter("TdcTOF", 0);
  static const auto MaxTdcTOF  = gUser.GetParameter("TdcTOF", 1);
  static const auto MinTdcWC = gUser.GetParameter("TdcWC", 0);
  static const auto MaxTdcWC = gUser.GetParameter("TdcWC", 1);

  static const auto MinTimeCFT = gUser.GetParameter("TimeCFT", 0);
  static const auto MaxTimeCFT = gUser.GetParameter("TimeCFT", 1);
  static const auto MinAdcCFT = gUser.GetParameter("AdcCFT", 0);
  static const auto MaxAdcCFT = gUser.GetParameter("AdcCFT", 1);

  static const auto MinTotSDC0 = gUser.GetParameter("MinTotSDC0");
  static const auto MinTotSDC1 = gUser.GetParameter("MinTotSDC1");
  static const auto MinTotSDC2 = gUser.GetParameter("MinTotSDC2");
  static const auto MinTotSDC3 = gUser.GetParameter("MinTotSDC3");

  // static const auto StopTimeDiffSdcOut = gUser.GetParameter("StopTimeDiffSdcOut");
  // static const auto MinStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 0);
  // static const auto MaxStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 1);
  // static const auto MinTotBcOut = gUser.GetParameter("MinTotBcOut");
  // static const auto MinTotSDC3 = gUser.GetParameter("MinTotSDC3");
  // static const auto MinTotSDC4 = gUser.GetParameter("MinTotSDC4");

  // static const Int_t IdBH2 = gGeom.GetDetectorId("BH2");
  static const Int_t IdTOF = gGeom.GetDetectorId("TOF-X");
  //static const Int_t IdWC = gGeom.GetDetectorId("WC");

  if (gEventMan.RunNumCheck(gUnpacker.get_run_number()))
    if (! gEventMan.IsGood(gUnpacker.get_event_number()))
      return true;

  static TString evinfo;
  evinfo = Form("Run# %5d%4sEvent# %6d",
                gUnpacker.get_run_number(), "",
                gUnpacker.get_event_number());
  hddaq::cout << "\033c" << TString('=', 80) << std::endl
              << "[Info] " << evinfo << std::endl;
  gEvDisp.DrawRunEvent(0.04, 0.5, evinfo);

  RawData rawData;
  rawData.DecodeHits();
  HodoAnalyzer hodoAna(rawData);
  DCAnalyzer DCAna(rawData);

  //________________________________________________________
  //___ TrigRawHit
  /*
  std::bitset<NumOfSegTrig> trigger_flag;
  for(auto& hit: rawData.GetHodoRawHitContainer("TFlag")){
    if(hit->GetTdc(0) > 0) trigger_flag.set(hit->SegmentId());
  }
  if(trigger_flag[trigger::kSpillOnEnd] || trigger_flag[trigger::kSpillOffEnd])
    return true;
  // if(!trigger_flag[trigger::kTrigBPS]) return true;
  hddaq::cout << "[Info] TrigPat = " << trigger_flag << std::endl;
  */
  //________________________________________________________
  //___ AFT
  /*
  hodoAna.DecodeHits<FiberHit>("AFT");
  for(Int_t i=0, n=hodoAna.GetNHits("AFT"); i<n; ++i){
    const auto& hit = hodoAna.GetHit<FiberHit>("AFT", i);
    Int_t plane = hit->PlaneId();
    Int_t seg = hit->SegmentId();
    // Int_t m = hit->GetEntries();
    // for(Int_t j=0; j<m; ++j){
    //   auto mt = hit->MeanTime(j);
    //   if( MinTimeAFT < mt && mt < MaxTimeAFT ){
    // 	auto de_high = hit->DeltaEHighGain();
    // 	gEvDisp.FillAFT(plane, seg, de_high);
    // 	break;
    //   }
    // }
    auto de_high = hit->DeltaEHighGain();
    gEvDisp.FillAFT(plane, seg, de_high);
  }
  */
  //________________________________________________________
  //___ BH2RawHit
  /*
  std::vector<Int_t> BH2SegCont;
  for(const auto& hit: rawData.GetHodoRawHitContainer("BH2")){
    Int_t seg = hit->SegmentId();
    Bool_t is_hit_u = false;
    Bool_t is_hit_d = false;
    for(const auto& tdc: hit->GetArrayTdcLeading(0)){
      if(MinTdcBH2 < tdc && tdc < MaxTdcBH2){
        is_hit_u = true;
      }
      // gEvDisp.FillBH2(seg, Tu);
    }
    for(const auto& tdc: hit->GetArrayTdcLeading(1)){
      if(MinTdcBH2 < tdc && tdc < MaxTdcBH2){
        is_hit_d = true;
      }
      // gEvDisp.FillBH2(seg, Td);
    }
    if(is_hit_u && is_hit_d){
      hddaq::cout << "[Info] Bh2Seg = " << seg << std::endl;
      BH2SegCont.push_back(seg);
    }
  }
  */
  //________________________________________________________
  //___ BH1RawHit
  /*
  std::vector<Int_t> BH1SegCont;
  for(const auto& hit: rawData.GetHodoRawHitContainer("BH1")){
    Int_t seg = hit->SegmentId();
    Bool_t is_hit_u = false;
    Bool_t is_hit_d = false;
    for(const auto& tdc: hit->GetArrayTdcLeading(0)){
      if(MinTdcBH1 < tdc && tdc < MaxTdcBH1){
        is_hit_u = true;
      }
      // gEvDisp.FillBH1(seg, Tu);
    }
    for(const auto& tdc: hit->GetArrayTdcLeading(1)){
      if(MinTdcBH1 < tdc && tdc < MaxTdcBH1){
        is_hit_d = true;
      }
      // gEvDisp.FillBH1(seg, Td);
    }
    if(is_hit_u && is_hit_d){
      hddaq::cout << "[Info] Bh1Seg = " << seg << std::endl;
      BH1SegCont.push_back(seg);
    }
  }
  */
  //________________________________________________________
  //___ BH2HodoCluster
  /*
  hodoAna.DecodeHits("BH2");
  const auto Time0Cl = hodoAna.GetTime0BH2Cluster();
  Double_t ctime0 = 0.;
  if(!Time0Cl){
    hddaq::cout << "[Warning] Time0Cl is null!" << std::endl;
    // gEvDisp.GetCommand();
    // return true;
  } else {
    ctime0 = Time0Cl->CTime0();
  }
  // Double_t t0seg = Time0Cl->MeanSeg();
  */
  //________________________________________________________
  //___ BH1HodoCluster
  //hodoAna.DecodeHits("BH1");
  //Double_t btof = hodoAna.Btof0();

  //////////////T0 Analysis
  hodoAna.DecodeHits("T0");
  Int_t nhT0 = hodoAna.GetNHits("T0");
  Double_t time0 = -999.;
  //////////////T0 Analysis
  Double_t min_time = -999.;
  for(Int_t i=0; i<nhT0; ++i){
    const auto& hit = hodoAna.GetHit("T0", i);
    if(!hit) continue;
    Double_t multi = hit->GetEntries();
    for (Int_t im = 0; im <multi; ++im) {    
      Double_t mt  = hit->MeanTime(im);
      Double_t cmt  = hit->CMeanTime(im);      
      if(std::abs(mt)<std::abs(min_time)){
	time0    = cmt;
	min_time = mt;
      }
    }
  }

  
  //________________________________________________________
  //___ TOFRawHit
  std::vector<Int_t> TOFSegCont;
  for(const auto& hit: rawData.GetHodoRawHitContainer("TOF")){
    Int_t seg = hit->SegmentId();
    Bool_t is_hit_u = false;
    Bool_t is_hit_d = false;
    for(const auto& tdc: hit->GetArrayTdcLeading(0)){
      if(MinTdcTOF < tdc && tdc < MaxTdcTOF){
        is_hit_u = true;
      }
      // gEvDisp.FillTOF(seg, Tu);
    }
    for(const auto& tdc: hit->GetArrayTdcLeading(1)){
      if(MinTdcTOF < tdc && tdc < MaxTdcTOF){
        is_hit_d = true;
      }
      // gEvDisp.FillTOF(seg, Td);
    }
    if(is_hit_u || is_hit_d){
      gEvDisp.DrawHitHodoscope(IdTOF, seg, is_hit_u, is_hit_d);
    }
    if(is_hit_u && is_hit_d){
      hddaq::cout << "[Info] TofSeg = " << seg << std::endl;
      TOFSegCont.push_back(seg);
    }
  }

  //________________________________________________________
  //___ TOFHodoHit
  hodoAna.DecodeHits("TOF");
  //const auto& TOFCont = hodoAna.GetHitContainer("TOF");
  //if(TOFCont.empty()){
  //hddaq::cout << "[Warning] TOFCont is empty!" << std::endl;
  //gEvDisp.GetCommand();
  // return true;
  //}
  const auto& TOFCont = hodoAna.GetClusterContainer("TOF");
  Int_t nhTof = hodoAna.GetNClusters("TOF");
  {
    Double_t TofSeg[MaxHits];
    for (Int_t i=0; i<MaxHits; ++i)
      TofSeg[i] = -1.;
    
    for(Int_t i=0; i<nhTof; ++i){
      const auto& hit = hodoAna.GetCluster("TOF", i);
      Double_t seg = hit->MeanSeg()+1;
      //Double_t cmt = hit->CMeanTime();
      //Double_t dt  = hit->TimeDiff();
      //Double_t de   = hit->DeltaE();
      TofSeg[i] = seg;
    }

    Bool_t tofsingleflag=nhTof==1 && (TofSeg[0]<37 && TofSeg[0]>12);
    Bool_t tofdoubleflag=nhTof==2 && (TofSeg[0]<37 && TofSeg[0]>12 && TofSeg[1]<37 && TofSeg[1]>12);
      
    //if(!(tofsingleflag || tofdoubleflag)) return true;
    //if(!tofsingleflag) return true;
  }

  //Tag-Hodoana
  std::vector<int> PLCand;
  {  
    Int_t nc = 0;
#ifdef TAG_PL_FADC
    hodoAna.DecodeHits<HodoWaveformHit>("TAG-PL");    
#elif  
    hodoAna.DecodeHits("TAG-PL");
#endif
  
    Int_t nh=hodoAna.GetNHits("TAG-PL");
    Int_t nseg_goodtime=0;
    for(Int_t i=0;i<nh;++i){
#ifdef TAG_PL_FADC
      const auto& hit = hodoAna.GetHit<HodoWaveformHit>("TAG-PL", i);
#elif        
      const auto& hit = hodoAna.GetHit("TAG-PL",i);
#endif      
      if(!hit) continue;
      Int_t seg =hit->SegmentId();
      //std::cout<<"Tagger PL (seg) = ("<<seg<<std::endl;      
      bool is_hit_time =false;
      Int_t n_mhit =hit->GetEntries();
      for(Int_t m=0;m<n_mhit;++m){
	Double_t t=hit->GetTUp(m);
	//std::cout << t << "  ";
	//Double_t ct= hit->GetCTUp(m);
	if(fabs(t)<10.0) is_hit_time=true;
      }
      //std::cout << std::endl;
      if(is_hit_time){
	nseg_goodtime++;
	PLCand.push_back(seg);
	Double_t de = hit->DeltaE();
	gEvDisp.ShowHitTagger("TAG-PL", seg, de);	  	
      }

#ifdef TAG_PL_FADC
      Int_t ngr = hit->GetNGraph();
      if (ngr>0)
	nc++;
#endif    

    }

#ifdef TAG_PL_FADC
    gEvDisp.SetTagWaveformCanvas(nc);    
    nc = 1;
    for(Int_t i=0;i<nh;++i){
      const auto& hit = hodoAna.GetHit<HodoWaveformHit>("TAG-PL", i);
      if(!hit) continue;
      Int_t seg =hit->SegmentId();
      Int_t ngr = hit->GetNGraph();
      for (Int_t ig=0; ig<ngr; ig++) {
	TGraphErrors *gr = hit->GetTGraph(ig);
	gEvDisp.DrawTagWaveform(nc, ig, seg, gr);      
      }
      if (ngr>0)
	nc++;
    }
#endif
    
  }

  std::vector<double> SFFhit;
  std::vector<double> SFBhit;
  std::vector<double> SFFCand;
  std::vector<double> SFBCand;
  {  
    hodoAna.DecodeHits("TAG-SF");

    Int_t nh=hodoAna.GetNHits("TAG-SF");
    Int_t nseg_goodtime=0;
    for(Int_t i=0;i<nh;++i){
      const auto& hit =hodoAna.GetHit("TAG-SF",i);
      if(!hit) continue;
      Int_t seg =hit->SegmentId();
      Int_t plane =hit->PlaneId();
      TString planename=hit->PlaneName();
      //std::cout<<"(seg, plane) = ("<<seg<<", "<<plane<<")"<<std::endl;
      /*
	Double_t a =hit->GetAUp();
	
      */
      bool is_hit_time =false;
      Int_t n_mhit =hit->GetEntries();
      for(Int_t m=0;m<n_mhit;++m){
	Double_t t=hit->GetTUp(m);
	//std::cout << t << "  ";
	//Double_t ct= hit->GetCTUp(m);
	if(fabs(t)<10) is_hit_time=true;
      }
      //std::cout << std::endl;
      if(is_hit_time){
	//if(plane==0) SFFhit.push_back(seg);
	//if(plane==1) SFBhit.push_back(seg);
	Double_t de = hit->DeltaE();
	if (plane==0)
	  gEvDisp.ShowHitTagger("TAG-SFF", seg, de);
	else if (plane==1)
	  gEvDisp.ShowHitTagger("TAG-SFB", seg, de);	  
      }
    }
    
    hodoAna.TimeCut("TAG-SF",-10,10);
  
    Int_t nc=hodoAna.GetNClusters("TAG-SF");
    Int_t ncl1=0, ncl2=0;
    for(Int_t i=0;i<nc;++i){
      const auto& cl =hodoAna.GetCluster("TAG-SF",i);
      if(!cl) continue;
      Int_t plane=cl->PlaneId();
      if(plane==0) ncl1++;
      if(plane==1) ncl2++;
      Double_t ms =cl->MeanSeg();
      Double_t mt =cl->MeanTime();
      
      if(plane==0 && cl->ClusterSize()<4) SFFhit.push_back(ms);
      if(plane==1 && cl->ClusterSize()<4) SFBhit.push_back(ms);
      
      bool is_plmth=false;
      if(PLCand.size()==0) is_plmth=true;
      for(Int_t j=0;j<PLCand.size();j++){
	if(gTAGPLMth.Judge(ms,PLCand[j])) is_plmth=true;	  
      }
      if(is_plmth && cl->ClusterSize()<4){
	if(plane==0) SFFCand.push_back(ms);
	if(plane==1) SFBCand.push_back(ms);
      }
    }
  }
  
  std::vector<double> SFFCand_final;
  std::vector<double> SFBCand_final;
  double Egamma=0.0;

  const double eparf[3]={1.486,0.03312,-0.0001588};
  const double eparb[3]={1.49797,0.0327588,-0.000152769};
  const double offset_b=0.6421;

  if(SFFCand.size()>0 && SFBCand.size()>0){
    for(int i=0;i<SFFCand.size();i++){
      for(int j=0;j<SFBCand.size();j++){
	if(fabs(SFFCand[i]-SFBCand[j])<3){
	  SFFCand_final.push_back(SFFCand[i]);
	  SFBCand_final.push_back(SFBCand[j]);
	}
      }
    }
  }else if(SFFCand.size()==1 && PLCand.size()>0){
    Egamma=eparf[0]+eparf[1]*SFFCand[0]+eparf[2]*SFFCand[0]*SFFCand[0];
  }else if(SFBCand.size()==1 && PLCand.size()>0){
    double SFBpos=SFBCand[0]+offset_b;
    Egamma=eparb[0]+eparb[1]*SFBpos+eparb[2]*SFBpos*SFBpos;
  }
  
  if(SFFCand_final.size()==1){
    
    double egamf=eparf[0]+eparf[1]*SFFCand_final[0]+eparf[2]*SFFCand_final[0]*SFFCand_final[0];
    double SFBpos=SFBCand_final[0]+offset_b;
    double egamb=eparb[0]+eparb[1]*SFBpos+eparb[2]*SFBpos*SFBpos;

    Egamma=(egamf+egamb)*0.5;

  }

  if(Egamma>0){
  std::cout << "Egamma = " << Egamma << " GeV"<< std::endl;
  }

  //________________________________________________________
  //___ WCRawHit
  /*
  std::vector<Int_t> WCSegCont;
  for(const auto& hit: rawData.GetHodoRawHitContainer("WC")){
    Int_t seg = hit->SegmentId();
    Bool_t is_hit = false;
    for(const auto T: hit->GetArrayTdcLeading(HodoRawHit::kExtra)){
      if(MinTdcWC < T && T < MaxTdcWC){
        is_hit = true;
      }
    }
    if(is_hit){
      gEvDisp.DrawHitHodoscope(IdWC, seg, is_hit, is_hit);
      hddaq::cout << "[Info] WcSeg = " << seg << std::endl;
      WCSegCont.push_back(seg);
    }
  }
  */
#if 0
  static const Int_t IdSDC0 = gGeom.DetectorId("SDC0-X");
  static const Int_t IdSDC1 = gGeom.DetectorId("SDC1-X1");
  // static const Int_t IdSDC2 = gGeom.DetectorId("SDC2-X1");
  static const Int_t IdSDC3 = gGeom.DetectorId("SDC3-X1");
  static const Int_t IdSDC4 = gGeom.DetectorId("SDC4-X1");

  //________________________________________________________
  //___ SdcInRawHit
  for(const auto& hit: rawData.GetSdcInRawHC(IdSDC1)){
    Int_t wire = hit->WireId();
    for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
      Int_t tdc = hit->GetTdc(j);
      if(tdc>0) gEvDisp.FillSDC1(wire, tdc);
    }
  }
  for(const auto& hit: rawData.GetSdcInRawHC(IdSDC1+1)){
    Int_t wire = hit->WireId();
    for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
      Int_t tdc = hit->GetTdc(j);
      if(tdc>0) gEvDisp.FillSDC1p(wire, tdc);
    }
  }
  //________________________________________________________
  //___ SdcOutRawHit
  for(const auto& hit: rawData.GetSdcOutRawHC(IdSDC3-30)){
    Int_t wire = hit->WireId();
    for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
      Int_t tdc = hit->GetTdc(j);
      if(tdc>0) gEvDisp.FillSDC3_Leading(wire, tdc);
    }
    for(Int_t j=0, m=hit->GetTrailingSize(); j<m; ++j){
      Int_t tdc = hit->GetTrailing(j);
      if(tdc>0) gEvDisp.FillSDC3_Trailing(wire, tdc);
    }
  }
  for(const auto& hit: rawData.GetSdcOutRawHC(IdSDC3+1-30)){
    Int_t wire = hit->WireId();
    for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
      Int_t tdc = hit->GetTdc(j);
      if(tdc>0) gEvDisp.FillSDC3p_Leading(wire, tdc);
    }
    for(Int_t j=0, m=hit->GetTrailingSize(); j<m; ++j){
      Int_t tdc = hit->GetTrailing(j);
      if(tdc>0) gEvDisp.FillSDC3p_Trailing(wire, tdc);
    }
  }
  for(const auto& hit: rawData.GetSdcOutRawHC(IdSDC4-30)){
    Int_t wire = hit->WireId();
    for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
      Int_t tdc = hit->GetTdc(j);
      if(tdc>0) gEvDisp.FillSDC4_Leading(wire, tdc);
    }
    for(Int_t j=0, m=hit->GetTrailingSize(); j<m; ++j){
      Int_t tdc = hit->GetTrailing(j);
      if(tdc>0) gEvDisp.FillSDC4_Trailing(wire, tdc);
    }
  }
  for(const auto& hit: rawData.GetSdcOutRawHC(IdSDC4+1-30)){
    Int_t wire = hit->WireId();
    for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
      Int_t tdc = hit->GetTdc(j);
      if(tdc>0) gEvDisp.FillSDC4p_Leading(wire, tdc);
    }
    for(Int_t j=0, m=hit->GetTrailingSize(); j<m; ++j){
      Int_t tdc = hit->GetTrailing(j);
      if(tdc>0) gEvDisp.FillSDC4p_Trailing(wire, tdc);
    }
  }
  //________________________________________________________
  //___ BFTRawHit
  for(const auto& hit: rawData.GetHodoRawHitContainer("BFT")){
    Int_t plane = hit->PlaneId();
    Int_t seg = hit->SegmentId();
    for(Int_t j=0, m=hit->GetSizeTdcUp(); j<m; ++j){
      Int_t Tu = hit->GetTdcUp(j);
      // hddaq::cout << "[Info] BFTxSeg = " << seg << std::endl;
      if(Tu>0) gEvDisp.FillBFT(layer, seg, Tu);
    }
  }

#endif

  //________________________________________________________
  //___ SdcInDCHit
  rawData.TdcCutSDCIn();
  DCAna.DecodeSdcInHits();
  DCAna.TotCutSDC0(MinTotSDC0);
  DCAna.TotCutSDC1(MinTotSDC1);

  Double_t multi_SdcIn = 0.;
  for(Int_t plane=0; plane<NumOfLayersSdcIn; ++plane){
    const auto& cont = DCAna.GetSdcInHC(plane);
    Int_t n = cont.size();
    multi_SdcIn += n;
    if(n > MaxMultiHitSdcIn) continue;
    for(const auto& hit: cont){
      Int_t layer = hit->LayerId();
      Int_t wire = hit->GetWire();
      Int_t mhit = hit->GetEntries();

      Bool_t is_good = false;
      for(Int_t j=0; j<mhit; j++){
        if(hit->IsGood(j)){
          is_good = true;
	  Double_t dtime = hit->GetDriftTime(j);
	  //hddaq::cout << "SdcIn layer = " << layer << ", wire = " << wire
	  //<< ", dt = " << dtime << std::endl;	
          break;
        }
      }
      gEvDisp.DrawHitWire(layer, wire, is_good, is_good);
    }
  }
  multi_SdcIn /= (Double_t)NumOfLayersSdcIn;
  if(multi_SdcIn > MaxMultiHitSdcIn){
    hddaq::cout << "[Warning] SdcInHits exceed MaxMultiHit "
                << multi_SdcIn << "/" << MaxMultiHitSdcIn << std::endl;
    // gEvDisp.GetCommand();
    return true;
  }

  //________________________________________________________
  //___ SdcInTracking
  DCAna.TrackSearchSdcIn();
  Int_t ntSdcIn = DCAna.GetNtracksSdcIn();
  hddaq::cout << "[Info] ntSdcIn = " << ntSdcIn << std::endl;
  for(Int_t it=0; it<ntSdcIn; ++it){
    auto track = DCAna.GetTrackSdcIn(it);
    auto chisqr = track->GetChiSquare();
    hddaq::cout << "       " << it << "-th track, chi2 = "
                << chisqr << std::endl;
    // auto nh = track->GetNHit();
    // for(Int_t ih=0; ih<nh; ++ih){
    //   auto hit = track->GetHit(ih);
    //   Int_t layerId = hit->GetLayer();
    //   Double_t wire = hit->GetWire();
    //   Double_t res = hit->GetResidual();
    //   hddaq::cout << "       layer = " << layerId << ", wire = "
    //             << wire << ", res = " << res << std::endl;
    // }
    gEvDisp.DrawSdcInLocalTrack(track);
  }

  if(ntSdcIn == 0){
    hddaq::cout << "[Warning] SdcInTrack is empty!" << std::endl;
    return true;
  }

  //________________________________________________________
  //___ SdcOutDCHit
  DCAna.DecodeSdcOutHits();
  DCAna.TotCutSDC2(MinTotSDC2);
  DCAna.TotCutSDC3(MinTotSDC3);
  Double_t multi_SdcOut = 0.;
  for(Int_t plane=0; plane<NumOfLayersSdcOut; ++plane){
    const auto& cont = DCAna.GetSdcOutHC(plane);
    Int_t n = cont.size();
    multi_SdcOut += n;
    if(n > MaxMultiHitSdcOut) continue;
    for(const auto& hit: cont){
      Int_t layer = hit->LayerId();
      Int_t wire = hit->GetWire();
      Int_t mhit = hit->GetEntries();

      Bool_t is_good = false;
      for(Int_t j=0; j<mhit && !is_good; ++j){
        is_good = hit->IsGood(j);
	Double_t dtime = hit->GetDriftTime(j);
	//hddaq::cout << "SdcOut layer = " << layer << ", wire = " << wire
	//<< ", dt = " << dtime << std::endl;	
      }
      gEvDisp.DrawHitWire(layer, wire, is_good, is_good);
    }
  }
  multi_SdcOut /= (Double_t)NumOfLayersSdcOut;
  if(multi_SdcOut > MaxMultiHitSdcOut){
    hddaq::cout << "[Warning] SdcOutHits exceed MaxMultiHit "
                << multi_SdcOut << "/" << MaxMultiHitSdcOut << std::endl;
    // gEvDisp.GetCommand();
    return true;
  }

  //________________________________________________________
  //___ SdcOutTracking
#if UseTOF
  DCAna.TrackSearchSdcOut(TOFCont);
#else
  DCAna.TrackSearchSdcOut();
#endif

  Int_t ntSdcOut = DCAna.GetNtracksSdcOut();
  hddaq::cout << "[Info] ntSdcOut = " << ntSdcOut << std::endl;
  for(Int_t it=0; it<ntSdcOut; ++it){
    auto track = DCAna.GetTrackSdcOut(it);
    auto chisqr = track->GetChiSquare();
    hddaq::cout << "       " << it << "-th track, chi2 = "
                << chisqr << std::endl;
    // auto nh = track->GetNHit();
    // for(Int_t ih=0; ih<nh; ++ih){
    //   auto hit = track->GetHit(ih);
    //   Int_t layerId = hit->GetLayer();
    //   Double_t wire = hit->GetWire();
    //   Double_t res = hit->GetResidual();
    //   hddaq::cout << "       layer = " << layerId << ", wire = "
    //             << wire << ", res = " << res << std::endl;
    // }
    gEvDisp.DrawSdcOutLocalTrack(track);
  }
  if(ntSdcOut == 0){
    hddaq::cout << "[Warning] SdcOutTrack is empty!" << std::endl;
    return true;
  }

  //gEvDisp.Update();


  std::vector<ThreeVector> KpPCont, KpXCont;
  std::vector<Double_t> M2Cont;
  std::vector<Double_t> Chi2HypsCont;

  //________________________________________________________
  //___ HYPS Tracking
  static const auto StofOffset = gUser.GetParameter("StofOffset");
  // DCAna.SetMaxV0Diff(10.);
  DCAna.TrackSearchHyps();
  Bool_t through_target = false;
  Int_t ntHyps = DCAna.GetNTracksHyps();
  hddaq::cout << "[Info] ntHyps = " << ntHyps << std::endl;
  for(Int_t it=0; it<ntHyps; ++it){
    auto track = DCAna.GetHypsTrack(it);
    // track->Print();
    auto chisqr = track->GetChiSquare();
    hddaq::cout << "       " << it << "-th track, chi2 = "
                << chisqr << std::endl;
    const auto& postgt = track->PrimaryPosition();
    const auto& momtgt = track->PrimaryMomentum();
    Double_t path = track->PathLengthToTOF();
    Double_t p = momtgt.Mag();
    std::cout << "p = " << p << std::endl;
    gEvDisp.FillMomentum(p);
    if(chisqr > 20.) continue;
    if(TMath::Abs(postgt.x()) < 30.
       && TMath::Abs(postgt.y()) < 20.){
      through_target = true;
    }
    // MassSquare
    Double_t tofseg = track->TofSeg() + 1; // 1-origin
    for(const auto& hit: hodoAna.GetHitContainer("TOF")){
      Double_t seg = hit->SegmentId()+1;
      if(tofseg != seg) continue;
      Double_t stof = hit->CMeanTime()- time0+StofOffset;
      Int_t multi = hit->GetEntries();
      Double_t tdiff = 9999.;
      for (Int_t im = 0; im<multi; ++im) {
	Double_t mt=hit->CMeanTime(im);
	if (std::abs(mt-10)<std::abs(tdiff)) {
	  stof = mt- time0+StofOffset;
	  tdiff = mt-10;
	}
      }
      if(stof <= 0) continue;
      Double_t m2 = Kinematics::MassSquare(p, path, stof);
      std::cout << "m2 = " << m2 << std::endl;
      gEvDisp.FillMassSquare(m2);
      KpPCont.push_back(momtgt);
      KpXCont.push_back(postgt);
      M2Cont.push_back(m2);
      Chi2HypsCont.push_back(chisqr);
    }
  }
  if(KpPCont.size() == 0){
    hddaq::cout << "[Warning] Kp is empty!" << std::endl;
    // gEvDisp.GetCommand();
    // return true;
  }

  //________________________________________________________
  //___ DrawText
  TString buf;
  //std::stringstream ss; ss << trigger_flag;
  //buf = ss.str();
  //buf.ReplaceAll("0", ".").ReplaceAll("1", "!");
  //gEvDisp.DrawText(0.040, 0.960, Form("TrigFlag   %s", buf.Data()));
  //buf = "BH1Seg  ";
  //for(const auto& seg: BH1SegCont){
  //buf += Form(" %d", seg);
  //}
  //gEvDisp.DrawText(0.040, 0.920, buf);
  //buf = "BH2Seg  ";
  //for(const auto& seg: BH2SegCont){
  //buf += Form(" %d", seg);
  //}
  //gEvDisp.DrawText(0.040, 0.880, buf);
  // buf = "HTOFSeg  ";
  // for(const auto& seg: HTOFSegCont){
  //   buf += Form(" %d", seg);
  // }
  //gEvDisp.DrawText(0.040, 0.840, buf);
  buf = "TOFSeg  ";
  for(const auto& seg: TOFSegCont){
    buf += Form(" %d", seg);
  }
  gEvDisp.DrawText(0.040, 0.800, buf);
  // buf = "WCSeg  ";
  // for(const auto& seg: WCSegCont){
  //   buf += Form(" %d", seg);
  // }
  gEvDisp.DrawText(0.040, 0.760, buf);
  //buf = "BcOut"; gEvDisp.DrawText(0.040, 0.280, buf);
  //buf= "#chi^{2} = ";
  //for(Int_t i=0; i<ntBcOut; ++i){
  //buf += Form(" %.3f", DCAna.GetTrackBcOut(i)->GetChiSquare());
  //}
  //gEvDisp.DrawText(0.130, 0.280, buf);
  buf = "SdcIn"; gEvDisp.DrawText(0.040, 0.240, buf);
  buf = "#chi^{2} = ";
  for(Int_t i=0; i<ntSdcIn; ++i){
    buf += Form(" %.3f", DCAna.GetTrackSdcIn(i)->GetChiSquare());
  }
  gEvDisp.DrawText(0.130, 0.240, buf);
  buf = "SdcOut"; gEvDisp.DrawText(0.040, 0.20, buf);
  buf = "#chi^{2} = ";
  for(Int_t i=0; i<ntSdcOut; ++i){
    buf += Form(" %.3f", DCAna.GetTrackSdcOut(i)->GetChiSquare());
  }
  gEvDisp.DrawText(0.130, 0.200, buf);
  buf = "Hyps"; gEvDisp.DrawText(0.040, 0.16, buf);
  buf = "#chi^{2} = ";
  for(Int_t i=0; i<ntHyps; ++i){
    buf += Form(" %.3f", DCAna.GetHypsTrack(i)->GetChiSquare());
  }
  gEvDisp.DrawText(0.130, 0.160, buf);
  //gEvDisp.DrawText(0.680, 0.960, "BTOF");
  //gEvDisp.DrawText(0.860, 0.960, Form("%.3f", btof));

  hodoAna.DecodeHits<CFTFiberHit>("CFT");
  hodoAna.TimeCut("CFT", MinTimeCFT, MaxTimeCFT);
  hodoAna.AdcCut("CFT",  MinAdcCFT,  MaxAdcCFT);  

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
      Double_t adcLow = hit->GetRawHit()->GetAdcLow();
      Double_t mipHi  = hit->GetMipHigh();
      Double_t mipLow = hit->GetMipLow();      
      Double_t deLow  = hit->DeltaELowGain();

      //ADC 
      gEvDisp.DrawCFT_AdcCor(plane, seg, 0, (int)adccorHi);
      
      Int_t NhitT = hit->GetEntries(U);
      for(Int_t m = 0; m<NhitT; ++m){
	Double_t ctime = hit->GetCTUp(m);
	Double_t time = hit->GetTUp(m);

	if (adccorHi>0) {
	  gEvDisp.DrawCFT_Time(plane, seg, 0, ctime);
	  //gEvDisp.DrawCFT_Time(plane, seg, 1, time);
	}
	
	if (time>MinTimeCFT && time<MaxAdcCFT) {
	  gEvDisp.ShowHitFiber(plane, seg, adccorLow, -999);
	}
      }
    }
  }
  
  hodoAna.DecodeHits<HodoWaveformHit>("BGO");
  Int_t nhBGO=hodoAna.GetNHits("BGO");
  {
    int nc = 0;    
    const auto& U = HodoRawHit::kUp;
    for(Int_t i=0; i<nhBGO; ++i){
      const auto& hit = hodoAna.GetHit<HodoWaveformHit>("BGO", i);
      if(!hit) continue;
      Int_t seg   = hit->SegmentId();
      Int_t plane = hit->PlaneId();

      Int_t NhitT = hit->GetEntries(U);
      bool flagTime = false;
      for(Int_t m = 0; m<NhitT; ++m){
	Double_t time = hit->GetTUp(m);
	if (std::abs(time)<100)
	  flagTime = true;
      }
      
      Int_t NhitWF = hit->GetWaveformEntries(U);
      for(Int_t m = 0; m<NhitWF; ++m){
	std::pair<Double_t, Double_t> waveform = hit->GetWaveform(m);
      }


      Int_t Npulse = hit->GetNPulse(U);
      for(Int_t m = 0; m<Npulse; ++m){
	Double_t pulse_height = hit->GetPulseHeight(m);
	Double_t pulse_time   = hit->GetPulseTime(m);
	Double_t de           = hit->DeltaE(m);		

	if (flagTime && std::abs(pulse_time)<100) {
	  gEvDisp.ShowHitBGO(seg, de);	  
	}
      }

      Int_t ngr = hit->GetNGraph();
      if (ngr>0)
	nc++;

      //if (flagTime)
      //HF2 (1000*(plane+1)+206, seg, adccorHi);	
    }

    gEvDisp.SetBGOWaveformCanvas(nc);
    nc=1;
    for(Int_t i=0; i<nhBGO; ++i){
      const auto& hit = hodoAna.GetHit<HodoWaveformHit>("BGO", i);
      if(!hit) continue;
      Int_t seg   = hit->SegmentId();

      Int_t ngr = hit->GetNGraph();
      for (Int_t ig=0; ig<ngr; ig++) {
	TGraphErrors *gr = hit->GetTGraph(ig);
	gEvDisp.DrawBGOWaveform(nc, ig, seg, gr);      
      }
      Int_t Npulse = hit->GetNPulse(U);
      if (Npulse>0) {
	TF1 *func = hit->GetFitTF1();
	gEvDisp.DrawBGOFitFunc(nc, seg, func);      
      }
    }      
  }

 
  if(nhBGO == 0){
    hddaq::cout << "[Warning] BGO is no hit!" << std::endl;
    return true;
  }

  hodoAna.DecodeHits<CFTFiberHit>("PiID");
  {
    const auto& U = HodoRawHit::kUp;
    Int_t nh=hodoAna.GetNHits("PiID");
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit<CFTFiberHit>("PiID", i);
      if(!hit) continue;
      Int_t seg   = hit->SegmentId();
      Int_t NhitT = hit->GetEntries(U);
      bool flagTime = false;

      double time0 = -999;
      double ctime0 = -999;

      for(Int_t m = 0; m<NhitT; ++m){
	Double_t time = hit->GetTUp(m);
	if (time>-10 && time<20)
	  flagTime = true;
      }
      if (flagTime)
	gEvDisp.ShowHitPiID(seg);	  	
    }
  }

  const auto& CFTClCont = hodoAna.GetClusterContainer("CFT");
  DCAna.DecodeCFTHits(CFTClCont);
  DCAna.TrackSearchCFT();

  Int_t ntCFT=DCAna.GetNtracksCFT();
  
  for( Int_t i=0; i<ntCFT; ++i ){
    const auto& tp=DCAna.GetTrackCFT(i);
    Bool_t flagPTrack=false;
    gEvDisp.DrawCFTLocalTrack( tp, flagPTrack );
  }       

  std::vector <CFTParticle*> CFTPartCont;
  for( Int_t i=0; i<ntCFT; ++i ){
    const CFTLocalTrack *tp=DCAna.GetTrackCFT(i);
    CFTParticle * CFTPart = new CFTParticle(tp, &hodoAna);
    CFTPartCont.push_back(CFTPart);
  }

  for( Int_t i=0; i<ntCFT; ++i ){
    CFTParticle *CFTPart  = CFTPartCont[i];
    CFTPart->Calculate();

    bool flagPTrack = false;

    if (CFTPart->GetMass()>0.9) {
      //flagCFTProton = true;
      flagPTrack = true;
    }
    gEvDisp.DrawCFTLocalTrack_dE_E( CFTPart, flagPTrack );
    
  }

  
  //________________________________________________________
  //___ Reaction
  Bool_t is_good = false;
#if 0  
  if(KmPCont.size()==1 && KpPCont.size()==1){
    ThreeVector pkp = KpPCont[0];
    ThreeVector pkm = KmPCont[0];
    ThreeVector xkp = KpXCont[0];
    ThreeVector xkm = KmXCont[0];
    Double_t m2 = M2Cont[0];
    Double_t mass = TMath::QuietNaN();
    if(TMath::Abs(m2) < 0.15 && pkp.Mag() < 1.5) mass = PionMass;
    if(m2 > 0.15 && m2 < 0.35 && pkp.Mag() < 1.4) mass = KaonMass;
    if(m2 > 0.55) mass = ProtonMass;
    ThreeVector vertex = Kinematics::VertexPoint(xkm, xkp, pkm, pkp);
    Double_t closedist = Kinematics::CloseDist(xkm, xkp, pkm, pkp);
    LorentzVector LvKm(KmPCont[0], TMath::Sqrt(mass*mass+pkm.Mag2()));
    LorentzVector LvKp(KpPCont[0], TMath::Sqrt(mass*mass+pkp.Mag2()));
    LorentzVector LvP(0., 0., 0., ProtonMass);
    LorentzVector LvRp = LvKm+LvP-LvKp;
    ThreeVector MissMom = LvRp.Vect();
    Double_t MissMass = LvRp.Mag();
    hddaq::cout << "[Info] Vertex = " << vertex << std::endl;
    hddaq::cout << "[Info] MissingMomentum = " << MissMom << std::endl;
    gEvDisp.DrawVertex(vertex);
    gEvDisp.DrawMissingMomentum(MissMom, vertex);
    if(true
       && TMath::Abs(vertex.z()) < 200
       && TMath::Abs(vertex.x()) < 40.
       && TMath::Abs(vertex.y()) < 20.
       && closedist < 20.
       // && through_target
    ){
      gEvDisp.DrawText(0.680, 0.920, "pK18");
      gEvDisp.DrawText(0.860, 0.920, Form("%.3f", pkm.Mag()));
      gEvDisp.DrawText(0.680, 0.880, "pHyps");
      gEvDisp.DrawText(0.860, 0.880, Form("%.3f", pkp.Mag()));
      gEvDisp.DrawText(0.680, 0.840, "MassSquared");
      gEvDisp.DrawText(0.860, 0.840, Form("%.3f", m2));
      gEvDisp.DrawText(0.660, 0.280, "CloseDist");
      gEvDisp.DrawText(0.770, 0.280, Form("%.2f", closedist));
      gEvDisp.DrawText(0.660, 0.240, "Vertex");
      gEvDisp.DrawText(0.770, 0.240, Form("(%.2f, %.2f, %.2f)",
                                          vertex.X(), vertex.Y(), vertex.Z()));
      gEvDisp.DrawText(0.660, 0.200, "MissMom");
      gEvDisp.DrawText(0.770, 0.200, Form("(%.3f, %.3f, %.3f)",
                                          MissMom.X(), MissMom.Y(), MissMom.Z()));
      gEvDisp.DrawText(0.660, 0.160, "MissMass");
      gEvDisp.DrawText(0.770, 0.160, Form("%.4f", MissMass));
      if(true
         // && mass == KaonMass
         // && pkp.z() > 0
      ){
        // gEvDisp.GetCommand();
        is_good = true;
      }
    }
  }
#endif
  is_good = true;
  gEvDisp.Update();
  gEvDisp.GetCommand();
  hddaq::cout << "[Info] IsGood = " << is_good << std::endl;

  if(is_good){
#if SAVEPDF
    gEvDisp.Print(gUnpacker.get_run_number(),
                  gUnpacker.get_event_number());
#else
    //gSystem->Sleep(5000);
#endif
  }
#if DEBUG
  hddaq::cout << "Wait a moment ..." << std::endl;
  gSystem->Sleep(5000);
#endif
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingEnd()
{
  //gEvDisp.GetCommand();
  gEvDisp.EndOfEvent();
  // if(utility::UserStop()) gEvDisp.Run();
  // gEvDisp.Run();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan:: InitializeHistograms()
{
  gUnpacker.disable_istream_bookmark();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
#if SAVEPDF
  gEvDisp.SetSaveMode();
#endif
  return
    (InitializeParameter<DCGeomMan>("DCGEO")
     && InitializeParameter<DCDriftParamMan>("DCDRFT")
     && InitializeParameter<DCTdcCalibMan>("DCTDC")
     && InitializeParameter<HodoParamMan>("HDPRM")
     && InitializeParameter<HodoPHCMan>("HDPHC")
     && InitializeParameter<FieldMan>("FLDMAP")
     //&& InitializeParameter<K18TransMatrix>("K18TM")
     //&& InitializeParameter<BH2Filter>("BH2FLT")
     && InitializeParameter<UserParamMan>("USER")
     && InitializeParameter<CFTPedCorMan>("CFTPED") 
     && InitializeParameter<CFTPosParamMan>("CFTPOS")
     && InitializeParameter<TemplateFitMan>("BGOTEMP") 
     && InitializeParameter<BGOCalibMan>("BGOCALIB")
     && InitializeParameter<TAGPLMatch>("TAGPLMTH")
     && InitializeParameter<CATCHPidMan>("CATCHPID")      
     && InitializeParameter<EventSelectMan>("EVSELECT") 
     && InitializeParameter<EventDisplay>());
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
