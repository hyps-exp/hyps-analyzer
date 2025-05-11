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
#include "HodoWaveformHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "DCAnalyzer.hh"
#include "CFTLocalTrack.hh"
#include "CFTParticle.hh"
#include "HypsLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "DCGeomMan.hh"
#include "CFTPedCorMan.hh"
#include "CFTPosParamMan.hh"
#include "CATCHPidMan.hh"
#include "TemplateFitMan.hh"
#include "BGOCalibMan.hh"

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

  double dE[NumOfPlaneCFT][MaxDepth] , dE_max[NumOfPlaneCFT][MaxDepth] ;
  Int_t  MaxSegment[NumOfPlaneCFT][MaxDepth];
  Double_t  MaxAdcLow[NumOfPlaneCFT][MaxDepth];  
  
  int    ntCFT;
  double theta_cft[MaxDepth];
  int    nhit_phi[MaxDepth];
  int    nhit_uv[MaxDepth];   
  double Total_dE[MaxDepth],    Total_dE_max[MaxDepth];
  double Total_dEphi[MaxDepth], Total_dEphi_max[MaxDepth];
  double Total_dEuv[MaxDepth],  Total_dEuv_max[MaxDepth];

  // BGO
  int segBGOt[NumOfSegBGO];// matched to track
  double energybgo[NumOfSegBGO];

  // PiID counter
  int segPiIDt[MaxDepth];// matched to track


  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum     = 0;
  spill     = 0;
  ntCFT     = 0;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<NumOfPlaneCFT; ++it){
    for(Int_t m=0; m<MaxDepth; ++m){
      dE[it][m] = qnan;
      dE_max[it][m] = qnan;
      MaxSegment[it][m] = qnan;
      MaxAdcLow[it][m] = qnan;  

    }
  }

  for(Int_t m=0; m<MaxDepth; ++m){  
    theta_cft[m] = qnan;
    nhit_phi[m]  = -1;
    nhit_uv[m]   = -1;   

    Total_dE[m]        = qnan;
    Total_dE_max[m]    = qnan;
    Total_dEphi[m]     = qnan;
    Total_dEphi_max[m] = qnan;
    Total_dEuv[m]      = qnan;
    Total_dEuv_max[m]  = qnan;
  }

  for(Int_t m=0; m<NumOfSegBGO; ++m){    
    segBGOt[m] = -1;
    energybgo[m] = qnan;
  }
  for(Int_t m=0; m<MaxDepth; ++m){    
    segPiIDt[m];
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
  rawData.DecodeHits("PiID");
  hodoAna.DecodeHits<CFTFiberHit>("CFT");
  hodoAna.DecodeHits<CFTFiberHit>("PiID");
  
  hodoAna.TimeCut("CFT", MinTimeCFT, MaxTimeCFT);
  hodoAna.AdcCut("CFT",  MinAdcCFT,  MaxAdcCFT);  
  const auto& CFTClCont = hodoAna.GetClusterContainer("CFT");
  DCAna.DecodeCFTHits(CFTClCont);
  DCAna.TrackSearchCFT();

  Int_t ntCFT=DCAna.GetNtracksCFT();
  event.ntCFT = ntCFT;

  rawData.DecodeHits("BGO");
  hodoAna.DecodeHits<HodoWaveformHit>("BGO");

  std::vector <CFTParticle*> CFTPartCont;
  for( Int_t i=0; i<ntCFT; ++i ){
    const CFTLocalTrack *tp=DCAna.GetTrackCFT(i);
    CFTParticle * CFTPart = new CFTParticle(tp, &hodoAna);
    CFTPartCont.push_back(CFTPart);
  }

  for( Int_t i=0; i<ntCFT; ++i ){
    CFTParticle *CFTPart  = CFTPartCont[i];
    CFTPart->Calculate();

    const CFTLocalTrack *tp=CFTPart->GetTrack();

    Int_t nh   = tp->GetNHit();
    Int_t nhUV = tp->GetNHitUV();
    Double_t chisqrXY=tp->GetChiSquareXY();
    Double_t chisqrXYZ=tp->GetChiSquareZ();
    //Double_t vtx_z =tp->GetVtxZ();
    Double_t theta =tp->GetThetaCFT();
    Int_t xyFlag = tp->GetCFTxyFlag();
    Int_t zFlag  = tp->GetCFTzFlag() ;
       
    ThreeVector Pos0 = tp->GetPos0();
    ThreeVector Dir = tp->GetDir();
    Double_t A=(Dir.x()*Dir.x()+Dir.y()*Dir.y());

    // aka 
    Double_t D=(Dir.x()*Dir.x()+Dir.y()*Dir.y()+Dir.z()*Dir.z());

    Double_t phi = -999.;
    if(Dir.x()>=0 && Dir.y()>=0){
      phi = acos(Dir.x()/sqrt(A))*TMath::RadToDeg();
    }//0~90
    else if (Dir.x()<0 && Dir.y()>=0){
      phi = acos(Dir.x()/sqrt(A))*TMath::RadToDeg();
    }//90~180
    else if (Dir.x()<0 && Dir.y()<0){
      phi = 360. - acos(Dir.x()/sqrt(A))*TMath::RadToDeg();
    }//180~270
    else if (Dir.x()>=0 && Dir.y()<0){
      phi = 360. - acos(Dir.x()/sqrt(A))*TMath::RadToDeg(); ;
    }//270~360
    else{}

    // vertex
    ThreeVector  bPos(0., 0., 0.);
    ThreeVector  bdir(0., 0., 1.);
    Double_t dist = 1000.;
    if(Pos0.x()>-500.){
      ThreeVector vtx = Kinematics::VertexPoint3D(bPos, Pos0, bdir, Dir, dist);
      //event.vtx_x[i] = vtx.x();
      //event.vtx_y[i] = vtx.y();
      //event.vtx_z[i] = vtx.z();
    }

    event.nhit_phi[i]   = nh;
    event.nhit_uv[i]    = nhUV;
    event.theta_cft[i]  = theta;    
    event.Total_dE[i]   = tp->GetTotalSumdE()   ;
    event.Total_dEphi[i]= tp->GetTotalSumdEphi();
    event.Total_dEuv[i] = tp->GetTotalSumdEuv ();
    event.Total_dE_max[i]   = tp->GetTotalMaxdE()   ;
    event.Total_dEphi_max[i]= tp->GetTotalMaxdEphi();
    event.Total_dEuv_max[i] = tp->GetTotalMaxdEuv ();

    
    // straight layer
    for(Int_t ip=0; ip<nh; ip++){
      //CFTFiberCluster *hit = tp->GetHit(ip);
      const auto& hit = tp->GetHit(ip);
      //Int_t layer = hit->GetTrackingLayer() - layerId_U1;
      Int_t layer = hit->PlaneId();
      Int_t seg_max = hit->MaxSegment();
      event.dE[layer][i]  = hit->TotalDeltaE();
      event.MaxSegment[layer][i]  = hit->MaxSegment();
      event.MaxAdcLow[layer][i]  = hit->MaxAdcLow();      
    }

    // spiral layer
    Double_t xmin, xmax;
    Double_t ymin, ymax;
    for(Int_t ip=0; ip<nhUV; ip++){
      const auto& hit = tp->GetHitUV(ip);
      //CFTFiberCluster *hit = tp->GetHitUV(ip);
      //Int_t layer = hit->GetTrackingLayer() - layerId_U1;
      Int_t layer = hit->PlaneId();
      Int_t seg_max = hit->MaxSegment();
      event.dE[layer][i]  = hit->TotalDeltaE();
      event.MaxSegment[layer][i]  = hit->MaxSegment();
      event.MaxAdcLow[layer][i]  = hit->MaxAdcLow();      
    }    

    Int_t segBGOt  = CFTPart->GetTrackBGOSeg(); // BGO  track segment
    event.segBGOt[i] = segBGOt;
    Double_t bgo_energy = CFTPart->GetBGOEnergy();
    if (segBGOt>=0 && segBGOt < NumOfSegBGO)
      event.energybgo[segBGOt] = bgo_energy;
    
    Int_t segPiIDt = CFTPart->GetTrackPiIDSeg();// PiID track segment

    double dE = tp->GetTotalMaxdEphi()*sin(theta*TMath::DegToRad())/(double)(nh)
      + tp->GetTotalMaxdEuv()/(double)(nhUV);

    if (segPiIDt>=0)
      HF2(11, bgo_energy, dE);
    else
      HF2(10, bgo_energy, dE);
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

  HB2(10, "#Delta E - E (w/o PiID)", 100, 0., 200., 100, 0, 10);
  HB2(11, "#Delta E - E (w/ PiID)", 100, 0., 200., 100, 0, 10);  
  

  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spill",     &event.spill,     "spill/I");
  //Trig
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  tree->Branch("ntCFT",     &event.ntCFT,    "ntCFT/I");

  tree->Branch("theta_cft",    event.theta_cft,    "theta[ntCFT]/D");
  tree->Branch("nhit_phi",     event.nhit_phi,    "nhit_phi[ntCFT]/I");
  tree->Branch("nhit_uv",      event.nhit_uv,     "nhit_uv[ntCFT]/I");
  
  tree->Branch("dE",      event.dE,       Form("dE[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("dE_max",  event.dE_max,   Form("dE_max[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("MaxSegment",  event.MaxSegment,   Form("MaxSegment[%d][%d]/I", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("MaxAdcLow",   event.MaxAdcLow,    Form("MaxAdcLow[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );    
  
  tree->Branch("Total_dE",    event.Total_dE,    "totaldE[ntCFT]/D");
  tree->Branch("Total_dEphi", event.Total_dEphi, "totaldEphi[ntCFT]/D");
  tree->Branch("Total_dEuv",  event.Total_dEuv,  "totaldEuv[ntCFT]/D");
  tree->Branch("Total_dE_max",    event.Total_dE_max,    "totaldE_max[ntCFT]/D");
  tree->Branch("Total_dEphi_max", event.Total_dEphi_max, "totaldEphi_max[ntCFT]/D");
  tree->Branch("Total_dEuv_max",  event.Total_dEuv_max,  "totaldEuv_max[ntCFT]/D");

  // BGO
  //tree->Branch("nhBGO",     &event.nhBGO,    "nhBGO/I");
  //tree->Branch("segBGO",    event.segBGO,    "segBGO[nhBGO]/I");
  tree->Branch("segBGOt",   event.segBGOt,   "segBGOt[24]/I");
  //tree->Branch("adcBGO",    event.adcbgo,    "adcbgo[24]/D");
  tree->Branch("energyBGO", event.energybgo, "energybgo[24]/D");
  //tree->Branch("tdcBGO",    event.tdcbgo,    "adcbgo[24]/D");

  // PiID
  tree->Branch("segPiIDt",   event.segPiIDt,   "segPiIDt[ntCFT]/I");

  
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
     InitializeParameter<CFTPosParamMan>("CFTPOS") &&
     InitializeParameter<TemplateFitMan>("BGOTEMP") &&
     InitializeParameter<BGOCalibMan>("BGOCALIB")   &&
     InitializeParameter<CATCHPidMan>("CATCHPID") 
     );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
