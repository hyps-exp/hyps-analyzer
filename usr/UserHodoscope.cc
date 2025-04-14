// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <sstream>

#include <UnpackerManager.hh>

// #include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "TAGPLMatch.hh"
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
#include "HypsLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "DCGeomMan.hh"

// #define TimeCut    1 // in cluster analysis
#define FHitBranch 0 // make FiberHit branches (becomes heavy)
#define HodoHitPos 0
#define SDCoutreq 0

namespace
{
using namespace root;
const auto qnan = TMath::QuietNaN();
auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto& gRM       = RMAnalyzer::GetInstance();
auto& gUser     = UserParamMan::GetInstance();
auto& gTAGPLMth = TAGPLMatch::GetInstance();
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t spill;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Int_t rfnhits;
  Double_t rft[NumOfSegRF][MaxDepth];

  Int_t tagsffnhits;
  Int_t tagsffhitpat[MaxHits];
  Double_t tagsfft[NumOfSegTagSF][MaxDepth];

  Int_t tagsfbnhits;
  Int_t tagsfbhitpat[MaxHits];
  Double_t tagsfbt[NumOfSegTagSF][MaxDepth];

  Int_t tagplnhits;
  Int_t tagplhitpat[MaxHits];
  Double_t tagpla[NumOfSegTagPL];
  Double_t tagplt[NumOfSegTagPL][MaxDepth];


  Int_t t0nhits;
  Double_t t0la[NumOfSegT0];
  Double_t t0lt[NumOfSegT0][MaxDepth];
  Double_t t0ra[NumOfSegT0];
  Double_t t0rt[NumOfSegT0][MaxDepth];

  Int_t sacnhits;
  Int_t sachitpat[MaxHits];
  Double_t saca[NumOfSegSAC];
  Double_t sact[NumOfSegSAC][MaxDepth];

  Int_t sac1nhits;
  Int_t sac2nhits;

  Int_t e_vetonhits;
  Double_t e_vetola[NumOfSegE_Veto];
  Double_t e_vetolt[NumOfSegE_Veto][MaxDepth];
  Double_t e_vetora[NumOfSegE_Veto];
  Double_t e_vetort[NumOfSegE_Veto][MaxDepth];


  Int_t tofnhits;
  Int_t tofhitpat[MaxHits];
  Double_t tofua[NumOfSegTOF];
  Double_t tofut[NumOfSegTOF][MaxDepth];
  Double_t tofda[NumOfSegTOF];
  Double_t tofdt[NumOfSegTOF][MaxDepth];


  ////////// Normalized
  Double_t rf[NumOfSegRF][MaxDepth];
  Double_t crf[NumOfSegRF][MaxDepth];

  Double_t t0[NumOfSegT0][MaxDepth];
  Double_t ct0[NumOfSegT0][MaxDepth];

  Double_t sacmt[NumOfSegSAC][MaxDepth];
  Double_t sacde[NumOfSegSAC];

  Double_t e_veto[NumOfSegE_Veto][MaxDepth];
  Double_t ce_veto[NumOfSegE_Veto][MaxDepth];

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
  evnum       = 0;
  spill       = 0;
  rfnhits     = 0;
  tagsffnhits=0;
  tagsfbnhits=0;
  tagplnhits=0;
  t0nhits     = 0;
  sacnhits    = 0;
  sac1nhits   = 0;
  sac2nhits   = 0;
  e_vetonhits = 0;
  tofnhits    = 0;
  Time0Seg    = qnan;
  deTime0     = qnan;
  Time0       = qnan;
  CTime0      = qnan;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<MaxHits; ++it){
    tagsffhitpat[it]=-1;
    tagsfbhitpat[it]=-1;
    tagplhitpat[it]=-1;
    sachitpat[it]    = -1;
    tofhitpat[it]    = -1;
  }

  for(Int_t it=0; it<NumOfSegRF; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      rf[it][m]        = qnan;
      crf[it][m]       = qnan;
      rft[it][m]  = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegTagSF; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      tagsfft[it][m]  = qnan;
      tagsfbt[it][m]  = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegTagPL; it++){
    tagpla[it]=qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tagplt[it][m]  = qnan;
      tagplt[it][m]  = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegT0; it++){
    t0la[it]  = qnan;
    t0ra[it]  = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      t0lt[it][m]  = qnan;
      t0rt[it][m]  = qnan;
      t0[it][m]        = qnan;
      ct0[it][m]       = qnan;
    }
  }


  for(Int_t it=0; it<NumOfSegSAC; it++){
    saca[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      sact[it][m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegE_Veto; it++){
    e_vetola[it]  = qnan;
    e_vetora[it]  = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      e_vetolt[it][m]  = qnan;
      e_vetort[it][m]  = qnan;
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

  Int_t    nhRF;
  Int_t    csRF[NumOfSegRF*MaxDepth];
  Double_t RFSeg[NumOfSegRF*MaxDepth];
  Double_t tRF[NumOfSegRF*MaxDepth];

  // for HodoParam
  Double_t rft[NumOfSegRF][MaxDepth];
  Double_t tRFSeg[NumOfSegRF][MaxDepth];

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

  Int_t    nhE_veto;
  Int_t    csE_veto[NumOfSegE_Veto*MaxDepth];
  Double_t E_vetoSeg[NumOfSegE_Veto*MaxDepth];
  Double_t tE_veto[NumOfSegE_Veto*MaxDepth];

  // for HodoParam
  Double_t e_vetola[NumOfSegE_Veto];
  Double_t e_vetolt[NumOfSegE_Veto][MaxDepth];
  Double_t e_vetora[NumOfSegE_Veto];
  Double_t e_vetort[NumOfSegE_Veto][MaxDepth];
  Double_t ltE_vetoSeg[NumOfSegE_Veto][MaxDepth];
  Double_t rtE_vetoSeg[NumOfSegE_Veto][MaxDepth];


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
  nhRF     = 0;
  nhT0     = 0;
  nhSac    = 0;
  nhBh2    = 0;
  nhTof    = 0;
  nhE_veto = 0;
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

  for(Int_t it=0; it<NumOfSegRF; it++){
    RFSeg[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      rft[it][m]           = qnan;
      tRF[MaxDepth*it + m] = qnan;
      tRFSeg[it][m]        = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegT0; it++){
    t0la[it]  = qnan;
    t0ra[it]  = qnan;
    T0Seg[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      t0lt[it][m]          = qnan;
      t0rt[it][m]          = qnan;
      tT0[MaxDepth*it + m] = qnan;
      ltT0Seg[it][m]       = qnan;
      rtT0Seg[it][m]       = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegSAC; it++){
    SacSeg[it] = qnan;
    deSac[it]  = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tSac[MaxDepth*it + m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegE_Veto; it++){
    e_vetola[it]  = qnan;
    e_vetora[it]  = qnan;
    E_vetoSeg[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      e_vetolt[it][m]          = qnan;
      e_vetort[it][m]          = qnan;
      tE_veto[MaxDepth*it + m] = qnan;
      ltE_vetoSeg[it][m]       = qnan;
      rtE_vetoSeg[it][m]       = qnan;
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
  RFHid     = 10000,
  TagSFHid  = 20000,
  TagPLHid  = 30000,
  T0Hid     = 50000,
  SACHid    = 60000,
  E_VetoHid = 70000,
  TOFHid    = 80000,
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
  static const auto MinTdcRF     = gUser.GetParameter("TdcRF", 0);
  static const auto MaxTdcRF     = gUser.GetParameter("TdcRF", 1);
  static const auto MinTdcSF  = gUser.GetParameter("TdcSF", 0);
  static const auto MaxTdcSF  = gUser.GetParameter("TdcSF", 1);
  static const auto MinTdcPL  = gUser.GetParameter("TdcPL", 0);
  static const auto MaxTdcPL  = gUser.GetParameter("TdcPL", 1);
  static const auto MinTdcT0     = gUser.GetParameter("TdcT0", 0);
  static const auto MaxTdcT0     = gUser.GetParameter("TdcT0", 1);
  static const auto MinTdcSAC    = gUser.GetParameter("TdcSAC", 0);
  static const auto MaxTdcSAC    = gUser.GetParameter("TdcSAC", 1);
  static const auto MinTdcE_Veto = gUser.GetParameter("TdcE_Veto", 0);
  static const auto MaxTdcE_Veto = gUser.GetParameter("TdcE_Veto", 1);
  static const auto MinTdcTOF    = gUser.GetParameter("TdcTOF", 0);
  static const auto MaxTdcTOF    = gUser.GetParameter("TdcTOF", 1);
#if HodoHitPos
  static const auto PropVelBH2   = gUser.GetParameter("PropagationBH2");
#endif

  RawData rawData;
  HodoAnalyzer hodoAna(rawData);

#if SDCoutreq
  for(const auto& name: DCNameList.at("SdcIn")) rawData.DecodeHits(name);
  for(const auto& name: DCNameList.at("SdcOut")) rawData.DecodeHits(name);
  DCAnalyzer   DCAna(rawData);
  DCAna.DecodeSdcOutHits();
  DCAna.TotCutSDC2(30);
  DCAna.TotCutSDC3(40);
  DCAna.TrackSearchSdcOut();
  Int_t nt=DCAna.GetNtracksSdcOut();
  //std::cout<<nt<<std::endl;
  if(nt<1) return true;
#endif

    

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

  ///// RF
  rawData.DecodeHits("RF");
  {
    Int_t rf_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("RF");
    Int_t nh = cont.size();
    HF1(RFHid, Double_t(nh));
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId();


      Bool_t is_hit = false;
      Int_t m_l = 0;
      for(const auto& T: hit->GetArrayTdc()){
        HF1(RFHid +100*(seg+1) +3, T);
        if(m_l < MaxDepth){
          event.rft[seg][m_l] = T;
          dst.rft[seg][m_l]   = T;
          ++m_l;
        }
        if(MinTdcRF < T && T < MaxTdcRF) is_hit = true;
      }


    }
    event.rfnhits = rf_nhits;
  }

      ///// Tag-SF
  rawData.DecodeHits("TAG-SF");

  {
    Int_t tagsf_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("TAG-SF");
    const auto& U= HodoRawHit::kUp;
    Int_t nh = cont.size();
    //HF1(T0Hid, Double_t(nh));
    Int_t sffnhits = 0, sfbnhits = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t plane =hit->PlaneId();
      Int_t seg = hit->SegmentId();
      //HF1(TagSFHid+1, seg+0.5);
      // Left
      Bool_t is_hit_l = false;
      Int_t m_l = hit->GetSizeTdcLeading();
      for(Int_t j=0;j<m_l;++j){
        Double_t leading =hit->GetTdcLeading(U,j);
	HF1(TagSFHid + plane +3, leading);
	HF1(TagSFHid + 100*(seg+1) + 3 + plane,leading);
	if(plane ==0) event.tagsfft[seg][j]= leading;
	if(plane ==1) event.tagsfbt[seg][j]= leading;
	if(MinTdcSF<leading && leading<MaxTdcSF){
	  HF1(TagSFHid + plane + 5, seg+0.5);
	  if(plane ==0) event.tagsffhitpat[sffnhits++]= seg;
	  if(plane ==1) event.tagsfbhitpat[sfbnhits++]= seg;
	}

      }

    }
    event.tagsffnhits = sffnhits;
    event.tagsfbnhits = sfbnhits;
  }

      ///// Tag-PL

  rawData.DecodeHits("TAG-PL");

  {
    Int_t tagpl_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("TAG-PL");
    const auto& U= HodoRawHit::kUp;
    Int_t nh = cont.size();
    //HF1(TagPLHid+0, Double_t(nh));
    //HF1(T0Hid, Double_t(nh));
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId();
      HF1(TagPLHid+5, seg+0.5);
      Double_t a=hit->GetAdc();
      event.tagpla[seg]=a;
      Bool_t is_hit_l = false;
      Int_t m_l = hit->GetSizeTdcLeading();
      for(Int_t j=0;j<m_l;++j){
        Double_t leading =hit->GetTdcLeading(U,j);
	HF1(TagPLHid +3, leading);
	HF1(TagPLHid + 100*(seg+1) + 3 ,leading);
	event.tagplt[seg][j]= leading;
	if(MinTdcPL<leading && leading<MaxTdcPL){
	  HF1(TagPLHid + 6, seg+0.5);
	  event.tagplhitpat[tagpl_nhits++]= seg;
	}

      }

    }
    event.tagplnhits = tagpl_nhits;
    HF1(TagPLHid+0, tagpl_nhits);
  }

  //Tag-Hodoana
  std::vector<int> PLCand;
  {
  hodoAna.DecodeHits("TAG-PL");

    Int_t nh=hodoAna.GetNHits("TAG-PL");
    Int_t nseg_goodtime=0;
    for(Int_t i=0;i<nh;++i){
      const auto& hit =hodoAna.GetHit("TAG-PL",i);
      if(!hit) continue;
      Int_t seg =hit->SegmentId();
      /*
	Double_t a =hit->GetAUp();

      */
      bool is_hit_time =false;
      Int_t n_mhit =hit->GetEntries();
      for(Int_t m=0;m<n_mhit;++m){
	Double_t t=hit->GetTUp(m);
	//Double_t ct= hit->GetCTUp(m);
	//HF1(TagPLHid +100*(seg+1) +13, t);
	HF1(TagPLHid +13, t);
	if(fabs(t)<5.0) is_hit_time=true;
	//if(fabs(t)<300.0) is_hit_time=true;
      }
      if(is_hit_time){
	nseg_goodtime++;
	HF1(TagPLHid +16, seg+0.5);
	PLCand.push_back(seg);
      }

    }
    HF1(TagPLHid+10, Double_t(PLCand.size()));

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
	//std::cout<<"(seg, plane) = ("<<seg<<", "<<planename<<")"<<std::endl;
	/*
	  Double_t a =hit->GetAUp();

	*/
	bool is_hit_time =false;
	Int_t n_mhit =hit->GetEntries();
	for(Int_t m=0;m<n_mhit;++m){
	  Double_t t=hit->GetTUp(m);
	  //Double_t ct= hit->GetCTUp(m);
	  //HF1(TagPLHid +100*(seg+1) +13, t);
	  if(fabs(t)<5) is_hit_time=true;
	  HF1(TagSFHid +13 +plane, t);
	}
	if(is_hit_time){
	  //if(plane==0) SFFhit.push_back(seg);
	  //if(plane==1) SFBhit.push_back(seg);
	}

      }

      hodoAna.TimeCut("TAG-SF",-5,5);
      //hodoAna.TimeCut("TAG-SF",-100,100);

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
	HF1(TagSFHid+10+plane,cl->ClusterSize());

	if(plane==0 && cl->ClusterSize()<4) SFFhit.push_back(ms);
	if(plane==1 && cl->ClusterSize()<4) SFBhit.push_back(ms);

	HF1(TagSFHid +15 +plane, ms);
	bool is_plmth=false;
	if(PLCand.size()==0) is_plmth=true;
	for(Int_t j=0;j<PLCand.size();j++){
	  if(gTAGPLMth.Judge(ms,PLCand[j])) is_plmth=true;
	}
	if(is_plmth && cl->ClusterSize()<4){
	  HF1(TagSFHid +25 +plane, ms);
	  if(plane==0) SFFCand.push_back(ms);
	  if(plane==1) SFBCand.push_back(ms);
	}


      }
      HF1(TagSFHid+0, Double_t(ncl1));
      HF1(TagSFHid+1, Double_t(ncl2));
      HF1(TagSFHid+20, Double_t(SFFCand.size()));
      HF1(TagSFHid+21, Double_t(SFBCand.size()));

    }

    std::vector<double> SFFCand_final;
    std::vector<double> SFBCand_final;

    for(int i=0;i<SFFhit.size();i++){
      for(int j=0;j<SFBhit.size();j++){
    	HF2(TagSFHid+17,SFFhit[i],SFBhit[j]);
      }
    }

    for(int i=0;i<SFFCand.size();i++){
      for(int j=0;j<SFBCand.size();j++){
	HF2(TagSFHid+27,SFFCand[i],SFBCand[j]);
	if(fabs(SFFCand[i]-SFBCand[j])<3){
	  SFFCand_final.push_back(SFFCand[i]);
	  SFBCand_final.push_back(SFBCand[j]);
	}
      }
    }
    HF1(TagSFHid+30,Double_t(SFFCand_final.size()));
    HF1(TagSFHid+31,Double_t(SFBCand_final.size()));

    if(SFFCand_final.size()==1){
      HF2(TagSFHid+37,SFFCand[0],SFBCand[0]);
      const double eparf[3]={1.42893,0.0350856,-0.000132564};
      const double eparb[3]={1.42415,0.0351546,-0.000136274};
      const double offset_b=0.6421;

      double egamf=eparf[0]+eparf[1]*SFFCand[0]+eparf[2]*SFFCand[0]*SFFCand[0];
      double SFBpos=SFBCand[0]+offset_b;
      double egamb=eparb[0]+eparb[1]*SFBpos+eparb[2]*SFBpos*SFBpos;

      HF1(TagSFHid+35,egamf);
      HF1(TagSFHid+36,egamb);

    }

    for(int i=0;i<PLCand.size();i++){
      for(int j=0;j<SFFhit.size();j++){
	HF2(TagSFHid+18,PLCand[i],SFFhit[j]);
      }
      for(int k=0;k<SFBhit.size();k++){
	HF2(TagSFHid+19,PLCand[i],SFBhit[k]);
      }
    }

  ///// T0
  rawData.DecodeHits("T0");
  {
    Int_t t0_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("T0");
    Int_t nh = cont.size();
    HF1(T0Hid, Double_t(nh));
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId();

      // Left
      Int_t Al = hit->GetAdcLeft();
      HF1(T0Hid+100*(seg+1)+1, Al);
      event.t0la[seg] = Al;
      dst.t0la[seg]   = Al;

      Bool_t is_hit_l = false;
      Int_t m_l = 0;
      for(const auto& T: hit->GetArrayTdcLeft()){
        HF1(T0Hid +100*(seg+1) +3, T);
        if(m_l < MaxDepth){
          event.t0lt[seg][m_l] = T;
          dst.t0lt[seg][m_l]   = T;
          ++m_l;
        }
        if(MinTdcT0 < T && T < MaxTdcT0) is_hit_l = true;
      }
      if(is_hit_l) HF1(T0Hid+100*(seg+1)+5, Al);
      else         HF1(T0Hid+100*(seg+1)+7, Al);

      // Right
      Int_t Ar = hit->GetAdcRight();
      HF1(T0Hid+100*(seg+1)+2, Ar);
      event.t0ra[seg] = Ar;
      dst.t0ra[seg]   = Ar;

      Bool_t is_hit_r = false;
      Int_t m_r = 0;
      for(const auto& T: hit->GetArrayTdcRight()){
        HF1(T0Hid +100*(seg+1) +4, T);
        if(m_r < MaxDepth){
          event.t0rt[seg][m_r] = T;
          dst.t0rt[seg][m_r]   = T;
          ++m_r;
        }
        if(MinTdcT0 < T && T < MaxTdcT0)  is_hit_r = true;
      }
      if(is_hit_r) HF1(T0Hid+100*(seg+1)+6, Ar);
      else         HF1(T0Hid+100*(seg+1)+8, Ar);

      // Hitpat
      if(is_hit_l || is_hit_r){
        ++nh1;
      }
      if(is_hit_l && is_hit_r){
        ++nh2;
      }
    }
    HF1(T0Hid+2, nh1); HF1(T0Hid+4, nh2);
    event.t0nhits = t0_nhits;
  }


#if 1
  ///// SAC
  rawData.DecodeHits("SAC");
  {
    Int_t sac_nhits = 0;
    Int_t sac1_nhits = 0;
    Int_t sac2_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("SAC");
    Int_t nh = cont.size();
    HF1(SACHid, nh);
    Int_t nh1 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(SACHid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      HF1(SACHid+100*seg+1, A);
      event.saca[seg-1] = A;
      Bool_t is_hit = false;
      Int_t m = 0;
      for(const auto& T: hit->GetArrayTdcLeading()){
        HF1(SACHid+100*seg+3, T);
        if(m < MaxDepth) event.sact[seg-1][m++] = T;
        if(MinTdcSAC < T && T < MaxTdcSAC){
	  is_hit = true;
	  if(seg==1)++sac1_nhits;
	  if(seg==2)++sac2_nhits;
	}
      }
      if(is_hit) HF1(SACHid+100*seg+5, A);
      else       HF1(SACHid+100*seg+7, A);
      // Hitpat
      if(is_hit){
        event.sachitpat[sac_nhits++] = seg;
        ++nh1; HF1(SACHid+3, seg-0.5);
      }
    }
    HF1(SACHid+2, nh1);
    event.sacnhits = sac_nhits;
    event.sac1nhits = sac1_nhits;
    event.sac2nhits = sac2_nhits;
  }
#endif // BH2, SAC

  ///// E_Veto
  rawData.DecodeHits("E_Veto");
  {
    Int_t e_veto_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("E_Veto");
    Int_t nh = cont.size();
    HF1(E_VetoHid, Double_t(nh));
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId();

      // Left
      Int_t Al = hit->GetAdcLeft();
      HF1(E_VetoHid+100*(seg+1)+1, Al);
      event.e_vetola[seg] = Al;
      dst.e_vetola[seg]   = Al;

      Bool_t is_hit_l = false;
      Int_t m_l = 0;
      for(const auto& T: hit->GetArrayTdcLeft()){
        HF1(E_VetoHid +100*(seg+1) +3, T);
        if(m_l < MaxDepth){
          event.e_vetolt[seg][m_l] = T;
          dst.e_vetolt[seg][m_l]   = T;
          ++m_l;
        }
        if(MinTdcE_Veto < T && T < MaxTdcE_Veto) is_hit_l = true;
      }
      if(is_hit_l) HF1(E_VetoHid+100*(seg+1)+5, Al);
      else         HF1(E_VetoHid+100*(seg+1)+7, Al);

      // Right
      Int_t Ar = hit->GetAdcRight();
      HF1(E_VetoHid+100*(seg+1)+2, Ar);
      event.e_vetora[seg] = Ar;
      dst.e_vetora[seg]   = Ar;

      Bool_t is_hit_r = false;
      Int_t m_r = 0;
      for(const auto& T: hit->GetArrayTdcRight()){
        HF1(E_VetoHid +100*(seg+1) +4, T);
        if(m_r < MaxDepth){
          event.e_vetort[seg][m_r] = T;
          dst.e_vetort[seg][m_r]   = T;
          ++m_r;
        }
        if(MinTdcT0 < T && T < MaxTdcT0)  is_hit_r = true;
      }
      if(is_hit_r) HF1(E_VetoHid+100*(seg+1)+6, Ar);
      else         HF1(E_VetoHid+100*(seg+1)+8, Ar);

      // Hitpat
      if(is_hit_l || is_hit_r){
        ++nh1;
      }
      if(is_hit_l && is_hit_r){
        ++nh2;
      }
    }
    HF1(E_VetoHid+2, nh1); HF1(E_VetoHid+4, nh2);
    event.e_vetonhits = e_veto_nhits;
  }


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

#if 1 // Normalized

#if 0 // BH2, SAC
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
            // if(event.sacnhits > 0) continue;
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

  // SAC
  hodoAna.DecodeHits("SAC");
  {
    Int_t nh=hodoAna.GetNHits("SAC");
    dst.nhSac = nh;
    HF1(SACHid+10, Double_t(nh));
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit("SAC", i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId()+1;
      HF1(SACHid+11, seg-0.5);
      Double_t a = hit->GetAUp();
      event.sacde[seg-1]  = a;
      Int_t n_mhit = hit->GetEntries();
      for(Int_t m=0; m<n_mhit; ++m){
        Double_t t = hit->GetTUp(m);
        Double_t ct = hit->GetCTUp(m);
        event.sacmt[i][m] = ct;
        HF1(SACHid+100*seg+11, t);
        HF1(SACHid+100*seg+12, a);
        HF1(SACHid+100*seg+13, ct);
      }
    }
  }
#endif // BH2, SAC

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
    Int_t nc = hodoAna.GetNClusters("SAC");
    dst.nhSac = nc;
    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster("SAC", i);
      if(!cl) continue;
      dst.csSac[i]  = cl->ClusterSize();
      dst.SacSeg[i] = cl->MeanSeg()+1;
      dst.tSac[i]   = cl->CMeanTime();
      dst.deSac[i]  = cl->DeltaE();
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

  // RF
  HB1(RFHid +0, "#Hits RF",        NumOfSegRF+1, 0., Double_t(NumOfSegRF+1));

  for(Int_t i=0; i<NumOfSegRF; ++i){
    TString title2 = Form("RF-%d Tdc", i);
    HB1(RFHid +100*(i+1) +3, title2, NbinTdcHr, MinTdcHr, MaxTdcHr);
  }

    //Tag-SF

  HB1(TagSFHid +0, "#Cluster Tag-SFF", NumOfSegTagSF+1,0.,Double_t(NumOfSegTagSF+1));
  HB1(TagSFHid +1, "#Cluster Tag-SFB", NumOfSegTagSF+1,0.,Double_t(NumOfSegTagSF+1));
  HB1(TagSFHid +10, "#ClusterSize Tag-SFF", NumOfSegTagSF+1,0.,Double_t(NumOfSegTagSF+1));
  HB1(TagSFHid +11, "#ClusterSize Tag-SFB", NumOfSegTagSF+1,0.,Double_t(NumOfSegTagSF+1));
  HB1(TagSFHid +20, "#Cluster Tag-SFF PL matching", NumOfSegTagSF+1,0.,Double_t(NumOfSegTagSF+1));
  HB1(TagSFHid +21, "#Cluster Tag-SFB PL matching", NumOfSegTagSF+1,0.,Double_t(NumOfSegTagSF+1));
  HB1(TagSFHid +30, "#Cluster Tag-SFF all matching", NumOfSegTagSF+1,0.,Double_t(NumOfSegTagSF+1));
  HB1(TagSFHid +31, "#Cluster Tag-SFB all matching", NumOfSegTagSF+1,0.,Double_t(NumOfSegTagSF+1));
  HB1(TagSFHid +3, "SFF leading", NbinTdc, MinTdc, MaxTdc);
  HB1(TagSFHid +4, "SFB leading", NbinTdc, MinTdc, MaxTdc);
  HB1(TagSFHid +13, "Tag-SFF time", 20, -10, 10);
  HB1(TagSFHid +14, "Tag-SFB time", 20, -10, 10);
  HB1(TagSFHid +5, "HitPat SFF", NumOfSegTagSF, 0, Double_t(NumOfSegTagSF));
  HB1(TagSFHid +6, "HitPat SFB", NumOfSegTagSF, 0, Double_t(NumOfSegTagSF));
  HB1(TagSFHid +15, "HitPat SFF cluster", NumOfSegTagSF*2, 0, Double_t(NumOfSegTagSF));
  HB1(TagSFHid +16, "HitPat SFB cluster", NumOfSegTagSF*2, 0, Double_t(NumOfSegTagSF));
  HB2(TagSFHid +17, "HitPat 2D SFF-SFB cluster",  NumOfSegTagSF*2, 0, Double_t(NumOfSegTagSF), NumOfSegTagSF*2, 0, Double_t(NumOfSegTagSF));
  //HB2(TagSFHid +17, "HitPat 2D SFF-SFB cluster",  NumOfSegTagSF, 0, Double_t(NumOfSegTagSF), NumOfSegTagSF, 0, Double_t(NumOfSegTagSF));
  HB2(TagSFHid +18, "HitPat 2D PL-SFF cluster",  NumOfSegTagPL, 0, Double_t(NumOfSegTagPL), NumOfSegTagSF*2, 0, Double_t(NumOfSegTagSF));
  //HB2(TagSFHid +18, "HitPat 2D PL-SFF cluster",  NumOfSegTagPL, 0, Double_t(NumOfSegTagPL), NumOfSegTagSF, 0, Double_t(NumOfSegTagSF));
  HB2(TagSFHid +19, "HitPat 2D PL-SFB cluster",  NumOfSegTagPL, 0, Double_t(NumOfSegTagPL), NumOfSegTagSF*2, 0, Double_t(NumOfSegTagSF));
  //HB2(TagSFHid +19, "HitPat 2D PL-SFB cluster",  NumOfSegTagPL, 0, Double_t(NumOfSegTagPL), NumOfSegTagSF, 0, Double_t(NumOfSegTagSF));
  HB1(TagSFHid +25, "HitPat SFF cluster PL matching", NumOfSegTagSF*2, 0, Double_t(NumOfSegTagSF));
  HB1(TagSFHid +26, "HitPat SFB cluster PL matching", NumOfSegTagSF*2, 0, Double_t(NumOfSegTagSF));
  HB2(TagSFHid +27, "HitPat 2D SFF-SFB cluster PL matching",  NumOfSegTagSF*2, 0, Double_t(NumOfSegTagSF), NumOfSegTagSF*2, 0, Double_t(NumOfSegTagSF));
  HB1(TagSFHid +35, "egamma-f", 40, 1.0, 3.0);
  HB1(TagSFHid +36, "egamma-b", 40, 1.0, 3.0);
  HB2(TagSFHid +37, "HitPat 2D SFF-SFB cluster final cand",  NumOfSegTagSF*2, 0, Double_t(NumOfSegTagSF), NumOfSegTagSF*2, 0, Double_t(NumOfSegTagSF));
  for(Int_t i=1;i<=NumOfSegTagSF;++i){
    TString title3 = Form("SFF-%d Tdc", i);
    TString title4 = Form("SFB-%d Tdc", i);
    HB1(TagSFHid +100*i +3, title3, NbinTdc, MinTdc, MaxTdc);
    HB1(TagSFHid +100*i +4, title4, NbinTdc, MinTdc, MaxTdc);
  }

  //Tag-PL
  HB1(TagPLHid +0, "#Hits Tag-PL w/o cut", NumOfSegTagPL+1,0.,Double_t(NumOfSegTagPL+1));
  HB1(TagPLHid +10, "#Hits Tag-PL", NumOfSegTagPL+1,0.,Double_t(NumOfSegTagPL+1));
  HB1(TagPLHid +3, "Tag-PL leading", NbinTdc, MinTdc, MaxTdc);
  HB1(TagPLHid +13, "Tag-PL time", 20, -10, 10);
  HB1(TagPLHid +5, "HitPat Tag-PL (w/o TDC cut)", NumOfSegTagPL, 0, Double_t(NumOfSegTagPL));
  HB1(TagPLHid +6, "HitPat Tag-PL (w/ TDC cut)", NumOfSegTagPL, 0, Double_t(NumOfSegTagPL));
  HB1(TagPLHid +16, "HitPat Tag-PL (w/ Time cut)", NumOfSegTagPL, 0, Double_t(NumOfSegTagPL));
  for(Int_t i=1;i<=NumOfSegTagPL;++i){
    TString title3 = Form("TagPL-%d Tdc", i);
    HB1(TagPLHid +100*i +3, title3, NbinTdc, MinTdc, MaxTdc);
  }

  // T0
  HB1(T0Hid +0, "#Hits T0",        NumOfSegT0+1, 0., Double_t(NumOfSegT0+1));
  HB1(T0Hid +2, "#Hits T0(Tor)",   NumOfSegT0+1, 0., Double_t(NumOfSegT0+1));
  HB1(T0Hid +4, "#Hits T0(Tand)",  NumOfSegT0+1, 0., Double_t(NumOfSegT0+1));

  for(Int_t i=0; i<NumOfSegT0; ++i){
    TString title1 = Form("T0-%d LeftAdc", i);
    TString title2 = Form("T0-%d RightAdc", i);
    TString title3 = Form("T0-%d LeftTdc", i);
    TString title4 = Form("T0-%d RightTdc", i);
    TString title5 = Form("T0-%d LeftAdc(w Tdc)", i);
    TString title6 = Form("T0-%d RightAdc(w Tdc)", i);
    TString title7 = Form("T0-%d LeftAdc(w/o Tdc)", i);
    TString title8 = Form("T0-%d RightAdc(w/o Tdc)", i);
    HB1(T0Hid +100*(i+1) +1, title1, NbinAdc, MinAdc, MaxAdc);
    HB1(T0Hid +100*(i+1) +2, title2, NbinAdc, MinAdc, MaxAdc);
    HB1(T0Hid +100*(i+1) +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(T0Hid +100*(i+1) +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(T0Hid +100*(i+1) +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(T0Hid +100*(i+1) +6, title6, NbinAdc, MinAdc, MaxAdc);
    HB1(T0Hid +100*(i+1) +7, title7, NbinAdc, MinAdc, MaxAdc);
    HB1(T0Hid +100*(i+1) +8, title8, NbinAdc, MinAdc, MaxAdc);
  }

#if 1
  // SAC
  HB1(SACHid +0, "#Hits SAC",        NumOfSegSAC+1, 0., Double_t(NumOfSegSAC+1));
  HB1(SACHid +1, "Hitpat SAC",       NumOfSegSAC,   0., Double_t(NumOfSegSAC));
  HB1(SACHid +2, "#Hits SAC(Tor)",   NumOfSegSAC+1, 0., Double_t(NumOfSegSAC+1));
  HB1(SACHid +3, "Hitpat SAC(Tor)",  NumOfSegSAC,   0., Double_t(NumOfSegSAC));
  HB1(SACHid +4, "#Hits SAC(Tand)",  NumOfSegSAC+1, 0., Double_t(NumOfSegSAC+1));
  HB1(SACHid +5, "Hitpat SAC(Tand)", NumOfSegSAC,   0., Double_t(NumOfSegSAC));
  for(Int_t i=1; i<=NumOfSegSAC; ++i){
    TString title1 = Form("SAC-%d UpAdc", i);
    TString title3 = Form("SAC-%d UpTdc", i);
    TString title5 = Form("SAC-%d UpAdc(w Tdc)", i);
    TString title7 = Form("SAC-%d UpAdc(w/o Tdc)", i);
    HB1(SACHid +100*i +1, title1, NbinTdc, MinTdc, MaxTdc);
    HB1(SACHid +100*i +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(SACHid +100*i +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(SACHid +100*i +7, title7, NbinAdc, MinAdc, MaxAdc);
  }
  HB1(SACHid +10, "#Hits SAC[Hodo]",     NumOfSegSAC+1, 0., Double_t(NumOfSegSAC+1));
  HB1(SACHid +11, "Hitpat SAC[Hodo]",    NumOfSegSAC,   0., Double_t(NumOfSegSAC));
  for(Int_t i=1; i<=NumOfSegSAC; ++i){
    TString title1 = Form("SAC-%d Time", i);
    TString title3 = Form("SAC-%d dE", i);
    TString title5 = Form("SAC-%d CTime", i);
    HB1(SACHid +100*i +11, title1, 500, -5., 45.);
    HB1(SACHid +100*i +12, title3, 200, -0.5, 4.5);
    HB1(SACHid +100*i +13, title5, 500, -5., 45.);
  }
#endif // BH2, SAC

  // E_Veto
  HB1(E_VetoHid +0, "#Hits E_Veto",        NumOfSegE_Veto+1, 0., Double_t(NumOfSegE_Veto+1));
  HB1(E_VetoHid +2, "#Hits E_Veto(Tor)",   NumOfSegE_Veto+1, 0., Double_t(NumOfSegE_Veto+1));
  HB1(E_VetoHid +4, "#Hits E_Veto(Tand)",  NumOfSegE_Veto+1, 0., Double_t(NumOfSegE_Veto+1));

  for(Int_t i=0; i<NumOfSegE_Veto; ++i){
    TString title1 = Form("E_Veto-%d LeftAdc", i);
    TString title2 = Form("E_Veto-%d RightAdc", i);
    TString title3 = Form("E_Veto-%d LeftTdc", i);
    TString title4 = Form("E_Veto-%d RightTdc", i);
    TString title5 = Form("E_Veto-%d LeftAdc(w Tdc)", i);
    TString title6 = Form("E_Veto-%d RightAdc(w Tdc)", i);
    TString title7 = Form("E_Veto-%d LeftAdc(w/o Tdc)", i);
    TString title8 = Form("E_Veto-%d RightAdc(w/o Tdc)", i);
    HB1(E_VetoHid +100*(i+1) +1, title1, NbinAdc, MinAdc, MaxAdc);
    HB1(E_VetoHid +100*(i+1) +2, title2, NbinAdc, MinAdc, MaxAdc);
    HB1(E_VetoHid +100*(i+1) +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(E_VetoHid +100*(i+1) +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(E_VetoHid +100*(i+1) +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(E_VetoHid +100*(i+1) +6, title6, NbinAdc, MinAdc, MaxAdc);
    HB1(E_VetoHid +100*(i+1) +7, title7, NbinAdc, MinAdc, MaxAdc);
    HB1(E_VetoHid +100*(i+1) +8, title8, NbinAdc, MinAdc, MaxAdc);
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
    HB1(TOFHid +100*(i+1) +1, title1, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*(i+1) +2, title2, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*(i+1) +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(TOFHid +100*(i+1) +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(TOFHid +100*(i+1) +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*(i+1) +6, title6, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*(i+1) +7, title7, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*(i+1) +8, title8, NbinAdc, MinAdc, MaxAdc);
  }

#if 1 // TOF Normalized
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

  //RF
  tree->Branch("rfnhits",   &event.rfnhits,   "rfnhits/I");
  tree->Branch("rft",        event.rft,       Form("rft[%d][%d]/D", NumOfSegRF, MaxDepth));
  //TagSF
  tree->Branch("tagsffnhits", &event.tagsffnhits, "tagsffnhits/I");
  tree->Branch("tagsffhitpat", event.tagsffhitpat, Form("tagsffhitpat[%d]/I",NumOfSegTagSF));
  tree->Branch("tagsfft",event.tagsfft, Form("tagsfft[%d][%d]/D",NumOfSegTagSF,MaxDepth));
  tree->Branch("tagsfbnhits", &event.tagsfbnhits, "tagsfbnhits/I");
  tree->Branch("tagsfbhitpat", event.tagsfbhitpat, Form("tagsfbhitpat[%d]/I",NumOfSegTagSF));
  tree->Branch("tagsfbt",event.tagsfbt, Form("tagsfbt[%d][%d]/D",NumOfSegTagSF,MaxDepth));
  //TagPL
  tree->Branch("tagplnhits", &event.tagplnhits, "tagplnhits/I");
  tree->Branch("tagplhitpat", event.tagplhitpat, Form("tagplhitpat[%d]/I",NumOfSegTagPL));
  tree->Branch("tagpla", event.tagpla, Form("tagpla[%d]/D",NumOfSegTagPL));
  tree->Branch("tagplt",event.tagplt, Form("tagplt[%d][%d]/D",NumOfSegTagPL,MaxDepth));

  //T0
  tree->Branch("t0nhits",   &event.t0nhits,   "t0nhits/I");
  tree->Branch("t0la",       event.t0la,      Form("t0la[%d]/D", NumOfSegT0));
  tree->Branch("t0lt",       event.t0lt,      Form("t0lt[%d][%d]/D", NumOfSegT0, MaxDepth));
  tree->Branch("t0ra",       event.t0ra,      Form("t0ra[%d]/D", NumOfSegT0));
  tree->Branch("t0rt",       event.t0rt,      Form("t0rt[%d][%d]/D", NumOfSegT0, MaxDepth));
  // //SAC
  tree->Branch("sacnhits",   &event.sacnhits,   "sacnhits/I");
  tree->Branch("sachitpat",   event.sachitpat,  Form("sachitpat[%d]/I", NumOfSegSAC));
  tree->Branch("saca",        event.saca,       Form("saca[%d]/D", NumOfSegSAC));
  tree->Branch("sact",        event.sact,       Form("sact[%d][%d]/D", NumOfSegSAC, MaxDepth));
  tree->Branch("sac1nhits",   &event.sac1nhits,   "sac1nhits/I");
  tree->Branch("sac2nhits",   &event.sac2nhits,   "sac2nhits/I");
  //E_Veto
  tree->Branch("e_vetonhits",   &event.e_vetonhits,   "e_vetonhits/I");
  tree->Branch("e_vetola",       event.e_vetola,      Form("e_vetola[%d]/D", NumOfSegE_Veto));
  tree->Branch("e_vetolt",       event.e_vetolt,      Form("e_vetolt[%d][%d]/D", NumOfSegE_Veto, MaxDepth));
  tree->Branch("e_vetora",       event.e_vetora,      Form("e_vetora[%d]/D", NumOfSegE_Veto));
  tree->Branch("e_vetort",       event.e_vetort,      Form("e_vetort[%d][%d]/D", NumOfSegE_Veto, MaxDepth));
  //TOF
  tree->Branch("tofnhits",   &event.tofnhits,   "tofnhits/I");
  tree->Branch("tofhitpat",   event.tofhitpat,  Form("tofhitpat[%d]/I", NumOfSegTOF));
  tree->Branch("tofua",       event.tofua,      Form("tofua[%d]/D", NumOfSegTOF));
  tree->Branch("tofut",       event.tofut,      Form("tofut[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofda",       event.tofda,      Form("tofda[%d]/D", NumOfSegTOF));
  tree->Branch("tofdt",       event.tofdt,      Form("tofdt[%d][%d]/D", NumOfSegTOF, MaxDepth));

#if 1 // Normalized, Dst
  /*
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
  tree->Branch("sacmt",     event.sacmt,     Form("sacmt[%d][%d]/D", NumOfSegSAC, MaxDepth));
  tree->Branch("sacde",     event.sacde,     Form("sacde[%d]/D", NumOfSegSAC));
  */
  tree->Branch("tofmt",     event.tofmt,     Form("tofmt[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofde",     event.tofde,     Form("tofde[%d]/D", NumOfSegTOF));
  tree->Branch("tofude",     event.tofude,     Form("tofude[%d]/D", NumOfSegTOF));
  tree->Branch("tofdde",     event.tofdde,     Form("tofdde[%d]/D", NumOfSegTOF));
  tree->Branch("tofctu",     event.tofctu,     Form("tofctu[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofctd",     event.tofctd,     Form("tofctd[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofcmt",     event.tofcmt,     Form("tofcmt[%d][%d]/D", NumOfSegTOF, MaxDepth));
  /*
  tree->Branch("t0",        event.t0,        Form("t0[%d][%d]/D",  NumOfSegBH2, MaxDepth));
  tree->Branch("ct0",       event.ct0,       Form("ct0[%d][%d]/D", NumOfSegBH2, MaxDepth));

  tree->Branch("Time0Seg", &event.Time0Seg,  "Time0Seg/D");
  tree->Branch("deTime0",  &event.deTime0,   "deTime0/D");
  tree->Branch("Time0",    &event.Time0,     "Time0/D");
  tree->Branch("CTime0",   &event.CTime0,    "CTime0/D");
  */
#endif

#if 0
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

  hodo->Branch("nhSac",     &dst.nhSac,     "nhSac/I");
  hodo->Branch("SacSeg",     dst.SacSeg,    "SacSeg[nhSac]/D");
  hodo->Branch("tSac",       dst.tSac,      "tSac[nhSac]/D");
  hodo->Branch("deSac",      dst.deSac,     "deSac[nhSac]/D");

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
     InitializeParameter<TAGPLMatch>("TAGPLMTH") &&
#if SDCoutreq
     InitializeParameter<DCDriftParamMan>("DCDRFT") &&
     InitializeParameter<DCTdcCalibMan>("DCTDC")    &&
#endif
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
