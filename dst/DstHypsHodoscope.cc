// -*- C++ -*-

#include "DstHelper.hh"

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include <TMath.h>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DCGeomMan.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "HodoPHCMan.hh"

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const Bool_t USE_M2 = false;
const Bool_t USE_XYCut = false;
TString ClassName() { return TString("DstHypsHodoscope"); }
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kHypsTracking, kHodoscope, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[HypsTracking]",
  "[Hodoscope]", "[OutFile]" };
std::vector<TString> TreeName =
{ "", "", "hyps", "hodo", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________
struct Event
{
  Int_t status;

  // SdcOut
  Int_t ntSdcOut;
  Double_t chisqrSdcOut[MaxHits];
  Double_t x0SdcOut[MaxHits];
  Double_t y0SdcOut[MaxHits];
  Double_t u0SdcOut[MaxHits];
  Double_t v0SdcOut[MaxHits];

  // HypsTracking
  Int_t ntHyps;
  Double_t path[MaxHits];
  Double_t pHyps[MaxHits];
  Double_t qHyps[MaxHits];
  Double_t chisqrHyps[MaxHits];
  Double_t thetaHyps[MaxHits];
  Double_t xtgtHyps[MaxHits];
  Double_t ytgtHyps[MaxHits];
  Double_t utgtHyps[MaxHits];
  Double_t vtgtHyps[MaxHits];
  Double_t xtofHyps[MaxHits];
  Double_t ytofHyps[MaxHits];
  Double_t utofHyps[MaxHits];
  Double_t vtofHyps[MaxHits];
  Double_t lxtofHyps[MaxHits];
  Double_t lytofHyps[MaxHits];
  Double_t lutofHyps[MaxHits];
  Double_t lvtofHyps[MaxHits];
  Double_t tofsegHyps[MaxHits];
  Double_t vpx[NumOfLayersVP*MaxHits];
  Double_t vpy[NumOfLayersVP*MaxHits];
  Double_t vpu[NumOfLayersVP*MaxHits];
  Double_t vpv[NumOfLayersVP*MaxHits];

  // Hodoscope
  Int_t trigflag[NumOfSegTrig];
  Int_t trigpat[NumOfSegTrig];

  Int_t nhBh1;
  Int_t csBh1[NumOfSegBH1*MaxDepth];
  Double_t Bh1Seg[NumOfSegBH1*MaxDepth];
  Double_t tBh1[NumOfSegBH1*MaxDepth];
  Double_t dtBh1[NumOfSegBH1*MaxDepth];
  Double_t deBh1[NumOfSegBH1*MaxDepth];

  Int_t nhBh2;
  Int_t csBh2[NumOfSegBH2*MaxDepth];
  Double_t Bh2Seg[NumOfSegBH2*MaxDepth];
  Double_t tBh2[NumOfSegBH2*MaxDepth];
  Double_t dtBh2[NumOfSegBH2*MaxDepth];
  Double_t t0Bh2[NumOfSegBH2*MaxDepth];
  Double_t deBh2[NumOfSegBH2*MaxDepth];

  Double_t btof[NumOfSegBH1*MaxDepth];
  Double_t Time0Seg;

  Int_t nhTof;
  Int_t csTof[NumOfSegTOF];
  Double_t TofSeg[NumOfSegTOF];
  Double_t tTof[NumOfSegTOF];
  Double_t dtTof[NumOfSegTOF];
  Double_t deTof[NumOfSegTOF];

  Int_t nhAc1;
  Int_t csAc1[NumOfSegAC1];
  Double_t Ac1Seg[NumOfSegAC1];
  Double_t tAc1[NumOfSegAC1];

  // HypsHodoscope
  Int_t    m2Combi;
  Double_t beta[MaxHits];
  Double_t stof[MaxHits];
  Double_t cstof[MaxHits];
  Double_t m2[MaxHits];

  Double_t tofua[NumOfSegTOF];
  Double_t tofda[NumOfSegTOF];
  enum eParticle { Pion, Kaon, Proton, nParticle };
  Double_t tTofCalc[nParticle];
  Double_t utTofSeg[NumOfSegTOF][MaxDepth];
  Double_t dtTofSeg[NumOfSegTOF][MaxDepth];
  Double_t udeTofSeg[NumOfSegTOF];
  Double_t ddeTofSeg[NumOfSegTOF];

};

//_____________________________________________________________________
struct Src
{
  Int_t trigflag[NumOfSegTrig];
  Int_t trigpat[NumOfSegTrig];

  Int_t ntSdcOut;
  Double_t chisqrSdcOut[MaxHits];
  Double_t x0SdcOut[MaxHits];
  Double_t y0SdcOut[MaxHits];
  Double_t u0SdcOut[MaxHits];
  Double_t v0SdcOut[MaxHits];

  Int_t    ntHyps;
  Double_t path[MaxHits];
  Double_t pHyps[MaxHits];
  Double_t qHyps[MaxHits];
  Double_t chisqrHyps[MaxHits];
  Double_t thetaHyps[MaxHits];
  Double_t xtgtHyps[MaxHits];
  Double_t ytgtHyps[MaxHits];
  Double_t utgtHyps[MaxHits];
  Double_t vtgtHyps[MaxHits];
  Double_t xtofHyps[MaxHits];
  Double_t ytofHyps[MaxHits];
  Double_t utofHyps[MaxHits];
  Double_t vtofHyps[MaxHits];
  Double_t lxtofHyps[MaxHits];
  Double_t lytofHyps[MaxHits];
  Double_t lutofHyps[MaxHits];
  Double_t lvtofHyps[MaxHits];
  Double_t tofsegHyps[MaxHits];
  Double_t vpx[NumOfLayersVP];
  Double_t vpy[NumOfLayersVP];
  Double_t vpu[NumOfLayersVP];
  Double_t vpv[NumOfLayersVP];

  Int_t    nhBh1;
  Int_t    csBh1[NumOfSegBH1*MaxDepth];
  Double_t Bh1Seg[NumOfSegBH1*MaxDepth];
  Double_t tBh1[NumOfSegBH1*MaxDepth];
  Double_t dtBh1[NumOfSegBH1*MaxDepth];
  Double_t deBh1[NumOfSegBH1*MaxDepth];

  Int_t    nhBh2;
  Int_t    csBh2[NumOfSegBH2*MaxDepth];
  Double_t Bh2Seg[NumOfSegBH2*MaxDepth];
  Double_t tBh2[NumOfSegBH2*MaxDepth];
  Double_t dtBh2[NumOfSegBH2*MaxDepth];
  Double_t t0Bh2[NumOfSegBH2*MaxDepth];
  Double_t deBh2[NumOfSegBH2*MaxDepth];

  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  Int_t    nhTof;
  Int_t    csTof[NumOfSegTOF];
  Double_t TofSeg[NumOfSegTOF];
  Double_t tTof[NumOfSegTOF];
  Double_t dtTof[NumOfSegTOF];
  Double_t deTof[NumOfSegTOF];

  Int_t nhAc1;
  Int_t csAc1[NumOfSegAC1];
  Double_t Ac1Seg[NumOfSegAC1];
  Double_t tAc1[NumOfSegAC1];

  ////////// for HodoParam
  Double_t tofua[NumOfSegTOF];
  Double_t tofda[NumOfSegTOF];
  Double_t utTofSeg[NumOfSegTOF][MaxDepth];
  Double_t dtTofSeg[NumOfSegTOF][MaxDepth];
  Double_t udeTofSeg[NumOfSegTOF];
  Double_t ddeTofSeg[NumOfSegTOF];

};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
}

//_____________________________________________________________________
int
main(int argc, char **argv)
{
  std::vector<std::string> arg(argv, argv+argc);

  if(!CheckArg(arg))
    return EXIT_FAILURE;
  if(!DstOpen(arg))
    return EXIT_FAILURE;
  if(!gConf.Initialize(arg[kConfFile]))
    return EXIT_FAILURE;
  if(!gConf.InitializeUnpacker())
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries(TTreeCont);
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  Int_t ievent = skip;
  for(; ievent<nevent && !CatchSignal::Stop(); ++ievent){
    gCounter.check();
    InitializeEvent();
    if(DstRead(ievent)) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
	    << ievent << std::endl;

  DstClose();

  return EXIT_SUCCESS;
}

//_____________________________________________________________________
bool
dst::InitializeEvent()
{
  event.status   = 0;
  event.ntSdcOut = 0;
  event.ntHyps = 0;
  event.nhBh1    = 0;
  event.nhBh2    = 0;
  event.nhTof    = 0;
  event.nhAc1    = 0;
  event.m2Combi  = 0;
  event.Time0Seg = qnan;

  for(Int_t i=0; i<NumOfSegTrig; ++i){
    event.trigflag[i] = -1;
    event.trigpat[i]  = -1;
  }

  for(Int_t i=0; i<Event::nParticle; ++i){
    event.tTofCalc[i] = qnan;
  }

  for(Int_t i=0;i<MaxHits;++i){
    event.chisqrSdcOut[i] = qnan;
    event.x0SdcOut[i]     = qnan;
    event.y0SdcOut[i]     = qnan;
    event.u0SdcOut[i]     = qnan;
    event.v0SdcOut[i]     = qnan;
    event.path[i]         = qnan;
    event.pHyps[i]      = qnan;
    event.qHyps[i]      = qnan;
    event.chisqrHyps[i] = qnan;
    event.thetaHyps[i] = qnan;
    event.xtgtHyps[i]  = qnan;
    event.ytgtHyps[i]  = qnan;
    event.utgtHyps[i]  = qnan;
    event.vtgtHyps[i]  = qnan;
    event.xtofHyps[i]  = qnan;
    event.ytofHyps[i]  = qnan;
    event.utofHyps[i]  = qnan;
    event.vtofHyps[i]  = qnan;
    event.lxtofHyps[i]  = qnan;
    event.lytofHyps[i]  = qnan;
    event.lutofHyps[i]  = qnan;
    event.lvtofHyps[i]  = qnan;
    event.tofsegHyps[i]  = qnan;
  }

  for (Int_t l = 0; l < NumOfLayersVP; ++l) {
    event.vpx[l] = qnan;
    event.vpy[l] = qnan;
    event.vpu[l] = qnan;
    event.vpv[l] = qnan;
  }

  for(Int_t i=0;i<NumOfSegBH1*MaxDepth;++i){
    event.Bh1Seg[i] = -1;
    event.csBh1[i]  = 0;
    event.tBh1[i]   = qnan;
    event.dtBh1[i]  = qnan;
    event.deBh1[i]  = qnan;
    event.btof[i]   = qnan;
  }

  for(Int_t i=0;i<NumOfSegBH2*MaxDepth;++i){
    event.Bh2Seg[i] = -1;
    event.csBh2[i]  = 0;
    event.tBh2[i]   = qnan;
    event.dtBh2[i]  = qnan;
    event.t0Bh2[i]  = qnan;
    event.deBh2[i]  = qnan;
  }

  for(Int_t i=0;i<NumOfSegTOF;++i){
    event.TofSeg[i] = -1;
    event.csTof[i]  = 0;
    event.tTof[i]   = qnan;
    event.dtTof[i]  = qnan;
    event.deTof[i]  = qnan;
    event.tofua[i] = qnan;
    event.tofda[i] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      event.utTofSeg[i][m] = qnan;
      event.dtTofSeg[i][m] = qnan;
    }
    event.udeTofSeg[i] = qnan;
    event.ddeTofSeg[i] = qnan;
  }

  for(Int_t i=0;i<NumOfSegAC1;++i){
    event.Ac1Seg[i] = -1;
    event.csAc1[i]  = 0;
    event.tAc1[i]   = qnan;
  }

  for(Int_t i=0; i<Event::nParticle; ++i){
    event.tTofCalc[i] = qnan;
  }

  for(Int_t i=0; i<MaxHits; ++i){
    event.beta[i]  = qnan;
    event.stof[i]  = qnan;
    event.cstof[i] = qnan;
    event.m2[i]    = qnan;
  }
  return true;
}

//_____________________________________________________________________
bool
dst::DstOpen(std::vector<std::string> arg)
{
  Int_t open_file = 0;
  Int_t open_tree = 0;
  for(std::size_t i=0; i<nArgc; ++i){
    if(i==kProcess || i==kConfFile || i==kOutFile) continue;
    open_file += OpenFile(TFileCont[i], arg[i]);
    open_tree += OpenTree(TFileCont[i], TTreeCont[i], TreeName[i]);
  }

  if(open_file!=open_tree || open_file!=nArgc-3)
    return false;
  if(!CheckEntries(TTreeCont))
    return false;

  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");

  return true;
}

//_____________________________________________________________________
bool
dst::DstRead(Int_t ievent)
{
  static const auto StofOffset = gUser.GetParameter("StofOffset");

  if(ievent%10000==0){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  HF1(1, event.status++);

  event.ntSdcOut = src.ntSdcOut;
  event.ntHyps = src.ntHyps;
  event.nhBh1    = src.nhBh1;
  event.nhBh2    = src.nhBh2;
  event.nhTof    = src.nhTof;
  event.nhAc1    = src.nhAc1;

#if 0
  std::cout<<"[event]: "<<std::setw(6)<<ievent<<" ";
  std::cout<<"[ntHyps]: "<<std::setw(2)<<src.ntHyps<<" ";
  std::cout<<"[nhBh1]: "<<std::setw(2)<<src.nhBh1<<" ";
  std::cout<<"[nhBh2]: "<<std::setw(2)<<src.nhBh2<<" ";
  std::cout<<"[nhTof]: "<<std::setw(2)<<src.nhTof<<" "<<std::endl;
  std::cout<<"[nhAc1]: "<<std::setw(2)<<src.nhAc1<<" "<<std::endl;
#endif

  for(Int_t i=0;i<NumOfSegTrig;++i){
    Int_t tdc = src.trigflag[i];
    if(tdc<=0) continue;
    event.trigpat[i]  = i + 1;
    event.trigflag[i] = tdc;
  }

  // if(event.nhBh1<=0) return true;
  HF1(1, event.status++);

  // if(event.nhBh2<=0) return true;
  HF1(1, event.status++);

  // if(event.nhTof<=0) return true;
  HF1(1, event.status++);

  // if(event.ntHyps<=0) return true;
  // if(event.ntHyps>MaxHits)
  //   event.ntHyps = MaxHits;

  HF1(1, event.status++);

  Double_t time0 = src.CTime0;
  for(Int_t i=0; i<src.nhBh2; ++i){
    event.Bh2Seg[i] = src.Bh2Seg[i];
    event.tBh2[i]   = src.tBh2[i];
    event.t0Bh2[i]  = src.t0Bh2[i];
    event.deBh2[i]  = src.deBh2[i];
  }
  event.Time0Seg = src.Time0Seg;

  ////////// for BeamTof
  // Double_t btof = qnan;
  for(Int_t i=0; i<src.nhBh1; ++i){
    event.Bh1Seg[i] = src.Bh1Seg[i];
    event.tBh1[i]   = src.tBh1[i];
    event.deBh1[i]  = src.deBh1[i];
    event.btof[i]   = src.tBh1[i] - time0;
    // if(i==0) btof = src.tBh1[i] - time0;
  }

  for(Int_t it=0; it<src.ntSdcOut; ++it){
    event.chisqrSdcOut[it] = src.chisqrSdcOut[it];
    event.x0SdcOut[it] = src.x0SdcOut[it];
    event.y0SdcOut[it] = src.y0SdcOut[it];
    event.u0SdcOut[it] = src.u0SdcOut[it];
    event.v0SdcOut[it] = src.v0SdcOut[it];
  }

  for(Int_t it=0; it<src.nhTof; ++it){
    event.csTof[it]  = src.csTof[it];
    event.TofSeg[it] = src.TofSeg[it];
    event.tTof[it]   = src.tTof[it];
    event.dtTof[it]  = src.dtTof[it];
    event.deTof[it]  = src.deTof[it];
  }
  for(Int_t it=0; it<NumOfSegTOF; ++it){
    event.tofua[it] = src.tofua[it];
    event.tofda[it] = src.tofda[it];
    for(Int_t m=0; m<MaxDepth; ++m){
      event.utTofSeg[it][m] = src.utTofSeg[it][m];
      event.dtTofSeg[it][m] = src.dtTofSeg[it][m];
    }
    event.udeTofSeg[it] = src.udeTofSeg[it];
    event.ddeTofSeg[it] = src.ddeTofSeg[it];
    // TOF ADC-Pedestal
    Bool_t is_hit = false;
    for(Int_t i=0; i<event.nhTof; ++i){
      if(it == event.TofSeg[i]-1) is_hit = true;
    }
    if(!is_hit){
      HF1(30000+(it+1)*100+1,  event.tofua[it]);
      HF1(30000+(it+1)*100+11, event.tofda[it]);
    }
  }

  for(Int_t it=0; it<src.nhAc1; ++it){
    event.csAc1[it]  = src.csAc1[it];
    event.Ac1Seg[it] = src.Ac1Seg[it];
    event.tAc1[it]   = src.tAc1[it];
  }

  Int_t m2Combi = event.nhTof*event.ntHyps;
  if(m2Combi>MaxHits || m2Combi<0){
    std::cout << FUNC_NAME << " too much m2Combi : " << m2Combi << std::endl;
    return false;
  }
  event.m2Combi = m2Combi;

  HF1(1, event.status++);
  Int_t mm=0;
  for(Int_t it=0; it<src.ntHyps; ++it){
    event.path[it] = src.path[it];
    event.pHyps[it] = src.pHyps[it];
    event.qHyps[it] = src.qHyps[it];
    event.chisqrHyps[it] = src.chisqrHyps[it];
    event.thetaHyps[it]  = src.thetaHyps[it];
    event.xtgtHyps[it] = src.xtgtHyps[it];
    event.ytgtHyps[it] = src.ytgtHyps[it];
    event.utgtHyps[it] = src.utgtHyps[it];
    event.vtgtHyps[it] = src.vtgtHyps[it];
    event.xtofHyps[it] = src.xtofHyps[it];
    event.ytofHyps[it] = src.ytofHyps[it];
    event.utofHyps[it] = src.utofHyps[it];
    event.vtofHyps[it] = src.vtofHyps[it];
    event.lxtofHyps[it] = src.lxtofHyps[it];
    event.lytofHyps[it] = src.lytofHyps[it];
    event.lutofHyps[it] = src.lutofHyps[it];
    event.lvtofHyps[it] = src.lvtofHyps[it];
    event.tofsegHyps[it] = src.tofsegHyps[it];
    Double_t xtgt  = event.xtgtHyps[it];
    Double_t ytgt  = event.ytgtHyps[it];
    Double_t lytof = event.lytofHyps[it];
    Double_t pHyps = event.pHyps[it];
    Double_t qHyps = event.qHyps[it];
    if(event.chisqrHyps[it] < 200){
      HF1(10, pHyps);
      HF1(13, event.path[it]);
      HF1(50, event.chisqrHyps[it]);
      HF1(51, event.xtofHyps[it]);
      HF1(52, event.ytofHyps[it]);
      HF1(53, event.utofHyps[it]);
      HF1(54, event.vtofHyps[it]);
      HF1(55, event.lxtofHyps[it]);
      HF1(56, event.lytofHyps[it]);
      HF1(57, event.lutofHyps[it]);
      HF1(58, event.lvtofHyps[it]);
      HF1(59, event.tofsegHyps[it]);
    }
    if(src.ntHyps == 1){
      for (Int_t l = 0; l < NumOfLayersVP; ++l) {
	event.vpx[l] = src.vpx[l];
	event.vpy[l] = src.vpy[l];
	event.vpu[l] = src.vpu[l];
	event.vpv[l] = src.vpv[l];
      }
    }
    if(it == 0){
      event.tTofCalc[Event::Pion] =
        Kinematics::CalcTimeOfFlight(event.pHyps[it],
                                     event.path[it],
                                     pdg::PionMass());
      event.tTofCalc[Event::Kaon] =
        Kinematics::CalcTimeOfFlight(event.pHyps[it],
                                     event.path[it],
                                     pdg::KaonMass());
      event.tTofCalc[Event::Proton] =
        Kinematics::CalcTimeOfFlight(event.pHyps[it],
                                     event.path[it],
                                     pdg::ProtonMass());
    }

    for(Int_t itof=0; itof<src.nhTof; ++itof){
      Int_t tofseg = (Int_t)event.TofSeg[itof];
      if( tofseg != event.tofsegHyps[0] ) continue;

      ////////// TimeCut
      Double_t stof = event.tTof[itof] - time0 + StofOffset;
      Double_t cstof = stof;
      // gPHC.DoStofCorrection(8, 0, tofseg-1, 2, stof, btof, cstof);
      Double_t beta = event.path[it]/cstof/MathTools::C();
      event.beta[mm] = beta;
      event.stof[mm] = stof;// - event.tTofCalc[0];
      event.cstof[mm] = cstof;
      Double_t m2 = Kinematics::MassSquare(pHyps, event.path[it], cstof);
      event.m2[mm] = m2;
      // Bool_t is_pion = (// qHyps > 0
      //                   // &&
      //                   TMath::Abs(m2-pdg::PionMass()*pdg::PionMass()) < 0.1);
      // Bool_t is_kaon = (qHyps > 0
      //                   && TMath::Abs(m2-pdg::KaonMass()*pdg::KaonMass()) < 0.12);
#if 0
      std::cout << "#D DebugPrint() Event : " << ievent << std::endl
		<< std::setprecision(3) << std::fixed
		<< "   time0   : " << time0 << std::endl
		<< "   offset  : " << StofOffset << std::endl
		<< "   tTof    : " << event.tTof[itof] << std::endl
		<< "   stof    : " << stof << std::endl
		<< "   pHyps    : " << pHyps << std::endl
		<< "   m2      : " << m2 << std::endl;
#endif

      if(event.chisqrHyps[it] < 200.){
        for(Int_t ip=0; ip<Event::nParticle; ++ip){
          if(USE_M2){
            if(ip == Event::Pion && TMath::Abs(m2-0.0194) > 0.1) continue;
            if(ip == Event::Proton && TMath::Abs(m2-0.88) > 0.2) continue;
          }
          HF1(33000+tofseg*100+ip+1,  event.tofua[tofseg-1]);
          HF1(33000+tofseg*100+ip+11, event.tofda[tofseg-1]);
        }
        if(TMath::Abs(lytof) < 5.0){
          HF1(20000+tofseg*100+1, event.dtTof[itof]);
          HF2(20000+1, event.dtTof[itof], tofseg-1);
        }

	Bool_t xy_ok = USE_XYCut
	  ? (TMath::Abs(xtgt) < 25. && TMath::Abs(ytgt) < 20.)
	  : true;
	if(xy_ok){
	  if( it == 0 ){
	    for(Int_t ip=0; ip<Event::nParticle; ++ip){
	      if(USE_M2){
		if(ip == Event::Pion && m2 > 0.2) continue;
		if(ip == Event::Kaon) continue;
		if(ip == Event::Proton && (m2 < 0.5 || !xy_ok)) continue;
	      }
	      else{
		Bool_t flagAc1 = false;
		for( Int_t it=0, nhit=event.nhAc1; it<nhit; it++ ){
		  if( event.Ac1Seg[it] == 21 ){ flagAc1=true; break; }
		}
		if(ip == Event::Pion && !flagAc1) continue; // select Pion
		if(ip == Event::Kaon) continue;
		if(ip == Event::Proton && (m2 < 0.5 || !xy_ok)) continue;
	      }
	      HF2(10000+ip+1, cstof-event.tTofCalc[ip], tofseg-1);
	      HF1(10000+tofseg*100+ip+1, cstof-event.tTofCalc[ip]);
	      if(TMath::Abs(event.dtTof[itof] < 0.1)){
		for(Int_t mh=0; mh<MaxDepth; ++mh){
		  HF2(40000+tofseg*100+ip+1,
		      event.udeTofSeg[tofseg-1],
		      event.tTofCalc[ip] - stof);
		  HF2(40000+tofseg*100+ip+11,
		      event.ddeTofSeg[tofseg-1],
		      event.tTofCalc[ip] - stof);
		  HF2(40000+tofseg*100,
		      event.udeTofSeg[tofseg-1],
		      event.tTofCalc[ip] - stof);
		  HF2(40000+tofseg*100+10,
		      event.ddeTofSeg[tofseg-1],
		      event.tTofCalc[ip] - stof);
		}
	      }
	    }
	  }
	  HF1(11, qHyps*m2);
	  HF1(12, beta);
	  HF1(14, stof);
	  HF2(20, qHyps*m2, pHyps);
	  if(!TMath::IsNaN(event.Time0Seg)){
	    HF1(1100+(Int_t)event.Time0Seg, qHyps*m2);
	    HF2(1200+(Int_t)event.Time0Seg, qHyps*m2, pHyps);
	  }
	  HF1(2100+tofseg, qHyps*m2);
	  HF2(2200+tofseg, qHyps*m2, pHyps);
	}
      }
      ++mm;
    }
  }

  HF1(1, event.status++);

  return true;
}

//_____________________________________________________________________
bool
dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  for(Int_t i=0, n=TFileCont.size(); i<n; ++i){
    if(TTreeCont[i]) delete TTreeCont[i];
    if(TFileCont[i]) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeHistograms()
{
  const Int_t NBinP = 220;
  const Double_t MinP =   0;
  const Double_t MaxP = 2.2;
  const Int_t NBinM2 = 340;
  const Double_t MinM2 = -1.4;
  const Double_t MaxM2 =  2.0;

  HB1(1, "Status", 21, 0., 21.);

  TString name[Event::nParticle] = { "Pion", "Kaon", "Proton" };

  HB1(10, "pHyps",    NBinP,  MinP, MaxP);
  HB1(11, "ChargexMassSquare", NBinM2, MinM2, MaxM2);
  HB1(12, "beta", 500, 0., 1.);
  HB1(13, "path", 500, 3000., 8000.);
  HB1(14, "stof", 500, 0., 100.);
  HB2(20, "pHyps % ChargexMassSquare", NBinM2, MinM2, MaxM2, NBinP, MinP, MaxP);

  HB1(50, "Chisqr HypsTrack", 500, 0., 50.);
  HB1(51, "Xtof HypsTrack", 200, -100., 100.);
  HB1(52, "Ytof HypsTrack", 400, -200., 200.);
  HB1(53, "Utof HypsTrack", 300, -0.30, 0.30);
  HB1(54, "Vtof HypsTrack", 300, -0.30, 0.30);
  HB1(55, "Local Xtof HypsTrack", 700, -700., 700.);
  HB1(56, "Local Ytof HypsTrack", 350, -350., 350.);
  HB1(57, "Local Utof HypsTrack", 300, -0.30, 0.30);
  HB1(58, "Local Vtof HypsTrack", 300, -0.30, 0.30);
  HB1(59, "tofseg HypsTrack", 21, 0, 21);

  for(Int_t i=0;i<NumOfSegBH2;++i){
    HB1(1100+i+1, Form("ChargexMassSquare [BH2-%d]", i+1),
        NBinM2, MinM2, MaxM2);
    HB2(1200+i+1, Form("pHyps %% ChargexMassSquare [BH2-%d]", i+1),
        NBinM2, MinM2, MaxM2, NBinP, MinP, MaxP);
  }
  for(Int_t i=0;i<NumOfSegTOF;++i){
    HB1(2100+i+1, Form("ChargexMassSquare [TOF-%d]", i+1),
        NBinM2, MinM2, MaxM2);
    HB2(2200+i+1, Form("pHyps %% ChargexMassSquare [TOF-%d]", i+1),
        NBinM2, MinM2, MaxM2, NBinP, MinP, MaxP);
  }
  // for TOF Param
  for(Int_t ip=0; ip<Event::nParticle; ++ip){
    HB2(10000+ip+1,
        Form("TofTime-%sTime %% TofSeg", name[ip].Data()),
        500, -25., 25., NumOfSegTOF, 0., (Double_t)NumOfSegTOF);
  }
  HB2(20001, "Tof TimeDiff U-D",
      1000, -25., 25., NumOfSegTOF, 0., (Double_t)NumOfSegTOF);
  for(Int_t i=0; i<NumOfSegTOF; ++i){
    HB1(20000+(i+1)*100+1,
	Form("Tof TimeDiff U-D %d", i+1),   1000, -25., 25.);
    HB1(30000+(i+1)*100+1,
        Form("TOF ADC-Pedestal %d-U", i+1), 4000, 0., 4000.);
    HB1(30000+(i+1)*100+11,
        Form("TOF ADC-Pedestal %d-D", i+1), 4000, 0., 4000.);
    HB2(40000+(i+1)*100,
        Form("tCalc-Time %% TOF De %d-U [All]", i+1), 100, 0., 4., 100, -3., 3.);
    HB2(40000+(i+1)*100+10,
        Form("tCalc-Time %% TOF De %d-D [All]", i+1), 100, 0., 4., 100, -3., 3.);
    for(Int_t ip=0; ip<Event::nParticle; ++ip){
      HB1(10000+(i+1)*100+ip+1,
          Form("Tof-%d TofTime-%sTime", i+1, name[ip].Data()),
          500, -25., 25.);
      HB1(33000+(i+1)*100+ip+1,
          Form("TOF ADC-Signal %d-U [%s]", i+1, name[ip].Data()),
          4000, 0., 4000.);
      HB1(33000+(i+1)*100+ip+11,
          Form("TOF ADC-Signal %d-D [%s]", i+1, name[ip].Data()),
          4000, 0., 4000.);
      HB2(40000+(i+1)*100+ip+1,
          Form("tCalc-Time %% TOF De %d-U [%s]", i+1, name[ip].Data()),
          100, 0., 4., 100, -3., 3.);
      HB2(40000+(i+1)*100+ip+11,
          Form("tCalc-Time %% TOF De %d-D [%s]", i+1, name[ip].Data()),
          100, 0., 4., 100, -3., 3.);
    }
  }

  HBTree("shodo", "tree of DstHypsHodoscope");
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("status",     &event.status,      "status/I");

  tree->Branch("nhBh1", &event.nhBh1, "nhBh1/I");
  tree->Branch("csBh1",  event.csBh1, "csBh1[nhBh1]/D");
  tree->Branch("Bh1Seg", event.Bh1Seg,"Bh1Seg[nhBh1]/D");
  tree->Branch("tBh1",   event.tBh1,  "tBh1[nhBh1]/D");
  tree->Branch("dtBh1",  event.dtBh1, "dtBh1[nhBh1]/D");
  tree->Branch("deBh1",  event.deBh1, "deBh1[nhBh1]/D");

  tree->Branch("nhBh2", &event.nhBh2, "nhBh2/I");
  tree->Branch("csBh2",  event.csBh2, "csBh2[nhBh2]/D");
  tree->Branch("Bh2Seg", event.Bh2Seg,"Bh2Seg[nhBh2]/D");
  tree->Branch("tBh2",   event.tBh2,  "tBh2[nhBh2]/D");
  tree->Branch("dtBh2",  event.dtBh2, "dtBh2[nhBh2]/D");
  tree->Branch("t0Bh2",  event.t0Bh2, "t0Bh2[nhBh2]/D");
  tree->Branch("deBh2",  event.deBh2, "deBh2[nhBh2]/D");

  tree->Branch("btof",   event.btof, "btof[nhBh1]/D");
  tree->Branch("Time0Seg", &event.Time0Seg, "Time0Seg/D");

  tree->Branch("nhTof",  &event.nhTof, "nhTof/I");
  tree->Branch("csTof",   event.csTof, "csTof[nhTof]/D");
  tree->Branch("TofSeg",  event.TofSeg,"TofSeg[nhTof]/D");
  tree->Branch("tTof",    event.tTof,  "tTof[nhTof]/D");
  tree->Branch("dtTof",   event.dtTof, "dtTof[nhTof]/D");
  tree->Branch("deTof",   event.deTof, "deTof[nhTof]/D");
  tree->Branch("tofua",   event.tofua, Form("tofua[%d]/D", NumOfSegTOF));
  tree->Branch("tofda",   event.tofda, Form("tofda[%d]/D", NumOfSegTOF));
  tree->Branch("utTofSeg",  event.utTofSeg,  Form("utTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("dtTofSeg",  event.dtTofSeg,  Form("dtTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("udeTofSeg", event.udeTofSeg, Form("udeTofSeg[%d]/D", NumOfSegTOF));
  tree->Branch("ddeTofSeg", event.ddeTofSeg, Form("ddeTofSeg[%d]/D", NumOfSegTOF));
  tree->Branch("nhAc1",     &event.nhAc1,     "nhAc1/I");
  tree->Branch("csAc1",      event.csAc1,     "csAc1[nhAc1]/I");
  tree->Branch("Ac1Seg",     event.Ac1Seg,    "Ac1Seg[nhAc1]/D");
  tree->Branch("tAc1",       event.tAc1,      "tAc1[nhAc1]/D");

  tree->Branch("ntSdcOut", &event.ntSdcOut,  "ntSdcOut/I");
  tree->Branch("chisqrSdcOut", event.chisqrSdcOut, "chisqrSdcOut[ntSdcOut]/D");
  tree->Branch("x0SdcOut", event.x0SdcOut, "x0SdcOut[ntSdcOut]/D");
  tree->Branch("y0SdcOut", event.y0SdcOut, "y0SdcOut[ntSdcOut]/D");
  tree->Branch("u0SdcOut", event.u0SdcOut, "u0SdcOut[ntSdcOut]/D");
  tree->Branch("v0SdcOut", event.v0SdcOut, "v0SdcOut[ntSdcOut]/D");
  tree->Branch("ntHyps", &event.ntHyps,  "ntHyps/I");
  tree->Branch("path",      event.path,      "path[ntHyps]/D");
  tree->Branch("pHyps",   event.pHyps,   "pHyps[ntHyps]/D");
  tree->Branch("qHyps",   event.qHyps,   "qHyps[ntHyps]/D");
  tree->Branch("chisqrHyps", event.chisqrHyps, "chisqrHyps[ntHyps]/D");
  tree->Branch("thetaHyps",   event.thetaHyps,  "thetaHyps[ntHyps]/D");
  tree->Branch("xtgtHyps",    event.xtgtHyps,   "xtgtHyps[ntHyps]/D");
  tree->Branch("ytgtHyps",    event.ytgtHyps,   "ytgtHyps[ntHyps]/D");
  tree->Branch("utgtHyps",    event.utgtHyps,   "utgtHyps[ntHyps]/D");
  tree->Branch("vtgtHyps",    event.vtgtHyps,   "vtgtHyps[ntHyps]/D");
  tree->Branch("xtofHyps",    event.xtofHyps,   "xtofHyps[ntHyps]/D");
  tree->Branch("ytofHyps",    event.ytofHyps,   "ytofHyps[ntHyps]/D");
  tree->Branch("utofHyps",    event.utofHyps,   "utofHyps[ntHyps]/D");
  tree->Branch("vtofHyps",    event.vtofHyps,   "vtofHyps[ntHyps]/D");
  tree->Branch("lxtofHyps",    event.lxtofHyps,   "lxtofHyps[ntHyps]/D");
  tree->Branch("lytofHyps",    event.lytofHyps,   "lytofHyps[ntHyps]/D");
  tree->Branch("lutofHyps",    event.lutofHyps,   "lutofHyps[ntHyps]/D");
  tree->Branch("lvtofHyps",    event.lvtofHyps,   "lvtofHyps[ntHyps]/D");
  tree->Branch("tofsegHyps",  event.tofsegHyps, "tofsegHyps[ntHyps]/D");
  tree->Branch("vpx", event.vpx, Form("vpx[%d]/D", NumOfLayersVP));
  tree->Branch("vpy", event.vpy, Form("vpy[%d]/D", NumOfLayersVP));
  tree->Branch("vpu", event.vpu, Form("vpu[%d]/D", NumOfLayersVP));
  tree->Branch("vpv", event.vpv, Form("vpv[%d]/D", NumOfLayersVP));

  tree->Branch("tTofCalc",  event.tTofCalc,  Form("tTofCalc[%d]/D", Event::nParticle));

  tree->Branch("m2Combi", &event.m2Combi, "m2Combi/I");
  tree->Branch("beta",     event.beta,    "beta[m2Combi]/D");
  tree->Branch("stof",     event.stof,    "stof[m2Combi]/D");
  tree->Branch("cstof",    event.cstof,   "cstof[m2Combi]/D");
  tree->Branch("m2",       event.m2,      "m2[m2Combi]/D");

  ////////// Bring Address From Dst
  TTreeCont[kHodoscope]->SetBranchStatus("*", 0);
  TTreeCont[kHodoscope]->SetBranchStatus("trigflag",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("trigpat",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhBh1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("csBh1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("Bh1Seg",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("tBh1",      1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtBh1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("deBh1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhBh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("csBh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("Bh2Seg",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("tBh2",      1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtBh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("t0Bh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("deBh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("Time0Seg",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("deTime0",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("Time0",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("CTime0",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("csTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("TofSeg",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("tTof",      1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("deTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("tofua",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("tofda",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("utTofSeg",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtTofSeg",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("udeTofSeg", 1);
  TTreeCont[kHodoscope]->SetBranchStatus("ddeTofSeg", 1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhAc1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("csAc1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("Ac1Seg",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("tAc1",      1);

  TTreeCont[kHodoscope]->SetBranchAddress("trigflag", src.trigflag);
  TTreeCont[kHodoscope]->SetBranchAddress("trigpat",  src.trigpat);
  TTreeCont[kHodoscope]->SetBranchAddress("nhBh1",  &src.nhBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("csBh1",  src.csBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("Bh1Seg", src.Bh1Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("tBh1",   src.tBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("dtBh1",  src.dtBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("deBh1",  src.deBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("nhBh2", &src.nhBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("csBh2",  src.csBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("Bh2Seg",src.Bh2Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("tBh2",  src.tBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("dtBh2", src.dtBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("t0Bh2", src.t0Bh2);
  TTreeCont[kHodoscope]->SetBranchAddress("deBh2", src.deBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("Time0Seg", &src.Time0Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("deTime0",  &src.deTime0);
  TTreeCont[kHodoscope]->SetBranchAddress("Time0",    &src.Time0);
  TTreeCont[kHodoscope]->SetBranchAddress("CTime0",   &src.CTime0);
  TTreeCont[kHodoscope]->SetBranchAddress("nhTof", &src.nhTof);
  TTreeCont[kHodoscope]->SetBranchAddress("csTof", src.csTof);
  TTreeCont[kHodoscope]->SetBranchAddress("TofSeg",src.TofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("tTof",  src.tTof);
  TTreeCont[kHodoscope]->SetBranchAddress("dtTof", src.dtTof);
  TTreeCont[kHodoscope]->SetBranchAddress("deTof", src.deTof);
  TTreeCont[kHodoscope]->SetBranchAddress("tofua", src.tofua);
  TTreeCont[kHodoscope]->SetBranchAddress("tofda", src.tofda);
  TTreeCont[kHodoscope]->SetBranchAddress("utTofSeg", src.utTofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("dtTofSeg", src.dtTofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("udeTofSeg", src.udeTofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("ddeTofSeg", src.ddeTofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("nhAc1", &src.nhAc1);
  TTreeCont[kHodoscope]->SetBranchAddress("csAc1", src.csAc1);
  TTreeCont[kHodoscope]->SetBranchAddress("Ac1Seg",src.Ac1Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("tAc1",  src.tAc1);

  TTreeCont[kHypsTracking]->SetBranchStatus("*",      0);
  TTreeCont[kHypsTracking]->SetBranchStatus("ntSdcOut", 1);
  TTreeCont[kHypsTracking]->SetBranchStatus("chisqrSdcOut", 1);
  TTreeCont[kHypsTracking]->SetBranchStatus("x0SdcOut", 1);
  TTreeCont[kHypsTracking]->SetBranchStatus("y0SdcOut", 1);
  TTreeCont[kHypsTracking]->SetBranchStatus("u0SdcOut", 1);
  TTreeCont[kHypsTracking]->SetBranchStatus("v0SdcOut", 1);
  TTreeCont[kHypsTracking]->SetBranchStatus("ntHyps",     1);
  TTreeCont[kHypsTracking]->SetBranchStatus("path",         1);
  TTreeCont[kHypsTracking]->SetBranchStatus("pHyps",      1);
  TTreeCont[kHypsTracking]->SetBranchStatus("qHyps",      1);
  TTreeCont[kHypsTracking]->SetBranchStatus("chisqrHyps", 1);
  TTreeCont[kHypsTracking]->SetBranchStatus("thetaHyps",  1);
  TTreeCont[kHypsTracking]->SetBranchStatus("xtgtHyps",   1);
  TTreeCont[kHypsTracking]->SetBranchStatus("ytgtHyps",   1);
  TTreeCont[kHypsTracking]->SetBranchStatus("utgtHyps",   1);
  TTreeCont[kHypsTracking]->SetBranchStatus("vtgtHyps",   1);
  TTreeCont[kHypsTracking]->SetBranchStatus("xtofHyps",   1);
  TTreeCont[kHypsTracking]->SetBranchStatus("ytofHyps",   1);
  TTreeCont[kHypsTracking]->SetBranchStatus("utofHyps",   1);
  TTreeCont[kHypsTracking]->SetBranchStatus("vtofHyps",   1);
  TTreeCont[kHypsTracking]->SetBranchStatus("lxtofHyps",   1);
  TTreeCont[kHypsTracking]->SetBranchStatus("lytofHyps",   1);
  TTreeCont[kHypsTracking]->SetBranchStatus("lutofHyps",   1);
  TTreeCont[kHypsTracking]->SetBranchStatus("lvtofHyps",   1);
  TTreeCont[kHypsTracking]->SetBranchStatus("tofsegHyps", 1);
  TTreeCont[kHypsTracking]->SetBranchStatus("vpx",          1);
  TTreeCont[kHypsTracking]->SetBranchStatus("vpy",          1);
  TTreeCont[kHypsTracking]->SetBranchStatus("vpu",          1);
  TTreeCont[kHypsTracking]->SetBranchStatus("vpv",          1);

  TTreeCont[kHypsTracking]->SetBranchAddress("ntSdcOut", &src.ntSdcOut);
  TTreeCont[kHypsTracking]->SetBranchAddress("chisqrSdcOut", src.chisqrSdcOut);
  TTreeCont[kHypsTracking]->SetBranchAddress("x0SdcOut", src.x0SdcOut);
  TTreeCont[kHypsTracking]->SetBranchAddress("y0SdcOut", src.y0SdcOut);
  TTreeCont[kHypsTracking]->SetBranchAddress("u0SdcOut", src.u0SdcOut);
  TTreeCont[kHypsTracking]->SetBranchAddress("v0SdcOut", src.v0SdcOut);
  TTreeCont[kHypsTracking]->SetBranchAddress("ntHyps", &src.ntHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("path",     src.path);
  TTreeCont[kHypsTracking]->SetBranchAddress("pHyps",  src.pHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("qHyps",  src.qHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("chisqrHyps", src.chisqrHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("thetaHyps", src.thetaHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("xtgtHyps", src.xtgtHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("ytgtHyps", src.ytgtHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("utgtHyps", src.utgtHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("vtgtHyps", src.vtgtHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("xtofHyps",   src.xtofHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("ytofHyps",   src.ytofHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("utofHyps",   src.utofHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("vtofHyps",   src.vtofHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("lxtofHyps",   src.lxtofHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("lytofHyps",   src.lytofHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("lutofHyps",   src.lutofHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("lvtofHyps",   src.lvtofHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("tofsegHyps", src.tofsegHyps);
  TTreeCont[kHypsTracking]->SetBranchAddress("vpx",        src.vpx);
  TTreeCont[kHypsTracking]->SetBranchAddress("vpy",        src.vpy);
  TTreeCont[kHypsTracking]->SetBranchAddress("vpu",        src.vpu);
  TTreeCont[kHypsTracking]->SetBranchAddress("vpv",        src.vpv);

  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<UserParamMan>("USER") &&
     InitializeParameter<HodoPHCMan>("HDPHC"));
}

//_____________________________________________________________________
bool
ConfMan::FinalizeProcess()
{
  return true;
}
