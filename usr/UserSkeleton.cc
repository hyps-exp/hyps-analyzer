// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include <UnpackerManager.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "HodoRawHit.hh"
#include "HypsLib.hh"
#include "RawData.hh"

namespace
{
using namespace root;
auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
}

//_____________________________________________________________________________
struct Event
{
  Int_t runnum;
  Int_t evnum;
  Int_t spill;
  void clear()
    {
      runnum = -1;
      evnum = -1;
      spill = -1;
    }
};

//_____________________________________________________________________________
namespace root
{
Event  event;
TH1   *h[MaxHist];
TTree *tree;
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
  event.runnum = gUnpacker.get_run_number();
  event.evnum  = gUnpacker.get_event_number();
  // rawData = new RawData;
  // rawData->DecodeHits();

  // for(Int_t i=0; i<100; ++i){
  //   HF1(i, (double)i);
  // }

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
  for(Int_t i=0; i<100; ++i){
    HB1(i, Form("hist %d", i), 100, 0., 100.);
  }

  HBTree("skeleton","tree of Skeleton");
  tree->Branch("runnum", &event.runnum, "runnum/I");
  tree->Branch("evnum",  &event.evnum,  "evnum/I");
  tree->Branch("spill",  &event.spill,  "spill/I");

  HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
     InitializeParameter<HodoParamMan>("HDPRM"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
