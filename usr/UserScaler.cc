// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <filesystem_util.hh>
#include <lexical_cast.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "ScalerAnalyzer.hh"
#include "Unpacker.hh"
#include "UnpackerManager.hh"

#define USE_COMMA   0
#define SPILL_RESET 0
#define MAKE_LOG    0

namespace
{
using namespace root;
const auto qnan = TMath::QuietNaN();
// using hddaq::unpacker::GUnpacker;
auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto& gScaler   = ScalerAnalyzer::GetInstance();
// auto& gRM       = RMAnalyzer::GetInstance();
}

//_____________________________________________________________________________
struct Event
{
  Int_t    evnum;
  Int_t    spill;
  Int_t    L1req;
  Int_t    L1acc;
  Int_t    L2req;
  Int_t    L2acc;
  Int_t    realtime;
  Int_t    livetime;
  Double_t daqeff;
  Double_t l2eff;
  Double_t real_live;
  Double_t duty;
  Int_t    rf;
  Int_t    tagger_or_1;
  Int_t    tagger_or_2;
  Int_t    tagger_coin;
  Int_t    t0;
  Int_t    t0_L;
  Int_t    t0_R;
  Int_t    sac_sum;
  Int_t    eVeto;
  Int_t    eVeto_L;
  Int_t    eVeto_R;
  Int_t    tof;
  Int_t    tof_or_1;
  Int_t    tof_or_2;
  Int_t    coin_2nd;
  Int_t    bgo;
  Int_t    bgo_or_1;
  Int_t    bgo_or_2;
  Int_t    cft_phi1;
  Int_t    cft_phi2;
  Int_t    cft_phi3;
  Int_t    cft_phi4;

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum       = 0;
  spill       = 0;
  L1req       = 0;
  L1acc       = 0;
  L2req       = 0;
  L2acc       = 0;
  realtime    = 0;
  livetime    = 0;
  daqeff      = 0.;
  l2eff       = 0.;
  real_live   = 0.;
  duty        = 0.;
  rf          = 0;
  tagger_or_1 = 0;
  tagger_or_2 = 0;
  tagger_coin = 0;
  t0          = 0;
  t0_L        = 0;
  t0_R        = 0;
  sac_sum     = 0;
  eVeto       = 0;
  eVeto_L     = 0;
  eVeto_R     = 0;
  tof         = 0;
  tof_or_1    = 0;
  tof_or_2    = 0;
  coin_2nd    = 0;
  bgo         = 0;
  bgo_or_1    = 0;
  bgo_or_2    = 0;
  cft_phi1    = 0;
  cft_phi2    = 0;
  cft_phi3    = 0;
  cft_phi4    = 0;
}

//_____________________________________________________________________________
struct Spill
{
  Int_t    evnum;
  Int_t    spill;
  Int_t    L1req;
  Int_t    L1acc;
  Int_t    L2req;
  Int_t    L2acc;
  Int_t    realtime;
  Int_t    livetime;
  Double_t daqeff;
  Double_t l2eff;
  Double_t real_live;
  Double_t duty;
  Int_t    rf;
  Int_t    tagger_or_1;
  Int_t    tagger_or_2;
  Int_t    tagger_coin;
  Int_t    t0;
  Int_t    t0_L;
  Int_t    t0_R;
  Int_t    sac_sum;
  Int_t    eVeto;
  Int_t    eVeto_L;
  Int_t    eVeto_R;
  Int_t    tof;
  Int_t    tof_or_1;
  Int_t    tof_or_2;
  Int_t    coin_2nd;
  Int_t    bgo;
  Int_t    bgo_or_1;
  Int_t    bgo_or_2;
  Int_t    cft_phi1;
  Int_t    cft_phi2;
  Int_t    cft_phi3;
  Int_t    cft_phi4;

  void clear();
};

//_____________________________________________________________________________
void
Spill::clear()
{
  evnum       = 0;
  spill       = 0;
  L1req       = 0;
  L1acc       = 0;
  L2req       = 0;
  L2acc       = 0;
  realtime    = 0;
  livetime    = 0;
  daqeff      = 0.;
  l2eff       = 0.;
  real_live   = 0.;
  duty        = 0.;
  rf          = 0;
  tagger_or_1 = 0;
  tagger_or_2 = 0;
  tagger_coin = 0;
  t0          = 0;
  t0_L        = 0;
  t0_R        = 0;
  sac_sum     = 0;
  eVeto       = 0;
  eVeto_L     = 0;
  eVeto_R     = 0;
  tof         = 0;
  tof_or_1    = 0;
  tof_or_2    = 0;
  coin_2nd    = 0;
  bgo         = 0;
  bgo_or_1    = 0;
  bgo_or_2    = 0;
  cft_phi1    = 0;
  cft_phi2    = 0;
  cft_phi3    = 0;
  cft_phi4    = 0;
}

//_____________________________________________________________________________
// struct Dst
// {
//   Int_t evnum;
//   Int_t spill;
//   void clear();
// };

//_____________________________________________________________________________
// void
// Dst::clear()
// {
//   evnum = 0;
//   spill = 0;
// }

//_____________________________________________________________________________
namespace root
{
// Event  event;
Spill  spill;
// Dst    dst;
TH1*   h[MaxHist];
// TTree* event;
TTree* tree;
Int_t  spill_cnt;
}

//_____________________________________________________________________________
Bool_t
ProcessingBegin()
{
  // event.clear();
  // dst.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingNormal()
{
  // RawData rawData;
  // gRM.Decode();

  gScaler.Decode();

  if( gScaler.SpillIncrement() ){
    tree->Fill();
    spill.clear();
    spill_cnt++;
  }

  // event.evnum       = gUnpacker.get_event_number();
  // event.spill       = spill;
  // event.L1req       += gScaler.Get("L1-Req");
  // event.L1acc       += gScaler.Get("L1-Acc");
  // event.L2req       += gScaler.Get("L2-Req");
  // event.L2acc       += gScaler.Get("L2-Req"); // FastClear is not used
  // event.realtime    += gScaler.Get("Real-Time");
  // event.livetime    += gScaler.Get("Live-Time");
  // event.daqeff      = static_cast<Double_t>event.L1acc / event.L1req;
  // event.L2eff       = static_cast<Double_t>event.L2acc / event.L2req;
  // event.real_live   = static_cast<Double_t>event.livetime / event.realtime;
  // event.duty        = gScaler.Duty();
  // event.rf          += gScaler.Get("RF");
  // event.tagger_or_1 += gScaler.Get("Tagger_OR-1");
  // event.tagger_or_2 += gScaler.Get("Tagger_OR-2");
  // event.tagger_coin += gScaler.Get("Tagger_Coin");
  // event.t0          += gScaler.Get("T0");
  // event.t0_L        += gScaler.Get("T0-L");
  // event.t0_R        += gScaler.Get("T0-R");
  // event.sac_sum     += gScaler.Get("SAC-Sum");
  // event.eVeto       += gScaler.Get("e-Veto");
  // event.eVeto_L     += gScaler.Get("e-Veto-L");
  // event.eVeto_R     += gScaler.Get("e-Veto-R");
  // event.tof         += gScaler.Get("TOF");
  // event.tof_or_1    += gScaler.Get("MT-TOF_OR-1");
  // event.tof_or_2    += gScaler.Get("MT-TOF_OR-2");
  // event.coin_2nd    += gScaler.Get("Coin-2nd-Stage");
  // event.bgo         += gScaler.Get("BGO");
  // event.bgo_or_1    += gScaler.Get("BGO_OR-1");
  // event.bgo_or_2    += gScaler.Get("BGO_OR-2");
  // event.cft_phi1    += gScaler.Get("CFT-Phi1");
  // event.cft_phi2    += gScaler.Get("CFT-Phi2");
  // event.cft_phi3    += gScaler.Get("CFT-Phi3");
  // event.cft_phi4    += gScaler.Get("CFT-Phi4");

  spill.evnum       = gUnpacker.get_event_number();
  spill.spill       = spill_cnt;
  spill.L1req       = gScaler.Get("L1-Req");
  spill.L1acc       = gScaler.Get("L1-Acc");
  spill.L2req       = gScaler.Get("L2-Req");
  spill.L2acc       = gScaler.Get("L2-Req"); // FastClear is not used
  spill.realtime    = gScaler.Get("Real-Time");
  spill.livetime    = gScaler.Get("Live-Time");
  spill.daqeff      = gScaler.Fraction("L1-Acc", "L1-Req");
  spill.l2eff       = gScaler.Fraction("L2-Acc", "L2-Req");
  spill.real_live   = gScaler.Fraction("Live-Time", "Real-Time");
  spill.duty        = gScaler.Duty();
  spill.rf          = gScaler.Get("RF");
  spill.tagger_or_1 = gScaler.Get("Tagger_OR-1");
  spill.tagger_or_2 = gScaler.Get("Tagger_OR-2");
  spill.tagger_coin = gScaler.Get("Tagger_Coin");
  spill.t0          = gScaler.Get("T0");
  spill.t0_L        = gScaler.Get("T0-L");
  spill.t0_R        = gScaler.Get("T0-R");
  spill.sac_sum     = gScaler.Get("SAC-Sum");
  spill.eVeto       = gScaler.Get("e-Veto");
  spill.eVeto_L     = gScaler.Get("e-Veto-L");
  spill.eVeto_R     = gScaler.Get("e-Veto-R");
  spill.tof         = gScaler.Get("TOF");
  spill.tof_or_1    = gScaler.Get("MT-TOF_OR-1");
  spill.tof_or_2    = gScaler.Get("MT-TOF_OR-2");
  spill.coin_2nd    = gScaler.Get("Coin-2nd-Stage");
  spill.bgo         = gScaler.Get("BGO");
  spill.bgo_or_1    = gScaler.Get("BGO_OR-1");
  spill.bgo_or_2    = gScaler.Get("BGO_OR-2");
  spill.cft_phi1    = gScaler.Get("CFT-Phi1");
  spill.cft_phi2    = gScaler.Get("CFT-Phi2");
  spill.cft_phi3    = gScaler.Get("CFT-Phi3");
  spill.cft_phi4    = gScaler.Get("CFT-Phi4");

// #if !MAKE_LOG
//   if(event.evnum%400==0)
//     gScaler.Print();
// #endif

// #if SPILL_RESET
//   if(gScaler.SpillIncrement())
//     gScaler.Clear();
// #endif

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingEnd()
{
  // event->Fill();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{

  spill_cnt = 0;

  spill.clear();

  { // hul01_scr-2
    Int_t c = ScalerAnalyzer::kLeft;
    Int_t r = 0;
    gScaler.Set( c, r++, ScalerInfo("1M-CLK",           1, 80) );
    gScaler.Set( c, r++, ScalerInfo("RF",               1, 81) );
    gScaler.Set( c, r++, ScalerInfo("Tagger_Coin",      1, 82) );
    gScaler.Set( c, r++, ScalerInfo("T0-L",             1, 83) );
    gScaler.Set( c, r++, ScalerInfo("T0-R",             1, 84) );
    gScaler.Set( c, r++, ScalerInfo("SAC-Sum",          1, 85) );
    gScaler.Set( c, r++, ScalerInfo("e-Veto-L",         1, 86) );
    gScaler.Set( c, r++, ScalerInfo("e-Veto-R",         1, 87) );
    gScaler.Set( c, r++, ScalerInfo("MT-TOF_OR-1",      1, 88) );
    gScaler.Set( c, r++, ScalerInfo("MT-TOF_OR-2",      1, 89) );
    gScaler.Set( c, r++, ScalerInfo("MT-T0",            1, 90) );
    gScaler.Set( c, r++, ScalerInfo("MT-e-Veto",        1, 91) );
    gScaler.Set( c, r++, ScalerInfo("BGO_OR-1",         1, 92) );
    gScaler.Set( c, r++, ScalerInfo("BGO_OR-2",         1, 93) );
    gScaler.Set( c, r++, ScalerInfo("Tagger_OR-1",      1, 94) );
    gScaler.Set( c, r++, ScalerInfo("Tagger_OR-2",      1, 95) );

    // hul03_scr-1
    gScaler.Set( c, r++, ScalerInfo("Real-Time",        2,  0) );
    gScaler.Set( c, r++, ScalerInfo("L1-Req",           2,  1) );
    gScaler.Set( c, r++, ScalerInfo("L1-Acc",           2,  2) );
    gScaler.Set( c, r++, ScalerInfo("L2-Req",           2,  3) );
    gScaler.Set( c, r++, ScalerInfo("L2-Acc",           2,  3) ); // FastClear is not used
    gScaler.Set( c, r++, ScalerInfo("T0",               2,  4) );
    gScaler.Set( c, r++, ScalerInfo("TOF",              2,  5) );
    gScaler.Set( c, r++, ScalerInfo("SAC",              2,  6) );
    gScaler.Set( c, r++, ScalerInfo("e-Veto",           2,  7) );
    gScaler.Set( c, r++, ScalerInfo("Coin-2nd-Stage",   2,  8) );
    gScaler.Set( c, r++, ScalerInfo("BGO",              2,  9) );
    gScaler.Set( c, r++, ScalerInfo("CFT-Phi1",         2, 10) );
    gScaler.Set( c, r++, ScalerInfo("CFT-Phi2",         2, 11) );
    gScaler.Set( c, r++, ScalerInfo("CFT-Phi3",         2, 12) );
    gScaler.Set( c, r++, ScalerInfo("CFT-Phi4",         2, 13) );
    gScaler.Set( c, r++, ScalerInfo("Live-Time",        2, 14) );
  }

  //Tree
  // HBTree("event", "tree of Scaler");
  // tree->Branch("evnum",           &event.evnum,           "evnum/I");
  // tree->Branch("spill",           &event.spill,           "spill/I");
  // tree->Branch("L1req",           &event.L1req,           "L1req/I");
  // tree->Branch("L1acc",           &event.L1acc,           "L1acc/I");
  // tree->Branch("L2req",           &event.L2req,           "L2req/I");
  // tree->Branch("L2acc",           &event.L2acc,           "L2acc/I");
  // tree->Branch("realtime",        &event.realtime,        "realtime/I");
  // tree->Branch("livetime",        &event.livetime,        "livetime/I");
  // tree->Branch("daqeff",          &event.daqeff,          "daqeff/D");
  // tree->Branch("l2eff",           &event.l2eff,           "l2eff/D");
  // tree->Branch("real_live",       &event.real_live,       "real_live/D");
  // tree->Branch("duty",            &event.duty,            "duty/D");
  // tree->Branch("rf",              &event.rf,              "rf/I");
  // tree->Branch("tagger_or_1",     &event.tagger_or_1,     "tagger_or_1/I");
  // tree->Branch("tagger_or_2",     &event.tagger_or_2,     "tagger_or_2/I");
  // tree->Branch("tagger_coin",     &event.tagger_coin,     "tagger_coin/I");
  // tree->Branch("t0",              &event.t0,              "t0/I");
  // tree->Branch("t0_L",            &event.t0_L,            "t0_L/I");
  // tree->Branch("t0_R",            &event.t0_R,            "t0_R/I");
  // tree->Branch("sac_sum",         &event.sac_sum,         "sac_sum/I");
  // tree->Branch("eVeto",           &event.eVeto,           "eVeto/I");
  // tree->Branch("eVeto_L",         &event.eVeto_L,         "eVeto_L/I");
  // tree->Branch("eVeto_R",         &event.eVeto_R,         "eVeto_R/I");
  // tree->Branch("tof",             &event.tof,             "tof/I");
  // tree->Branch("tof_or_1",        &event.tof_or_1,        "tof_or_1/I");
  // tree->Branch("tof_or_2",        &event.tof_or_2,        "tof_or_2/I");
  // tree->Branch("coin_2nd",        &event.coin_2nd,        "coin_2nd/I");
  // tree->Branch("bgo",             &event.bgo,             "bgo/I");
  // tree->Branch("bgo_or_1",        &event.bgo_or_1,        "bgo_or_1/I");
  // tree->Branch("bgo_or_2",        &event.bgo_or_2,        "bgo_or_2/I");
  // tree->Branch("cft_phi1",        &event.cft_phi1,        "cft_phi1/I");
  // tree->Branch("cft_phi2",        &event.cft_phi2,        "cft_phi2/I");
  // tree->Branch("cft_phi3",        &event.cft_phi3,        "cft_phi3/I");
  // tree->Branch("cft_phi4",        &event.cft_phi4,        "cft_phi4/I");


  HBTree("spill", "Scaler data in each spill");
  tree->Branch("evnum",           &spill.evnum,           "evnum/I");
  tree->Branch("spill",           &spill.spill,           "spill/I");
  tree->Branch("L1req",           &spill.L1req,           "L1req/I");
  tree->Branch("L1acc",           &spill.L1acc,           "L1acc/I");
  tree->Branch("L2req",           &spill.L2req,           "L2req/I");
  tree->Branch("L2acc",           &spill.L2acc,           "L2acc/I");
  tree->Branch("realtime",        &spill.realtime,        "realtime/I");
  tree->Branch("livetime",        &spill.livetime,        "livetime/I");
  tree->Branch("daqeff",          &spill.daqeff,          "daqeff/D");
  tree->Branch("l2eff",           &spill.l2eff,           "l2eff/D");
  tree->Branch("real_live",       &spill.real_live,       "real_live/D");
  tree->Branch("duty",            &spill.duty,            "duty/D");
  tree->Branch("rf",              &spill.rf,              "rf/I");
  tree->Branch("tagger_or_1",     &spill.tagger_or_1,     "tagger_or_1/I");
  tree->Branch("tagger_or_2",     &spill.tagger_or_2,     "tagger_or_2/I");
  tree->Branch("tagger_coin",     &spill.tagger_coin,     "tagger_coin/I");
  tree->Branch("t0",              &spill.t0,              "t0/I");
  tree->Branch("t0_L",            &spill.t0_L,            "t0_L/I");
  tree->Branch("t0_R",            &spill.t0_R,            "t0_R/I");
  tree->Branch("sac_sum",         &spill.sac_sum,         "sac_sum/I");
  tree->Branch("eVeto",           &spill.eVeto,           "eVeto/I");
  tree->Branch("eVeto_L",         &spill.eVeto_L,         "eVeto_L/I");
  tree->Branch("eVeto_R",         &spill.eVeto_R,         "eVeto_R/I");
  tree->Branch("tof",             &spill.tof,             "tof/I");
  tree->Branch("tof_or_1",        &spill.tof_or_1,        "tof_or_1/I");
  tree->Branch("tof_or_2",        &spill.tof_or_2,        "tof_or_2/I");
  tree->Branch("coin_2nd",        &spill.coin_2nd,        "coin_2nd/I");
  tree->Branch("bgo",             &spill.bgo,             "bgo/I");
  tree->Branch("bgo_or_1",        &spill.bgo_or_1,        "bgo_or_1/I");
  tree->Branch("bgo_or_2",        &spill.bgo_or_2,        "bgo_or_2/I");
  tree->Branch("cft_phi1",        &spill.cft_phi1,        "cft_phi1/I");
  tree->Branch("cft_phi2",        &spill.cft_phi2,        "cft_phi2/I");
  tree->Branch("cft_phi3",        &spill.cft_phi3,        "cft_phi3/I");
  tree->Branch("cft_phi4",        &spill.cft_phi4,        "cft_phi4/I");


#if USE_COMMA
  gScaler.SetFlag(ScalerAnalyzer::kSeparateComma);
#endif

  // gScaler.SetFlag(ScalerAnalyzer::kSpillOn);

  gScaler.SetFlag(ScalerAnalyzer::kSpillBySpill);

  gScaler.PrintFlags();

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
#if 0
  if(event.evnum==0) return true;

  gScaler.Print();

#if MAKE_LOG
  const Int_t run_number = gUnpacker.get_root()->get_run_number();
  const TString& bin_dir(hddaq::dirname(hddaq::selfpath()));
  const TString& data_dir(hddaq::dirname(gUnpacker.get_istream()));

  std::stringstream run_number_ss; run_number_ss << run_number;
  const TString& recorder_log(data_dir+"/recorder.log");
  std::ifstream ifs(recorder_log);
  if(!ifs.is_open()){
    std::cerr << "#E " << FUNC_NAME << " "
	      << "cannot open recorder.log : "
	      << recorder_log << std::endl;
    return false;
  }

  const TString& scaler_dir(bin_dir+"/../auto_scaler");
  const TString& scaler_txt = Form("%s/scaler_%05d.txt",
				    scaler_dir.Data(), run_number);

  std::ofstream ofs(scaler_txt);
  if(!ofs.is_open()){
    std::cerr << "#E " << FUNC_NAME << " "
	      << "cannot open scaler.txt : "
	      << scaler_txt << std::endl;
    return false;
  }

  Int_t recorder_event_number = 0;
  Bool_t found_run_number = false;
  std::string line;
  while(ifs.good() && std::getline(ifs,line)){
    if(line.empty()) continue;
    std::istringstream input_line(line);
    std::istream_iterator<std::string> line_begin(input_line);
    std::istream_iterator<std::string> line_end;
    std::vector<std::string> log_column(line_begin, line_end);
    if(log_column.at(0) != "RUN") continue;
    if(log_column.at(1) != run_number_ss.str()) continue;
    recorder_event_number = hddaq::a2i(log_column.at(15));
    ofs << line << std::endl;
    found_run_number = true;
  }

  if(!found_run_number){
    std::cerr << "#E " << FUNC_NAME << " "
	      << "not found run# " << run_number
	      << " in " << recorder_log << std::endl;
    return false;
  }

  ofs << std::endl;
  ofs << std::left  << std::setw(15) << "" << "\t"
      << std::right << std::setw(15) << "Integral" << std::endl;
  ofs << std::left  << std::setw(15) << "Event"    << "\t"
      << std::right << std::setw(15) << event.evnum << std::endl;

  if(recorder_event_number != event.evnum){
    std::cerr << "#W " << FUNC_NAME << " "
	      << "event number mismatch" << std::endl
	      << "   recorder : " << recorder_event_number << std::endl
	      << "   decode   : " << event.evnum << std::endl;
  }

  {
    std::vector<Int_t> order = {
      ScalerAnalyzer::kRight,
      ScalerAnalyzer::kLeft,
      ScalerAnalyzer::kCenter
    };
    for(auto&& c : order){
      for(Int_t i=0; i<ScalerAnalyzer::MaxRow; i++){
	TString name = gScaler.GetScalerName(c, i);
	if(name=="n/a") continue;
	ofs << std::left  << std::setw(15) << name << "\t"
	    << std::right << std::setw(15) << gScaler.Get(name) << std::endl;
      }
      ofs << std::endl;
    }
  }

  Double_t reallive = gScaler.Fraction("Live-Time", "Real-Time");
  Double_t daqeff   = gScaler.Fraction("L1-Acc", "L1-Req");
  Double_t l2eff    = gScaler.Fraction("L2-Acc", "L1-Acc");
  Double_t beamtm   = gScaler.Fraction("Beam", "TM");
  Double_t kbeamtm  = gScaler.Fraction("K-Beam", "TM");
  Double_t pibeamtm = gScaler.Fraction("Pi-Beam", "TM");
  Double_t l1reqbeam = gScaler.Fraction("L1-Req", "Beam");
  Double_t beamrate   = gScaler.Fraction("Beam", "Spill");
  Double_t kbeamrate  = gScaler.Fraction("K-Beam", "Spill");
  Double_t pibeamrate = gScaler.Fraction("Pi-Beam", "Spill");
  Double_t l1rate   = gScaler.Fraction("L1-Req", "Spill");
  Double_t l2rate   = gScaler.Fraction("L2-Acc", "Spill");

  ofs << std::fixed << std::setprecision(6)
      << std::left  << std::setw(18) << "Live/Real"     << "\t"
      << std::right << std::setw(12)<<  reallive        << std::endl
      << std::left  << std::setw(18) << "DAQ-Eff"       << "\t"
      << std::right << std::setw(12) << daqeff          << std::endl
      << std::left  << std::setw(18) << "L2-Eff"        << "\t"
      << std::right << std::setw(12) << l2eff           << std::endl
      << std::left  << std::setw(18) << "Duty-Factor"   << "\t"
      << std::right << std::setw(12) << gScaler.Duty()  << std::endl
      << std::left  << std::setw(18) << "Beam/TM"       << "\t"
      << std::right << std::setw(12) << beamtm          << std::endl
      << std::left  << std::setw(18) << "K-Beam/TM"     << "\t"
      << std::right << std::setw(12) << kbeamtm         << std::endl
      << std::left  << std::setw(18) << "Pi-Beam/TM"    << "\t"
      << std::right << std::setw(12) << pibeamtm        << std::endl
      << std::left  << std::setw(18) << "L1-Req/Beam"   << "\t"
      << std::right << std::setw(12) << l1reqbeam       << std::endl
      << std::setprecision(0)
      << std::left  << std::setw(18) << "Beam/Spill"    << "\t"
      << std::right << std::setw(12) << beamrate        << std::endl
      << std::left  << std::setw(18) << "K-Beam/Spill"  << "\t"
      << std::right << std::setw(12) << kbeamrate       << std::endl
      << std::left  << std::setw(18) << "Pi-Beam/Spill" << "\t"
      << std::right << std::setw(12) << pibeamrate      << std::endl
      << std::left  << std::setw(18) << "L1-Req/Spill"  << "\t"
      << std::right << std::setw(12) << l1rate          << std::endl
      << std::left  << std::setw(18) << "L2-Acc/Spill"  << "\t"
      << std::right << std::setw(12) << l2rate          << std::endl;

#endif

#endif

  return true;
}
