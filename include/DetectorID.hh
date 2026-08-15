// -*- C++ -*-

#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH

#include <array>
#include <iostream>
#include <map>
#include <vector>

#include <TString.h>


// ======= Counters ============================================================
constexpr Int_t NumOfSegRF     =  1;
constexpr Int_t NumOfSegTagSF  = 55;
constexpr Int_t NumOfSegTagPL  = 10;
constexpr Int_t NumOfSegUpVeto =  1;
constexpr Int_t NumOfSegT0     =  1;
constexpr Int_t NumOfSegSAC    =  5;
constexpr Int_t NumOfSegE_Veto =  1;
constexpr Int_t NumOfSegTOF    = 48;

// ======= Trackers ============================================================
const std::map<TString, std::vector<TString>> DCNameList
  {
    {"SdcIn",  { "SDC0", "SDC1" }},
    {"SdcOut", { "SDC2", "SDC3" }},
  };

constexpr Int_t NumOfLayersSDC0 = 4;
constexpr Int_t NumOfLayersSDC1 = 6;
constexpr Int_t NumOfLayersSDC2 = 5;
constexpr Int_t NumOfLayersSDC3 = 5;

constexpr std::array<Int_t, NumOfLayersSDC0> NumOfWiresSDC0
  {
    112, // X
    112, // XP
    112, // U
    112  // UP
  };
constexpr std::array<Int_t, NumOfLayersSDC1> NumOfWiresSDC1
  {
    48, // XPP
    48, // V
    48, // UP
    48, // U
    48, // XP
    48  // X
  };
constexpr std::array<Int_t, NumOfLayersSDC2> NumOfWiresSDC2
  {
    79,  // V
    78,  // UP
    78,  // U
    104, // XP
    104  // X
  };
constexpr std::array<Int_t, NumOfLayersSDC3> NumOfWiresSDC3
  {
    79,  // V
    78,  // UP
    78,  // U
    104, // XP
    104  // X
  };
constexpr Int_t MaxWireSDC0 = 112;
constexpr Int_t MaxWireSDC1 =  48;
constexpr Int_t MaxWireSDC2 = 104;
constexpr Int_t MaxWireSDC3 = 104;

constexpr Double_t CellSizeSDC0 =  6.0;
constexpr Double_t CellSizeSDC1 =  6.0;
constexpr Double_t CellSizeSDC2 = 10.0;
constexpr Double_t CellSizeSDC3 = 10.0;

// ------- For compatibility with DCGeom ---------------------------------------
constexpr Int_t PlOffsSdcIn  =  0;
constexpr Int_t PlMinSdcIn   =  1;
constexpr Int_t PlMaxSdcIn   = 10;
constexpr Int_t PlOffsSdcOut = 30;
constexpr Int_t PlMinSdcOut  = 31;
constexpr Int_t PlMaxSdcOut  = 40;
constexpr Int_t PlOffsTOF    = 52;
constexpr Int_t PlMinTOF     = 53;
constexpr Int_t PlMaxTOF     = 58;
constexpr Int_t PlOffsVP     = 15;
constexpr Int_t PlMinVP      = 16;
constexpr Int_t PlMaxVP      = 16;
constexpr Int_t NumOfLayersSdcIn  = PlMaxSdcIn  - PlMinSdcIn  + 1;
constexpr Int_t NumOfLayersSdcOut = PlMaxSdcOut - PlMinSdcOut + 1;
constexpr Int_t NumOfLayersSdcInOut = NumOfLayersSdcIn + NumOfLayersSdcOut;
constexpr Int_t NumOfLayersTOF    = PlMaxTOF    - PlMinTOF    + 1;
constexpr Int_t NumOfLayersVP     = PlMaxVP     - PlMinVP     + 1;

// ======= CATCH ===============================================================
// ------- CFT -----------------------------------------------------------------
constexpr Int_t NumOfPlaneCFT = 8;
constexpr std::array<Int_t, NumOfPlaneCFT> NumOfSegCFT
  {
    426, // U1
    584, // PHI1
    472, // V2
    692, // PHI2
    510, // U3
    800, // PHI3
    538, // V4
    910  // PHI4
  };

constexpr Int_t NumOfSegCFT_UV1  = 426;
constexpr Int_t NumOfSegCFT_PHI1 = 584;
constexpr Int_t NumOfSegCFT_UV2  = 472;
constexpr Int_t NumOfSegCFT_PHI2 = 692;
constexpr Int_t NumOfSegCFT_UV3  = 510;
constexpr Int_t NumOfSegCFT_PHI3 = 800;
constexpr Int_t NumOfSegCFT_UV4  = 538;
constexpr Int_t NumOfSegCFT_PHI4 = 910;

enum class kCFTPlane : Int_t {
  U1,
  PHI1,
  V2,
  PHI2,
  U3,
  PHI3,
  V4,
  PHI4,
};

// ------- BGO -----------------------------------------------------------------
constexpr Int_t    NumOfSegBGO = 24;
constexpr Int_t    NumOfBGOUnit = 8;
constexpr Double_t BGO_X =  30.;
constexpr Double_t BGO_Y =  25.;
constexpr Double_t BGO_Z = 400.;
// Pair Unit
constexpr Int_t    NumOfBGOInOneUnit   =   2;
constexpr Double_t RadiusOfBGOSurface  = 100.;
// Single Unit
constexpr Int_t    NumOfBGOInOneUnit2  =   1;
constexpr Double_t RadiusOfBGOSurface2 = 120.;


// ------- PiID ----------------------------------------------------------------
constexpr Int_t NumOfSegPiID = 32;
constexpr Int_t NumOfPiIDUnit = 8;
// Pair Unit
constexpr Int_t    NumOfPiIDInOneUnit = 3;
constexpr Double_t PiID_X = 30.;
constexpr Double_t PiID_Y = 10.;
constexpr Double_t PiID_Z = 400.;
constexpr Double_t RadiusOfPiIDSurface = 164.;
// Single Unit
constexpr Int_t    NumOfPiIDInOneUnit2 = 1;
constexpr Double_t PiID2_X = 40.;
constexpr Double_t PiID2_Y = 10.;
constexpr Double_t PiID2_Z = 400.;
constexpr Double_t RadiusOfPiID2Surface = 180.;


// ======= Misc =================================================================
// ------- FEE -----------------------------------------------------------------
constexpr Int_t NumOfPlaneVmeEasiroc =  96;
constexpr Int_t NumOfSegVmeEasiroc   =  64;
constexpr Int_t NumOfPlaneHulRm      =   2;
constexpr Int_t NumOfPlaneVmeRm      =   2;
constexpr Int_t NumOfPlaneScaler     =   2;
constexpr Int_t NumOfSegScaler       =  96;
// constexpr Int_t DetIdVmeRm           =  81; // DIGIT before 2024.10.27
constexpr Int_t DetIdVmeRm           = 101; // DIGIT after 2024.10.27

// ------- Trigger Flag --------------------------------------------------------
enum class kTriggerFlag : Int_t {
  TrigAPS,
  TrigBPS,
  TrigCPS,
  TrigDPS,
  TrigEPS,
  TrigFPS,
  TrigPSORA,
  TrigPSORB,
  ClockPS,
  Reserve2PS,
  Level1OR,
  Clock10M,
  Clock1M,
  Clock100k,
  Clock10k,
  Clock1k,
  SpillOnStart,
  NTriggerFlag,
};
constexpr Int_t NumOfSegTrig = static_cast<Int_t>(kTriggerFlag::NTriggerFlag);

constexpr std::array<std::string_view, NumOfSegTrig> kTriggerFlagNames
  {
    "TrigA-PS",
    "TrigB-PS",
    "TrigC-PS",
    "TrigD-PS",
    "TrigE-PS",
    "TrigF-PS",
    "Trig-PSOR-A",
    "Trig-PSOR-B",
    "Clock-PS",
    "Reserve2-PS",
    "Level1-OR",
    "Clock-10M",
    "Clock-1M",
    "Clock-100k",
    "Clock-10k",
    "Clock-1k",
    "SpillOnStart"
  };

// ------- Unused in HYPS, but necessary to compile ----------------------------
constexpr Int_t NumOfSegBH1 = 11;
constexpr Int_t NumOfSegBH2 =  8;
constexpr Int_t NumOfSegSCH = 64;
constexpr Int_t NumOfSegWC  = 12;
constexpr Int_t NumOfSegAC1 = 30;

constexpr Int_t MaxWireBC3  = 64;
constexpr Int_t MaxWireBC4  = 64;

constexpr Int_t PlOffsBc    = 100;
constexpr Int_t PlMinBcIn   =   1;
constexpr Int_t PlMaxBcIn   =  12;
constexpr Int_t PlMinBcOut  = 113;
constexpr Int_t PlMaxBcOut  = 124;
constexpr Int_t NumOfLayersBcIn  = PlMaxBcIn  - PlMinBcIn  + 1;
constexpr Int_t NumOfLayersBcOut = PlMaxBcOut - PlMinBcOut + 1;

constexpr Int_t NumOfSegSFT_Mtx = 48;

#endif
