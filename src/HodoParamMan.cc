// -*- C++ -*-

#include "HodoParamMan.hh"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <TMath.h>

#include <std_ostream.hh>

#include "DeleteUtility.hh"
#include "FuncName.hh"

namespace
{
const Int_t SegMask  = 0x03FF;
const Int_t CidMask  = 0x00FF;
const Int_t PlidMask = 0x00FF;
const Int_t UdMask   = 0x0003;
const Int_t SegShift  =  0;
const Int_t CidShift  = 11;
const Int_t PlidShift = 19;
const Int_t UdShift   = 27;
const auto qnan = TMath::QuietNaN();
}

//_____________________________________________________________________________
HodoParamMan::HodoParamMan()
  : m_is_ready(false),
    m_file_name()
{
}

//_____________________________________________________________________________
HodoParamMan::~HodoParamMan()
{
  ClearAHCont();
  ClearALCont();
  ClearTCont();
}

//_____________________________________________________________________________
void
HodoParamMan::ClearAHCont()
{
  del::ClearMap(m_AHPContainer);
}

//_____________________________________________________________________________
void
HodoParamMan::ClearALCont()
{
  del::ClearMap(m_ALPContainer);
}

//_____________________________________________________________________________
void
HodoParamMan::ClearTCont()
{
  del::ClearMap(m_TPContainer);
}

//_____________________________________________________________________________
inline Int_t
KEY(Int_t cid, Int_t pl, Int_t seg, Int_t ud)
{
  return (((cid&CidMask) << CidShift )|
          ((pl&PlidMask) << PlidShift)|
          ((seg&SegMask) << SegShift )|
          ((ud&UdMask)   << UdShift  ));
}

//_____________________________________________________________________________
Bool_t
HodoParamMan::Initialize()
{
  if(m_is_ready){
    hddaq::cerr << FUNC_NAME
		<< " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    hddaq::cerr << FUNC_NAME << " file open fail : "
		<< m_file_name << std::endl;
    return false;
  }

  ClearAHCont();
  ClearALCont();
  ClearTCont();

  Int_t invalid=0;
  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    ++invalid;
    if(line.IsNull() || line[0]=='#') continue;
    std::istringstream input_line(line.Data());
    Int_t cid=-1, plid=-1, seg=-1, at=-1, ud=-1;
    Double_t p0=qnan, p1=qnan;
    Double_t p2=qnan, p3=qnan, p4=qnan, p5=qnan;
    if(input_line >> cid >> plid >> seg >> at >> ud >> p0 >> p1){
      Int_t key = KEY(cid, plid, seg, ud);
      if(at == kAdcHigh){
	auto pre_param = m_AHPContainer[key];
	auto param = new HodoAParam(p0, p1);
	m_AHPContainer[key] = param;
	if(pre_param){
	  hddaq::cerr << FUNC_NAME << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " key = " << key << std::endl;
	  delete pre_param;
	}
      }
      else if(at == kTdc){
	auto pre_param = m_TPContainer[key];
	auto param = new HodoTParam(p0, p1);
	m_TPContainer[key] = param;
	if(pre_param){
	  hddaq::cerr << FUNC_NAME << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " key = " << key << std::endl;
	  delete pre_param;
	}
      }
      else if(at == kAdcLow){
	auto pre_param = m_ALPContainer[key];
	auto param = new HodoAParam(p0, p1);
	m_ALPContainer[key] = param;
	if(pre_param){
	  hddaq::cerr << FUNC_NAME << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " key = " << key << std::endl;
	  delete pre_param;
	}
        // }else if(at == 3){// for fiber position correction
        //   if(input_line  >> p2 >> p3>> p4 >> p5){
        //     HodoFParam *pre_param = m_FPContainer[key];
        //     HodoFParam *param = new HodoFParam(p0, p1, p2, p3, p4, p5);
        //     m_FPContainer[key] = param;
        //     if(pre_param){
        //       hddaq::cerr << FUNC_NAME << ": duplicated key "
        //   		<< " following record is deleted." << std::endl
        //   		<< " key = " << key << std::endl;
        //       delete pre_param;
        //     }
        //   }
      }
      else if(at == kTdc2){
	auto pre_param = m_T2PContainer[key];
	auto param = new HodoTParam(p0, p1);
	m_T2PContainer[key] = param;
	if(pre_param){
	  hddaq::cerr << FUNC_NAME << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " key = " << key << std::endl;
	  delete pre_param;
	}
      }
 else{
	hddaq::cerr << FUNC_NAME << ": Invalid Input" << std::endl
		    << " ===> (" << invalid << "a)" << line << " " << std::endl;
      } /* if(at) */
    }
    else {
      hddaq::cerr << FUNC_NAME << ": Invalid Input" << std::endl
                  << " ===> (" << invalid << "b)" << line << " " << std::endl;
    } /* if(input_line >>) */
  } /* while(std::getline) */

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoParamMan::Initialize(const TString& file_name)
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
HodoParamMan::GetTime(Int_t cid, Int_t plid, Int_t seg,
                      Int_t ud, Int_t tdc, Double_t& time) const
{
  HodoTParam* map = GetTmap(cid, plid, seg, ud);
  if(!map) return false;
  time = map->Time(tdc);
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoParamMan::GetTime2(Int_t cid, Int_t plid, Int_t seg,
		       Int_t ud, Int_t tdc, Double_t& time) const
{
  HodoTParam* map = GetT2map(cid, plid, seg, ud);
  if(!map) return false;
  time = map->Time(tdc);
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoParamMan::GetDeHighGain(Int_t cid, Int_t plid, Int_t seg,
                            Int_t ud, Int_t adc, Double_t& de) const
{
  const auto& map = GetAHmap(cid, plid, seg, ud);
  if(!map) return false;
  de = map->DeltaE(adc);
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoParamMan::GetDeLowGain(Int_t cid, Int_t plid, Int_t seg,
                           Int_t ud, Int_t adc, Double_t& de) const
{
  const auto& map = GetALmap(cid, plid, seg, ud);
  if(!map) return false;
  de = map->DeltaE(adc);
  return true;
}

//____________________________________________________________________________
Double_t
HodoParamMan::GetP0(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const
{
  const auto& map = GetAHmap(cid, plid, seg, ud);
  if(!map) return -1;
  Double_t p0 = map->Pedestal();
  return p0;
}

//_____________________________________________________________________________
Double_t
HodoParamMan::GetP1(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const
{
  const auto& map = GetAHmap(cid, plid, seg, ud);
  if(!map) return -1;
  Double_t p1 = map->Gain();
  return p1;
}

//_____________________________________________________________________________
// Double_t
// HodoParamMan::GetPar(Int_t cid, Int_t plid, Int_t seg, Int_t ud, Int_t i) const
// {
//   HodoFParam *map = GetFmap(cid, plid, seg, ud);
//   if(!map) return -1;
//   Double_t par=0;
//   if(i==0) par=map->par0();
//   else if(i==1) par=map->par1();
//   else if(i==2) par=map->par2();
//   else if(i==3) par=map->par3();
//   else if(i==4) par=map->par4();
//   else if(i==5) par=map->par5();
//   return par;
// }

//_____________________________________________________________________________
Double_t
HodoParamMan::GetOffset(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const
{
  HodoTParam* map = GetTmap(cid, plid, seg, ud);
  if(!map) return qnan;
  return map->Offset();
}

//_____________________________________________________________________________
Double_t
HodoParamMan::GetGain(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const
{
  HodoTParam* map = GetTmap(cid, plid, seg, ud);
  if(!map) return qnan;
  return map->Gain();
}

//_____________________________________________________________________________
HodoTParam*
HodoParamMan::GetTmap(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const
{
  Int_t key = KEY(cid, plid, seg, ud);
  TIterator itr = m_TPContainer.find(key);
  if(itr != m_TPContainer.end())
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
HodoTParam*
HodoParamMan::GetT2map(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const
{
  Int_t key = KEY(cid, plid, seg, ud);
  TIterator itr = m_T2PContainer.find(key);
  if(itr != m_T2PContainer.end())
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
HodoAParam*
HodoParamMan::GetAHmap(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const
{
  Int_t key = KEY(cid, plid, seg, ud);
  AIterator itr = m_AHPContainer.find(key);
  if(itr != m_AHPContainer.end())
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
HodoAParam*
HodoParamMan::GetALmap(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const
{
  Int_t key = KEY(cid, plid, seg, ud);
  AIterator itr = m_ALPContainer.find(key);
  if(itr != m_ALPContainer.end())
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
HodoFParam*
HodoParamMan::GetFmap(Int_t cid, Int_t plid, Int_t seg, Int_t ud) const
{
  Int_t key = KEY(cid, plid, seg, ud);
  FIterator itr = m_FPContainer.find(key);
  if(itr != m_FPContainer.end())
    return itr->second;
  else
    return nullptr;
}
