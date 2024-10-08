// -*- C++ -*-

#include "RawData.hh"

#include <algorithm>
#include <iostream>
#include <string>

#include <TF1.h>

#include <std_ostream.hh>
#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "Exception.hh"
#include "FuncName.hh"
#include "HodoRawHit.hh"
#include "UserParamMan.hh"

namespace
{
using namespace hddaq::unpacker;
const auto& gUnpacker     = GUnpacker::get_instance();
const auto& gUser         = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
RawData::RawData()
  : m_is_decoded(),
    m_hodo_raw_hit_collection(),
    m_dc_raw_hit_collection()
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
RawData::~RawData()
{
  Clear();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
RawData::Clear(const TString& name)
{
  if(name.IsNull()){
    for(auto& elem: m_hodo_raw_hit_collection)
      del::ClearContainer(elem.second);
    for(auto& elem: m_dc_raw_hit_collection)
      del::ClearContainer(elem.second);
    m_hodo_raw_hit_collection.clear();
    m_dc_raw_hit_collection.clear();
  }else{
    del::ClearContainer(m_hodo_raw_hit_collection[name]);
    del::ClearContainer(m_dc_raw_hit_collection[name]);
  }
}

//_____________________________________________________________________________
Bool_t
RawData::DecodeHits(const TString& name)
{
  static const auto& digit_info = GConfig::get_instance().get_digit_info();

  if(m_is_decoded[name]){
    hddaq::cerr << FUNC_NAME << " " << name << " is already decoded"
                << std::endl;
    return false;
  }

  if(name.IsNull()){
    Bool_t ret = true;
    for(const auto& n: digit_info.get_name_list()){
      if(!n.empty())
        ret &= DecodeHits(n);
    }
    return ret;
  }

  auto id = digit_info.get_device_id(name.Data());
  const TString type = digit_info.get_device_type(id);

#if 0
  hddaq::cout << FUNC_NAME << std::endl
              << id << " " << name << " " << type << std::endl;
#endif

  if(type.IsNull())
    return false;

  Clear(name);

  Bool_t is_hodo  = type.Contains("Hodo", TString::kIgnoreCase);
  Bool_t is_fiber = type.Contains("Fiber", TString::kIgnoreCase);
  Bool_t is_dc    = type.Contains("DC", TString::kIgnoreCase);

  if(!is_hodo && !is_fiber && !is_dc)
    return false;

  for(Int_t plane=0, n_plane=gUnpacker.get_n_plane(id);
      plane<n_plane; ++plane){
    for(Int_t seg=0, n_seg=gUnpacker.get_n_segment(id, plane);
        seg<n_seg; ++seg){
      for(Int_t ch=0, n_ch=gUnpacker.get_n_ch(id, plane, seg);
          ch<n_ch; ++ch){
        for(Int_t data=0, n_data=gUnpacker.get_n_data(id, plane, seg, ch);
            data<n_data; ++data){
          for(Int_t i=0, n=gUnpacker.get_entries(id, plane, seg, ch, data);
              i<n; ++i){
            UInt_t val = gUnpacker.get(id, plane, seg, ch, data, i);
            if(is_hodo)  AddHodoRawHit(name, plane, seg, ch, data, val);
            if(is_fiber) AddFiberRawHit(name, plane, seg, ch, data, val);
            if(is_dc)    AddDCRawHit(name, plane, seg, ch, data, val);
          }
        }
      }
    }
  }

  m_is_decoded[name] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddHodoRawHit(const TString& name, Int_t plane, Int_t seg,
		       Int_t ch, Int_t data, Double_t val)
{
  auto& cont = m_hodo_raw_hit_collection[name];
  HodoRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    HodoRawHit* q = cont[i];
    if(true
       && q->DetectorName() == name
       && q->PlaneId() == plane
       && q->SegmentId() == seg){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit(name, plane, seg);
    cont.push_back(p);
  }

  if(data == gUnpacker.get_data_id(name, "adc")){
    p->SetAdc(ch, val);
  }else if(data == gUnpacker.get_data_id(name, "tdc")){
    p->SetTdcLeading(ch, val);
  }else if(data == gUnpacker.get_data_id(name, "trailing")){
    p->SetTdcTrailing(ch, val);
  // }else if(data == gUnpacker.get_data_id(name, "cstop")){
  //   ;
  }else if(data == gUnpacker.get_data_id(name, "overflow")){
    p->SetTdcOverflow(ch, val);
  }
  else{
    hddaq::cerr << FUNC_NAME << " wrong data type " << std::endl
                << " Detector   = " << name  << std::endl
                << " Plane      = " << plane << std::endl
                << " Segment    = " << seg   << std::endl
                << " Channel    = " << ch    << std::endl
                << " Data       = " << data  << std::endl;
    return false;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddFiberRawHit(const TString& name, Int_t plane, Int_t seg,
		       Int_t ch, Int_t data, Double_t val)
{
  auto& cont = m_hodo_raw_hit_collection[name];
  HodoRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    HodoRawHit* q = cont[i];
    if(true
       && q->DetectorName() == name
       && q->PlaneId() == plane
       && q->SegmentId() == seg){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit(name, plane, seg);
    cont.push_back(p);
  }

  if(data == gUnpacker.get_data_id(name, "leading")){
    p->SetTdcLeading(ch, val);
  }else if(data == gUnpacker.get_data_id(name, "trailing")){
    p->SetTdcTrailing(ch, val);
  }else if(data == gUnpacker.get_data_id(name, "highgain")){
    p->SetAdcHigh(ch, val);
  }else if(data == gUnpacker.get_data_id(name, "lowgain")){
    p->SetAdcLow(ch, val);
  }else{
    hddaq::cerr << FUNC_NAME << " wrong data type " << std::endl
                << " Detector   = " << name  << std::endl
                << " Plane      = " << plane << std::endl
                << " Segment    = " << seg   << std::endl
                << " Channel    = " << ch    << std::endl
                << " Data       = " << data  << std::endl;
    return false;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddDCRawHit(const TString& name, Int_t plane, Int_t seg,
                     Int_t ch, Int_t data, Double_t val)
{
  auto& cont = m_dc_raw_hit_collection[name];
  Int_t wire = ch;
  DCRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    DCRawHit* q = cont[i];
    if(true
       && q->DetectorName() == name
       && q->PlaneId() == plane
       && q->WireId() == wire){
      p=q; break;
    }
  }
  if(!p){
    p = new DCRawHit(name, plane, wire);
    cont.push_back(p);
  }

  if(data == gUnpacker.get_data_id(name, "leading")){
    p->SetTdc(val);
  }else if(data == gUnpacker.get_data_id(name, "trailing")){
    p->SetTrailing(val);
  }else if(data == gUnpacker.get_data_id(name, "overflow")){
    p->SetTdcOverflow(val);
  }else{
    hddaq::cerr << FUNC_NAME << " unknown data type " << std::endl
                << " Detector = " << name  << std::endl
		<< " PlaneId  = " << plane << std::endl
		<< " WireId   = " << wire  << std::endl
		<< " DataType = " << data  << std::endl
		<< " Value    = " << val   << std::endl;
  }
  return true;
}

//_____________________________________________________________________________
void
RawData::TdcCutBCOut()
{
  for(const auto& name: DCNameList.at("BcOut")){
    const Double_t MinTdc = gUser.GetParameter("Tdc"+name, 0);
    const Double_t MaxTdc = gUser.GetParameter("Tdc"+name, 1);
    TdcCut(name, MinTdc, MaxTdc);
  }
}

//_____________________________________________________________________________
void
RawData::TdcCutSDCIn()
{
  for(const auto& name: DCNameList.at("SdcIn")){
    const Double_t MinTdc = gUser.GetParameter("Tdc"+name, 0);
    const Double_t MaxTdc = gUser.GetParameter("Tdc"+name, 1);
    TdcCut(name, MinTdc, MaxTdc);
  }
}

//_____________________________________________________________________________
void
RawData::TdcCutSDCOut()
{
  for(const auto& name: DCNameList.at("SdcOut")){
    const Double_t MinTdc = gUser.GetParameter("Tdc"+name, 0);
    const Double_t MaxTdc = gUser.GetParameter("Tdc"+name, 1);
    TdcCut(name, MinTdc, MaxTdc);
  }
}

//_____________________________________________________________________________
void
RawData::TdcCut(const TString& name,
		Double_t min_tdc, Double_t max_tdc)
{
  DCRHC& HitCont = m_dc_raw_hit_collection.at(name);
  for(auto& hit: HitCont){
    hit->TdcCut(min_tdc, max_tdc);
  }
  EraseEmptyHits(HitCont);
}

//_____________________________________________________________________________
void
RawData::EraseEmptyHits(DCRHC& HitCont)
{
  auto i = HitCont.begin();
  while(i != HitCont.end()){
    if((*i)->IsEmpty()){
      delete *i;
      i = HitCont.erase(i);
    }
    else ++i;
  }
}

//_____________________________________________________________________________
const HodoRHC&
RawData::GetHodoRawHitContainer(const TString& name) const
{
  auto itr = m_hodo_raw_hit_collection.find(name);
  if(itr == m_hodo_raw_hit_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static HodoRHC null_container;
    return null_container;
  }else{
    return itr->second;
  }
}

//_____________________________________________________________________________
const DCRHC&
RawData::GetDCRawHitContainer(const TString& name) const
{
  auto itr = m_dc_raw_hit_collection.find(name);
  if(itr == m_dc_raw_hit_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static DCRHC null_container;
    return null_container;
  }else{
    return itr->second;
  }
}

//_____________________________________________________________________________
void
RawData::Print(Option_t*) const
{
  for(const auto& p: m_hodo_raw_hit_collection){
    for(const auto& hit: p.second){
      hit->Print();
    }
  }
  for(const auto& p: m_dc_raw_hit_collection){
    for(const auto& hit: p.second){
      hit->Print();
    }
  }
}
