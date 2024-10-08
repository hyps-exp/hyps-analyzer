// -*- C++ -*-

#include "BH2Filter.hh"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <iterator>

#include <TString.h>

#include <UnpackerConfig.hh>
#include <UnpackerXMLReadDigit.hh>

#include "BH2Hit.hh"
#include "DetectorID.hh"
#include "DCHit.hh"
#include "DCAnalyzer.hh"
#include "Exception.hh"
#include "FuncName.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"

//_____________________________________________________________________________
BH2Filter::Param::Param()
  : m_xmin(NumOfLayersBcOut+1),
    m_xmax(NumOfLayersBcOut+1)
{
}

//_____________________________________________________________________________
BH2Filter::Param::~Param()
{
}

//_____________________________________________________________________________
void
BH2Filter::Param::Print(const TString& arg) const
{
  const Int_t w = 3;
  std::cout << arg << ":";
  for(Int_t i=0, n=m_xmin.size(); i<n; ++i){
    if(i == 6){
      hddaq::cout << std::endl << std::setw(2) << " ";
    }
    std::cout << " (" << std::right
	      << std::setw(w) << m_xmin[i] << " "
    	      << std::setw(w) << m_xmax[i] << ")";
    // std::cout << " iplane " << std::setw(3) << i
    // 	      << " (" << std::setw(6) << m_xmin[i] << " "
    // 	      << std::setw(6) << m_xmax[i] << ")" << std::endl;
  }
  std::cout << std::endl;
}

//_____________________________________________________________________________
BH2Filter::BH2Filter()
  : m_is_ready(false),
    m_verbose(false),
    m_param(),
    m_dc(),
    m_hodo()
{
}

//_____________________________________________________________________________
BH2Filter::~BH2Filter()
{
}

//_____________________________________________________________________________
void
BH2Filter::Apply(const HodoAnalyzer& hodo, const DCAnalyzer& dc,
                 FilterList& cands)
{
  if(!m_is_ready) throw Exception(FUNC_NAME + " not initialized");
  m_dc = &dc;
  m_hodo = &hodo;
  std::set<Int_t> seg;
  for(Int_t i=0, n=hodo.GetNHits("BH2"); i<n; ++i){
    auto h = hodo.GetHit("BH2", i);
    if(!h) continue;
    seg.insert(h->SegmentId());
  }
  BuildCandidates(seg, cands);
}

//_____________________________________________________________________________
void
BH2Filter::Apply(Int_t T0Seg, const DCAnalyzer& dc, FilterList& cands)
{
  if(!m_is_ready) throw Exception(FUNC_NAME + " not initialized");
  m_dc = &dc;
  std::set<Int_t> seg;
  seg.insert(T0Seg);
  BuildCandidates(seg, cands);
}

//_____________________________________________________________________________
void
BH2Filter::BuildCandidates(std::set<Int_t>& seg, FilterList& cands)
{
  if(m_verbose) std::cout << FUNC_NAME << std::endl;
  cands.resize(seg.size());
  FIterator itCont = cands.begin();
  for(auto itSeg=seg.begin(), itSegEnd=seg.end();
      itSeg != itSegEnd; ++itSeg, ++itCont){
    const Int_t iSeg = *itSeg;
    auto& c = *itCont;
    c.resize(NumOfLayersBcOut+2);
    if(m_verbose) std::cout << "  BH2 seg = " << iSeg << std::endl;
    for(Int_t iplane=0; iplane<NumOfLayersBcOut+1; ++iplane){
      Int_t iLayer = iplane + 1;
      auto & after = c[iLayer];
      const Double_t xmin = m_param[iSeg].m_xmin[iplane];
      const Double_t xmax = m_param[iSeg].m_xmax[iplane];
      const auto& before = m_dc->GetBcOutHC(iLayer);
      for(Int_t ih=0, nh=before.size(); ih<nh; ++ih){
	auto h = before[ih];
	if(!h) continue;
	const Double_t wpos = h->GetWirePosition();
	if(m_verbose){
	  Int_t layer = h->GetLayer();
	  std::cout << " layer = " << std::setw(2) << iplane
		    << "(" << std::setw(3) << layer << ") : "
		    << std::setw(6) << wpos
		    << " (" << std::setw(6) << xmin
		    << ", " << std::setw(6) << xmax << ")";
	}
	if(wpos<xmin || xmax<wpos){
	  if(m_verbose) std::cout << std::endl;
	  continue;
	}
	if(m_verbose) std::cout << " good " << std::endl;
	after.push_back(const_cast<DCHit*>(h));
      }
    }
  }
}

//_____________________________________________________________________________
const std::vector<Double_t>&
BH2Filter::GetXmax(Int_t seg) const
{
  return m_param.at(seg).m_xmax;
}

//_____________________________________________________________________________
const std::vector<Double_t>&
BH2Filter::GetXmin(Int_t seg) const
{
  return m_param.at(seg).m_xmin;
}

//_____________________________________________________________________________
Bool_t
BH2Filter::Initialize(const TString& file_name)
{
  std::ifstream f(file_name);
  if(!f.is_open()){
    std::cerr << FUNC_NAME << " file open fail " << file_name << std::endl;
    return false;
  }

  TString line;
  while(f.good() && line.ReadLine(f)){
    if(line[0]=='#') continue;
    std::istringstream iss(line.Data());
    std::istream_iterator<Double_t> begin(iss);
    std::istream_iterator<Double_t> end;
    std::vector<Double_t> v(begin, end);
    if(v.size() < kNParam){
      hddaq::cerr << FUNC_NAME
		  << " number of parameters = " << v.size()
		  << ": required = " << kNParam << std::endl;
      continue;
    }

    const Int_t bh2Seg  = static_cast<Int_t>(v[kBH2Segment]);
    const Int_t bcPlane = static_cast<Int_t>(v[kLayerID]);
    const Double_t xmin = v[kXMin];
    const Double_t xmax = v[kXMax];
    Int_t iplane = bcPlane - PlMinBcOut;
    if(m_verbose){
      hddaq::cout << " seg = "    << std::setw(2) << bh2Seg
		  << ", plane = " << std::setw(3) << bcPlane
		  << "(" << iplane << ")"
		  << ", xmin = "  << std::setw(5) << xmin
		  << ", xmax = "  << std::setw(5) << xmax
		  << std::endl;
    }
    m_param[bh2Seg].m_xmin[iplane] = xmin;
    m_param[bh2Seg].m_xmax[iplane] = xmax;
  }

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
void
BH2Filter::Print(Option_t*) const
{
  std::cout << FUNC_NAME << std::endl
	    << "   (xmin xmax)" << std::endl;
  for(const auto& elem: m_param){
    elem.second.Print(std::to_string(elem.first));
  }
}
