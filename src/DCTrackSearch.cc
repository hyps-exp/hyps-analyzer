// -*- C++ -*-

#include "DCTrackSearch.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <TH2D.h>
#include <TH3D.h>

#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "DCLTrackHit.hh"
#include "DCPairHitCluster.hh"
#include "DCParameters.hh"
#include "CFTLocalTrack.hh"
#include "DebugTimer.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "HodoHit.hh"
#include "MathTools.hh"
#include "MWPCCluster.hh"
#include "TrackMaker.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"
#include "ConfMan.hh"

#include "RootHelper.hh"

namespace
{
const auto& gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& zTarget    = gGeom.LocalZ("Target");
const auto& zK18Target = gGeom.LocalZ("K18Target");
const auto& zBH2       = gGeom.LocalZ("BH2");
const Double_t MaxChisquare       = 2000.; // Set to be More than 30
const Double_t MaxChisquareSdcIn  = 5000.; // Set to be More than 30
const Double_t MaxChisquareCFT    = 500.;   // Set to be More than 30
const Double_t MaxChi2CFT2nd_phi = 30.;
const Double_t MaxChi2CFT2nd_uv  = 150.;

const Double_t MaxNumOfCluster = 20.;    // Set to be Less than 30
const Double_t MaxNumOfClusterCFT = 10.;    // Set to be Less than 30
const Double_t MaxCombi = 1.0e6;    // Set to be Less than 10^6
// SdcIn & BcOut for XUV Tracking routine
const Double_t MaxChisquareVXU = 50.;//
const Double_t ChisquareCutVXU = 50.;//

const Double_t Bh2SegX[NumOfSegBH2]      = {35./2., 10./2., 7./2., 7./2., 7./2., 7./2., 10./2., 35./2.};
const Double_t Bh2SegXAcc[NumOfSegBH2]   = {20., 6.5, 5., 5., 5., 5., 6.5, 20.};
const Double_t localPosBh2X_dX           = 0.;
const Double_t localPosBh2X[NumOfSegBH2] = {-41.5 + localPosBh2X_dX,
  -19.0 + localPosBh2X_dX,
  -10.5 + localPosBh2X_dX,
  -3.5  + localPosBh2X_dX,
  3.5   + localPosBh2X_dX,
  10.5  + localPosBh2X_dX,
  19.0  + localPosBh2X_dX,
  41.5  + localPosBh2X_dX};

//_____________________________________________________________________________
// Local Functions

//_____________________________________________________________________________
template <typename T> void
CalcTracks(std::vector<T*>& trackCont)
{
  for(auto& track: trackCont) track->Calculate();
}

//_____________________________________________________________________________
template <typename T> void
ClearFlags(std::vector<T*>& trackCont)
{
  for(const auto& track: trackCont){
    if(!track) continue;
    Int_t nh = track->GetNHit();
    for(Int_t j=0; j<nh; ++j) track->GetHit(j)->QuitTrack();
  }
}

//_____________________________________________________________________________
inline void
DeleteDuplicatedTracks(DCLocalTC& trackCont, Double_t ChisqrCut=0.)
{
  // evaluate container size in every iteration
  for(std::size_t i=0; i<trackCont.size(); ++i){
    const auto& tp = trackCont[i];
    if(!tp) continue;
    Int_t nh = tp->GetNHit();
    for(Int_t j=0; j<nh; ++j) tp->GetHit(j)->JoinTrack();
    for(std::size_t i2=trackCont.size()-1; i2>i; --i2){
      const DCLocalTrack* tp2 = trackCont[i2];
      Int_t nh2 = tp2->GetNHit(), flag=0;
      Double_t chisqr = tp2->GetChiSquare();
      for(Int_t j=0; j<nh2; ++j)
        if(tp2->GetHit(j)->BelongToTrack()) ++flag;
      if(flag>0 && chisqr>ChisqrCut){
        delete tp2;
        tp2 = 0;
        trackCont.erase(trackCont.begin()+i2);
      }
    }
  }
}

//_____________________________________________________________________________
inline void
DeleteDuplicatedTracks(DCLocalTC& trackCont,
                       Int_t first, Int_t second, Double_t ChisqrCut=0.)
{
  std::vector<Int_t> delete_index;
  // evaluate container size in every iteration
  for(std::size_t i=first; i<=second; ++i){

    auto itr = std::find(delete_index.begin(), delete_index.end(), i);
    if(itr != delete_index.end())
      continue;

    const DCLocalTrack* const tp = trackCont[i];
    if(!tp) continue;

    Int_t nh = tp->GetNHit();
    for(Int_t j=0; j<nh; ++j) tp->GetHit(j)->JoinTrack();

    for(std::size_t i2=second; i2>i; --i2){
      auto itr = std::find(delete_index.begin(), delete_index.end(), i2);
      if(itr != delete_index.end())
        continue;

      const DCLocalTrack* tp2 = trackCont[i2];
      Int_t nh2 = tp2->GetNHit(), flag=0;
      Double_t chisqr = tp2->GetChiSquare();
      for(Int_t j=0; j<nh2; ++j)
        if(tp2->GetHit(j)->BelongToTrack()) ++flag;
      if(flag>0 && chisqr>ChisqrCut){
        delete tp2;
        tp2 = 0;

        delete_index.push_back(i2);
      }
    }
  }

  // sort from bigger order
  std::sort(delete_index.begin(), delete_index.end(), std::greater<Int_t>());
  for(Int_t i=0; i<delete_index.size(); i++) {
    trackCont.erase(trackCont.begin()+delete_index[i]);
  }

  // reset hit record of DCHit
  for(std::size_t i=0; i<trackCont.size(); ++i){
    const DCLocalTrack* const tp = trackCont[i];
    if(!tp) continue;
    Int_t nh = tp->GetNHit();
    for(Int_t j=0; j<nh; ++j) tp->GetHit(j)->QuitTrack();
  }
}

//_____________________________________________________________________________
inline void // for SSD PreTrack
DeleteWideTracks(DCLocalTC& TrackContX,
                 DCLocalTC& TrackContY,
                 const std::size_t Nth=1)
{
  const Double_t xSize = 400./2*1.; // SDC1 X-Size
  const Double_t ySize = 250./2*1.; // SDC1 Y-Size
  const Double_t zPos  = -658.;  // SDC1 Z-Position "V1"
  // X
  if(TrackContX.size()>Nth){
    for(Int_t i=TrackContX.size()-1; i>=0; --i){
      const DCLocalTrack* tp = TrackContX[i];
      if(!tp) continue;
      Double_t x = tp->GetVXU(zPos);
      if(std::abs(x)>xSize){
        delete tp;
        tp = 0;
        TrackContX.erase(TrackContX.begin()+i);
      }
    }
  }
  // Y
  if(TrackContY.size()>Nth){
    for(Int_t i=TrackContY.size()-1; i>=0; --i){
      const DCLocalTrack* tp = TrackContY[i];
      if(!tp) continue;
      Double_t y = tp->GetVXU(zPos);
      if(std::abs(y)>ySize){
        delete tp;
        tp = 0;
        TrackContY.erase(TrackContY.begin()+i);
      }
    }
  }
}

//_____________________________________________________________________________
inline void // for SSD PreTrack
DeleteWideTracks(DCLocalTC& TrackContX,
                 const std::vector<DCHC>& SdcInHC,
                 const std::size_t Nth=1)
{
  const Double_t xSize =  100.;    // MaxDiff
  const Double_t zPos  = -635.467; // SDC1 Z-Position "X1/X2"
  // X
  if(TrackContX.size()>Nth){
    for(Int_t i=TrackContX.size()-1; i>=0; --i){
      const DCLocalTrack* tp = TrackContX[i];
      if(!tp) continue;
      Double_t x = tp->GetVXU(zPos);

      Bool_t accept = false;
      for(std::size_t j=0; j<NumOfLayersSDC1; ++j){
        if(j!=2 && j!=3) continue; // only X1/X2
        const std::size_t nh = SdcInHC[j+1].size();
        for(std::size_t k=0; k<nh; ++k){
          const DCHit *hit = SdcInHC[j+1][k];
          if(!hit) continue;
          if(std::abs(x-hit->GetWirePosition()) < xSize){
            accept = true;
          }
        }
      }
      if(!accept){
        delete tp;
        tp = 0;
        TrackContX.erase(TrackContX.begin()+i);
      }
    }
  }

  return;
}

//_____________________________________________________________________________
[[maybe_unused]] void
DebugPrint(const IndexList& nCombi,
           const TString& func_name="",
           const TString& msg="")
{
  Int_t n  =1;
  Int_t nn =1;
  Int_t sum=0;
  hddaq::cout << func_name << ":" ;
  IndexList::const_iterator itr, end = nCombi.end();
  for(itr=nCombi.begin(); itr!=end; ++itr){
    Int_t val = *itr;
    sum += val;
    nn *= (val+1);
    hddaq::cout << " " << val;
    if(val!=0) n *= val;
  }
  if(sum==0)
    n=0;
  hddaq::cout << ": total = " << n << ", " << nn << ", " << std::endl;
  return;
}

//_____________________________________________________________________________
[[maybe_unused]] void
DebugPrint(const DCLocalTC& trackCont,
           const TString& arg="")
{
  const Int_t nn = trackCont.size();
  hddaq::cout << arg << " " << nn << std::endl;
  for(Int_t i=0; i<nn; ++i){
    const DCLocalTrack * const track=trackCont[i];
    if(!track) continue;
    Int_t    nh     = track->GetNHit();
    Double_t chisqr = track->GetChiSquare();
    hddaq::cout << std::setw(4) << i
                << "  #Hits : " << std::setw(2) << nh
                << "  ChiSqr : " << chisqr << std::endl;
  }
  hddaq::cout << std::endl;
}

//_____________________________________________________________________________
[[maybe_unused]] void
DebugPrint(const IndexList& nCombi,
           const std::vector<ClusterList>& CandCont,
           const TString& arg="")
{
  hddaq::cout << arg << " #Hits of each group" << std::endl;
  Int_t np = nCombi.size();
  Int_t nn = 1;
  for(Int_t i=0; i<np; ++i){
    hddaq::cout << std::setw(4) << nCombi[i];
    nn *= nCombi[i] + 1;
  }
  hddaq::cout << " -> " << nn-1 << " Combinations" << std::endl;
  for(Int_t i=0; i<np; ++i){
    Int_t n=CandCont[i].size();
    hddaq::cout << "[" << std::setw(3) << i << "]: "
                << std::setw(3) << n << " ";
    for(Int_t j=0; j<n; ++j){
      hddaq::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire()
                  << "(" << CandCont[i][j]->NumberOfHits() << ")"
                  << " ";
    }
    hddaq::cout << std::endl;
  }
}

//_____________________________________________________________________________
template <class Functor>
inline void
FinalizeTrack(const TString& arg,
              DCLocalTC& trackCont,
              Functor comp,
              std::vector<ClusterList>& candCont,
              Bool_t delete_flag=true)
{
  ClearFlags(trackCont);

#if 0
  DebugPrint(trackCont, arg+" Before Sorting ");
#endif

  std::stable_sort(trackCont.begin(), trackCont.end(), DCLTrackComp_Nhit());


#if 0
  DebugPrint(trackCont, arg+" After Sorting (Nhit) ");
#endif

  typedef std::pair<Int_t, Int_t> index_pair;
  std::vector<index_pair> index_pair_vec;

  std::vector<Int_t> nhit_vec;

  for(Int_t i=0; i<trackCont.size(); i++) {
    Int_t nhit = trackCont[i]->GetNHit();
    nhit_vec.push_back(nhit);
  }

  if(!nhit_vec.empty()) {
    Int_t max_nhit = nhit_vec.front();
    Int_t min_nhit = nhit_vec.back();
    for(Int_t nhit=max_nhit; nhit>=min_nhit; nhit--) {
      auto itr1 = std::find(nhit_vec.begin(), nhit_vec.end(), nhit);
      if(itr1 == nhit_vec.end())
        continue;

      size_t index1 = std::distance(nhit_vec.begin(), itr1);

      auto itr2 = std::find(nhit_vec.rbegin(), nhit_vec.rend(), nhit);
      size_t index2 = nhit_vec.size() - std::distance(nhit_vec.rbegin(), itr2) - 1;

      index_pair_vec.push_back(index_pair(index1, index2));
    }
  }

  for(Int_t i=0; i<index_pair_vec.size(); i++) {
    std::stable_sort(trackCont.begin() + index_pair_vec[i].first,
                     trackCont.begin() +  index_pair_vec[i].second + 1, DCLTrackComp_Chisqr());
  }

#if 0
  DebugPrint(trackCont, arg+" After Sorting (chisqr)");
#endif

  if(delete_flag) {
    for(Int_t i = index_pair_vec.size()-1; i>=0; --i) {
      DeleteDuplicatedTracks(trackCont, index_pair_vec[i].first, index_pair_vec[i].second, 0.);
    }
  }

#if 0
  DebugPrint(trackCont, arg+" After Deleting in each hit number");
#endif


  std::stable_sort(trackCont.begin(), trackCont.end(), comp);

#if 0
  DebugPrint(trackCont, arg+" After Sorting with comp func ");
#endif

  if(delete_flag) DeleteDuplicatedTracks(trackCont);

#if 0
  DebugPrint(trackCont, arg+" After Deleting ");
#endif

  CalcTracks(trackCont);
  del::ClearContainerAll(candCont);
}

//_____________________________________________________________________________
// MakeCluster _________________________________________________________________

//_____________________________________________________________________________
Bool_t
MakePairPlaneHitCluster(const DCHC & HC1,
                        const DCHC & HC2,
                        Double_t CellSize,
                        ClusterList& Cont,
                        Bool_t honeycomb=false)
{
  Int_t nh1=HC1.size(), nh2=HC2.size();
  std::vector<Int_t> UsedFlag(nh2,0);
  for(Int_t i1=0; i1<nh1; ++i1){
    DCHit *hit1=HC1[i1];
    Double_t wp1=hit1->GetWirePosition();
    Bool_t flag=false;
    for(Int_t i2=0; i2<nh2; ++i2){
      DCHit *hit2=HC2[i2];
      Double_t wp2=hit2->GetWirePosition();
      if(std::abs(wp1-wp2)<CellSize){
        Int_t multi1 = hit1->GetEntries();
        Int_t multi2 = hit2->GetEntries();
        for(Int_t m1=0; m1<multi1; ++m1) {
          if(!hit1->IsGood(m1))
            continue;
          for(Int_t m2=0; m2<multi2; ++m2) {
            if(!hit2->IsGood(m2))
              continue;
            Double_t x1,x2;
            if(wp1<wp2){
              x1=wp1+hit1->DriftLength(m1);
              x2=wp2-hit2->DriftLength(m2);
            }
            else {
              x1=wp1-hit1->DriftLength(m1);
              x2=wp2+hit2->DriftLength(m2);
            }
            DCPairHitCluster *cluster =
              new DCPairHitCluster(new DCLTrackHit(hit1,x1,m1),
                                   new DCLTrackHit(hit2,x2,m2));
            cluster->SetHoneycomb(honeycomb);
            Cont.push_back(cluster);
            flag=true; ++UsedFlag[i2];
          }
        }
      }
    }
#if 1
    if(!flag){
      Int_t multi1 = hit1->GetEntries();
      for(Int_t m1=0; m1<multi1; m1++) {
        if(!(hit1->IsGood(m1))) continue;
        Double_t dl=hit1->DriftLength(m1);
        DCPairHitCluster *cluster1 = new DCPairHitCluster(new DCLTrackHit(hit1,wp1+dl,m1));
        DCPairHitCluster *cluster2 = new DCPairHitCluster(new DCLTrackHit(hit1,wp1-dl,m1));
        cluster1->SetHoneycomb(honeycomb);
        cluster2->SetHoneycomb(honeycomb);
        Cont.push_back(cluster1);
        Cont.push_back(cluster2);
      }
    }
#endif
  }
#if 1
  for(Int_t i2=0; i2<nh2; ++i2){
    if(UsedFlag[i2]==0) {
      DCHit *hit2=HC2[i2];
      Int_t multi2 = hit2->GetEntries();
      for(Int_t m2=0; m2<multi2; m2++) {
        if(!(hit2->IsGood(m2))) continue;
        Double_t wp=hit2->GetWirePosition();
        Double_t dl=hit2->DriftLength(m2);
        DCPairHitCluster *cluster1 = new DCPairHitCluster(new DCLTrackHit(hit2,wp+dl,m2));
        DCPairHitCluster *cluster2 = new DCPairHitCluster(new DCLTrackHit(hit2,wp-dl,m2));
        cluster1->SetHoneycomb(honeycomb);
        cluster2->SetHoneycomb(honeycomb);
        Cont.push_back(cluster1);
        Cont.push_back(cluster2);
      }
    }
  }
#endif
  return true;
}

//_____________________________________________________________________________
Bool_t
MakeUnPairPlaneHitCluster(const DCHC& HC,
                          ClusterList& Cont,
                          Bool_t honeycomb=false)
{
  const std::size_t nh = HC.size();
  for(std::size_t i=0; i<nh; ++i){
    DCHit *hit = HC[i];
    if(!hit) continue;
    std::size_t mh = hit->GetEntries();
    for(std::size_t m=0; m<mh; ++m) {
      if(!hit->IsGood(m)) continue;
      Double_t wp = hit->GetWirePosition();
      Double_t dl = hit->DriftLength(m);
      DCPairHitCluster *cluster1 =
        new DCPairHitCluster(new DCLTrackHit(hit,wp+dl,m));
      DCPairHitCluster *cluster2 =
        new DCPairHitCluster(new DCLTrackHit(hit,wp-dl,m));
      cluster1->SetHoneycomb(honeycomb);
      cluster2->SetHoneycomb(honeycomb);
      Cont.push_back(cluster1);
      Cont.push_back(cluster2);
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
MakeMWPCPairPlaneHitCluster(const DCHC& HC,
                            ClusterList& Cont)
{
  Int_t nh=HC.size();
  for(Int_t i=0; i<nh; ++i){
    DCHit *hit=HC[i];
    if(hit){
      Int_t multi = hit->GetEntries();
      for(Int_t m=0; m<multi; m++) {
        if(!(hit->IsGood(m))) continue;
        Double_t wp=hit->GetWirePosition();
        // Double_t dl=hit->DriftLength(m);
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit, wp, m)));
      }
    }
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
MakeTOFHitCluster(const DCHC& HitCont,
                  ClusterList& Cont,
                  Int_t xy)
{
  Int_t nh = HitCont.size();
  for(Int_t i=0; i<nh; ++i){
    if(i%2!=xy) continue;
    DCHit *hit = HitCont[i];
    if(!hit) continue;
    Double_t wp = hit->GetWirePosition();
    Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit, wp, 0)));
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
MakePairPlaneHitClusterVUX(const DCHC& HC1,
                           const DCHC& HC2,
                           Double_t CellSize,
                           ClusterList& Cont,
                           Bool_t honeycomb=false)
{
  Int_t nh1=HC1.size(), nh2=HC2.size();
  std::vector<Int_t> UsedFlag(nh2,0);

  for(Int_t i1=0; i1<nh1; ++i1){
    DCHit *hit1=HC1[i1];

    Double_t wp1=hit1->GetWirePosition();
    Bool_t flag=false;
    for(Int_t i2=0; i2<nh2; ++i2){
      DCHit *hit2=HC2[i2];
      Double_t wp2=hit2->GetWirePosition();
      if(std::abs(wp1-wp2)<CellSize){

        Int_t multi1 = hit1->GetEntries();
        Int_t multi2 = hit2->GetEntries();
        for(Int_t m1=0; m1<multi1; m1++) {
          if(!(hit1->IsGood(m1))) continue;
          for(Int_t m2=0; m2<multi2; m2++) {
            if(!(hit2->IsGood(m2))) continue;
            Double_t dl1=hit1->DriftLength(m1);
            Double_t dl2=hit2->DriftLength(m2);

            Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1+dl1,m1),
                                                new DCLTrackHit(hit2,wp2+dl2,m2)));
            Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1+dl1,m1),
                                                new DCLTrackHit(hit2,wp2-dl2,m2)));
            Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1-dl1,m1),
                                                new DCLTrackHit(hit2,wp2+dl2,m2)));
            Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1-dl1,m1),
                                                new DCLTrackHit(hit2,wp2-dl2,m2)));

            flag=true; ++UsedFlag[i2];
          }
        }
      }
    }
    if(!flag){
      Int_t multi1 = hit1->GetEntries();
      for(Int_t m1=0; m1<multi1; m1++) {
        if(!(hit1->IsGood(m1))) continue;
        Double_t dl=hit1->DriftLength(m1);
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1+dl,m1)));
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1-dl,m1)));
      }
    }
  }
  for(Int_t i2=0; i2<nh2; ++i2){
    if(UsedFlag[i2]==0) {
      DCHit *hit2=HC2[i2];
      Int_t multi2 = hit2->GetEntries();
      for(Int_t m2=0; m2<multi2; m2++) {
        if(!(hit2->IsGood(m2))) continue;

        Double_t wp=hit2->GetWirePosition();
        Double_t dl=hit2->DriftLength(m2);
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit2,wp+dl,m2)));
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit2,wp-dl,m2)));
      }
    }
  }

  return true;
}

}


//_____________________________________________________________________________
namespace track
{

//_____________________________________________________________________________
std::vector<IndexList>
MakeIndex(Int_t ndim, const Int_t *index1, Bool_t& status)
{
  if(ndim==1){
    std::vector<IndexList> index2;
    for(Int_t i=-1; i<index1[0]; ++i){
      IndexList elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector<IndexList> index2 = MakeIndex(ndim-1, index1+1, status);

  std::vector<IndexList> index;
  Int_t n2=index2.size();
  for(Int_t j=0; j<n2; ++j){
    for(Int_t i=-1; i<index1[0]; ++i){
      IndexList elem;
      Int_t n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      for(Int_t k=0; k<n3; ++k)
        elem.push_back(index2[j][k]);
      index.push_back(elem);
      Int_t size1=index.size();
      if(size1>MaxCombi){
        status = false;
#if 1
        hddaq::cout << FUNC_NAME << " too much combinations... " << n2 << std::endl;
#endif
        return std::vector<IndexList>(0);
      }
    }
  }

  return index;
}

//_____________________________________________________________________________
std::vector<IndexList>
MakeIndex(Int_t ndim, const IndexList& index1, Bool_t& status)
{
  return MakeIndex(ndim, &(index1[0]), status);
}

//_____________________________________________________________________________
std::vector<IndexList>
MakeIndex_VXU(Int_t ndim,Int_t maximumHit, const Int_t *index1)
{
  if(ndim==1){
    std::vector<IndexList> index2;
    for(Int_t i=-1; i<index1[0]; ++i){
      IndexList elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector<IndexList>
    index2=MakeIndex_VXU(ndim-1, maximumHit, index1+1);

  std::vector<IndexList> index;
  Int_t n2=index2.size();
  for(Int_t j=0; j<n2; ++j){
    for(Int_t i=-1; i<index1[0]; ++i){
      IndexList elem;
      Int_t validHitNum=0;
      Int_t n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      if(i != -1)
        validHitNum++;
      for(Int_t k=0; k<n3; ++k){
        elem.push_back(index2[j][k]);
        if(index2[j][k] != -1)
          validHitNum++;
      }
      if(validHitNum <= maximumHit)
        index.push_back(elem);
      // Int_t size1=index.size();
    }
  }

  return index;
}

//_____________________________________________________________________________
std::vector<IndexList>
MakeIndex_VXU(Int_t ndim,Int_t maximumHit, const IndexList& index1)
{
  return MakeIndex_VXU(ndim, maximumHit, &(index1[0]));
}

//_____________________________________________________________________________
DCLocalTrack*
MakeTrack(const std::vector<ClusterList>& CandCont,
          const IndexList& combination)
{
  DCLocalTrack *tp = new DCLocalTrack;
  for(std::size_t i=0, n=CandCont.size(); i<n; ++i){
    Int_t m = combination[i];
    if(m<0) continue;
    DCPairHitCluster *cluster = CandCont[i][m];
    if(!cluster) continue;
    Int_t mm = cluster->NumberOfHits();
    for(Int_t j=0; j<mm; ++j){
      DCLTrackHit *hitp = cluster->GetHit(j);
      if(!hitp) continue;
      tp->AddHit(hitp);
    }
#if 0
    hddaq::cout << FUNC_NAME << ":" << std::setw(3)
                << i << std::setw(3) << m  << " "
                << CandCont[i][m] << " " << mm << std::endl;
#endif
  }
  return tp;
}

//______________________________________________________________________________
CFTLocalTrack
* MakeTrack(  const std::vector<CFTFiberCC>&  HC,
	      const IndexList& combination)
{
  static const std::string funcname = "[MakeTrack]";

  CFTLocalTrack *tp = new CFTLocalTrack();

  if(!tp){
    std::cerr << funcname << ": new fail" << std::endl;
    return 0;
  }

  for( std::size_t i=0, n=HC.size(); i<n; ++i ){
    Int_t m=combination[i];
    if (m>=0) {
      CFTFiberCluster *fcl = HC[i][m];

#if 0
      std::cout << funcname << ":" << std::setw(3)
		<< i << std::setw(3) << m  << " "
		<< CandCont[i][m] << std::endl;
#endif

      if(fcl){
	if (i==CFT_PHI1 || i==CFT_PHI2 || i==CFT_PHI3 || i==CFT_PHI4)
	  tp->AddHit( fcl );
	else if (i==CFT_U1 || i==CFT_V2 || i==CFT_U3 || i==CFT_V4)
	  tp->AddHitU( fcl );
      }
    }
  }

  return tp;
}

//_____________________________________________________________________________
Int_t /* Local Track Search without BH2Filter */
LocalTrackSearch(const std::vector<DCHC>& HC,
                 const DCPairPlaneInfo * PpInfo,
                 Int_t npp, DCLocalTC& TrackCont,
                 Int_t MinNumOfHits, Int_t T0Seg)
{
  std::vector<ClusterList> CandCont(npp);

  for(Int_t i=0; i<npp; ++i){
    Bool_t ppFlag    = PpInfo[i].pair;
    Bool_t honeycomb = PpInfo[i].honeycomb;
    Bool_t fiber     = PpInfo[i].fiber;
    Int_t  layer1    = PpInfo[i].id1;
    Int_t  layer2    = PpInfo[i].id2;
    if(ppFlag && !fiber){
      MakePairPlaneHitCluster(HC[layer1], HC[layer2],
                              PpInfo[i].CellSize, CandCont[i], honeycomb);
    }else if(!ppFlag && fiber){
      PpInfo[i].Print(FUNC_NAME+" invalid parameter", hddaq::cerr);
    }else{
      MakeUnPairPlaneHitCluster(HC[layer1], CandCont[i], honeycomb);
    }
  }

  IndexList nCombi(npp);
  for(Int_t i=0; i<npp; ++i){
    Int_t n = CandCont[i].size();
    nCombi[i] = n>MaxNumOfCluster ? 0 : n;
  }

#if 0
  DebugPrint(nCombi, CandCont, FUNC_NAME);
#endif

  Bool_t status = true;
  std::vector<IndexList> CombiIndex = MakeIndex(npp, nCombi, status);

  for(Int_t i=0, n=CombiIndex.size(); i<n; ++i){
    DCLocalTrack *track = MakeTrack(CandCont, CombiIndex[i]);
    if(!track) continue;
    if(track->GetNHit()>=MinNumOfHits
       && track->DoFit()
       && track->GetChiSquare()<MaxChisquare){

      if(T0Seg>=0 && T0Seg<NumOfSegBH2) {
        Double_t xbh2=track->GetX(zBH2), ybh2=track->GetY(zBH2);
        Double_t difPosBh2 = localPosBh2X[T0Seg] - xbh2;

        //   Double_t xtgt=track->GetX(zTarget), ytgt=track->GetY(zTarget);
        //   Double_t ytgt=track->GetY(zTarget);

        if(true
           && fabs(difPosBh2)<Bh2SegXAcc[T0Seg]
           && (-10 < ybh2 && ybh2 < 40)
           //       && fabs(ytgt)<21.
          ){
          TrackCont.push_back(track);
        }else{
          delete track;
        }

      }else{
        TrackCont.push_back(track);
      }
    }
    else{
      delete track;
    }
  }

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackComp(), CandCont);
  return status? TrackCont.size() : -1;
}

//_____________________________________________________________________________
Int_t /* Local Track Search with BH2Filter */
LocalTrackSearch(const std::vector<std::vector<DCHC>> &hcAssemble,
                 const DCPairPlaneInfo * PpInfo,
                 Int_t npp, DCLocalTC &trackCont,
                 Int_t MinNumOfHits, Int_t T0Seg)
{
  std::vector<std::vector<DCHC>>::const_iterator
    itr, itr_end = hcAssemble.end();

  Int_t status = 0;
  for(itr=hcAssemble.begin(); itr!=itr_end; ++itr){
    const std::vector<DCHC>& l = *itr;
    DCLocalTC tc;
    status = LocalTrackSearch(l, PpInfo, npp, tc, MinNumOfHits, T0Seg);
    trackCont.insert(trackCont.end(), tc.begin(), tc.end());
  }

  ClearFlags(trackCont);
  std::stable_sort(trackCont.begin(), trackCont.end(), DCLTrackComp());

  DeleteDuplicatedTracks(trackCont);
  //    CalcTracks(trackCont);

  return status < 0? status : trackCont.size();
}

//_____________________________________________________________________________
Int_t
LocalTrackSearchSdcOut(const std::vector<DCHC>& SdcOutHC,
                       const DCPairPlaneInfo* PpInfo,
                       Int_t npp, DCLocalTC& TrackCont,
                       Int_t MinNumOfHits /*=6*/)
{
  std::vector<ClusterList> CandCont(npp);

  for(Int_t i=0; i<npp; ++i){
    Bool_t ppFlag    = PpInfo[i].pair;
    Bool_t honeycomb = PpInfo[i].honeycomb;
    Int_t  layer1    = PpInfo[i].id1;
    Int_t  layer2    = PpInfo[i].id2;
    if(ppFlag){
      MakePairPlaneHitCluster(SdcOutHC[layer1], SdcOutHC[layer2],
                              PpInfo[i].CellSize, CandCont[i], honeycomb);
    }else{
      MakeUnPairPlaneHitCluster(SdcOutHC[layer1], CandCont[i], honeycomb);
    }
  }

  IndexList nCombi(npp);
  for(Int_t i=0; i<npp; ++i){
    Int_t n = CandCont[i].size();
    nCombi[i] = n>MaxNumOfCluster ? 0 : n;
  }

  Bool_t status = true;
  std::vector<IndexList> CombiIndex = MakeIndex(npp, nCombi, status);

#if 0
  DebugPrint(nCombi, CandCont, FUNC_NAME);
#endif

  for(Int_t i=0, n=CombiIndex.size(); i<n; ++i) {
    DCLocalTrack *track = MakeTrack(CandCont, CombiIndex[i]);
    if(!track) continue;
    if(track->GetNHit()>=MinNumOfHits     &&
       track->GetNHitY() >= 2             &&
       track->DoFit()                     &&
       track->GetChiSquare()<MaxChisquare)
    {
      TrackCont.push_back(track);
    }
    else
    {
      delete track;
    }
  }

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackComp(), CandCont);
  return status? TrackCont.size() : -1;
}

//_____________________________________________________________________________
Int_t
LocalTrackSearchSdcOut(const DCHC& TOFHC,
                       const std::vector<DCHC>& SdcOutHC,
                       const DCPairPlaneInfo* PpInfo,
                       Int_t npp,
                       DCLocalTC& TrackCont,
                       Int_t MinNumOfHits /*=6*/)
{
  std::vector<ClusterList> CandCont(npp);

  for(Int_t i=0; i<npp-2; ++i){
    Bool_t ppFlag    = PpInfo[i].pair;
    Bool_t honeycomb = PpInfo[i].honeycomb;
    Int_t  layer1    = PpInfo[i].id1;
    Int_t  layer2    = PpInfo[i].id2;
    if(ppFlag){
      MakePairPlaneHitCluster(SdcOutHC[layer1], SdcOutHC[layer2],
                              PpInfo[i].CellSize, CandCont[i], honeycomb);
    }else {
      MakeUnPairPlaneHitCluster(SdcOutHC[layer1], CandCont[i], honeycomb);
    }
  }

  // TOF
  MakeTOFHitCluster(TOFHC, CandCont[npp-2], 0);
  MakeTOFHitCluster(TOFHC, CandCont[npp-1], 1);

  IndexList nCombi(npp);
  for(Int_t i=0; i<npp; ++i){
    Int_t n = CandCont[i].size();
    nCombi[i] = n>MaxNumOfCluster ? 0 : n;
  }

  Bool_t status = true;
  std::vector<IndexList> CombiIndex = MakeIndex(npp, nCombi, status);

#if 0
  DebugPrint(nCombi, CandCont, FUNC_NAME);
#endif

  for(Int_t i=0, n=CombiIndex.size(); i<n; ++i){
    DCLocalTrack *track = MakeTrack(CandCont, CombiIndex[i]);
    if(!track) continue;
    /*
    static const Int_t IdTOF_UX = gGeom.GetDetectorId("TOF-UX");
    static const Int_t IdTOF_UY = gGeom.GetDetectorId("TOF-UY");
    static const Int_t IdTOF_DX = gGeom.GetDetectorId("TOF-DX");
    static const Int_t IdTOF_DY = gGeom.GetDetectorId("TOF-DY");
    */
    static const Int_t IdTOF_X = gGeom.GetDetectorId("TOF-X");
    static const Int_t IdTOF_Y = gGeom.GetDetectorId("TOF-Y");
    static const Int_t IdLTOF_X_TL = gGeom.GetDetectorId("LEPS-TOF-X-Tilt-L");
    static const Int_t IdLTOF_Y_TL = gGeom.GetDetectorId("LEPS-TOF-Y-Tilt-L");
    static const Int_t IdLTOF_X_TR = gGeom.GetDetectorId("LEPS-TOF-X-Tilt-R");
    static const Int_t IdLTOF_Y_TR = gGeom.GetDetectorId("LEPS-TOF-Y-Tilt-R");
    /*
    Bool_t TOFSegXYMatching =
      (track->GetWire(IdTOF_UX)==track->GetWire(IdTOF_UY)) ||
      (track->GetWire(IdTOF_DX)==track->GetWire(IdTOF_DY));
    */
    
    Int_t TOFmulti=0;
    if(track->GetWire(IdTOF_X)==track->GetWire(IdTOF_Y)) TOFmulti++;
    if(track->GetWire(IdLTOF_X_TL)==track->GetWire(IdLTOF_Y_TL)) TOFmulti++;
    if(track->GetWire(IdLTOF_X_TR)==track->GetWire(IdLTOF_Y_TR)) TOFmulti++;
    Bool_t TOFSegXYMatching=false;
    if(TOFmulti==1) TOFSegXYMatching=true;
    
    //Bool_t TOFSegXYMatching = (track->GetWire(IdTOF_X) == track->GetWire(IdTOF_Y)) ||(track->GetWire(IdLTOF_X_TL) == track->GetWire(IdLTOF_Y_TL)) || (track->GetWire(IdLTOF_X_TR) == track->GetWire(IdLTOF_Y_TR));

    // Int_t Track[20]={0};
    // Int_t layer;
    // for(Int_t i=0; i<(track->GetNHit()); ++i){
    //  layer=track->GetHit(i)->GetLayer();
    //  Track[layer]=1;
    // }

    // Bool_t FBT =
    //  (Track[80]==1 && Track[82]==1) || (Track[81]==1 && Track[83]==1) ||
    //  (Track[84]==1 && Track[86]==1) || (Track[85]==1 && Track[87]==1) ;
    // Bool_t DC23x_off =
    //  (Track[31]==0 && Track[32]==0 && Track[37]==0 && Track[38]==0);


    if(TOFSegXYMatching &&
       track->GetNHit()>=MinNumOfHits+2   &&
       track->GetNHitY() >= 2             &&
       track->DoFit()                     &&
       track->GetChiSquare()<MaxChisquare)
    {
      TrackCont.push_back(track);
    }
    else
    {
      delete track;
    }
  }

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackComp(), CandCont);
  return status? TrackCont.size() : -1;
}

//_____________________________________________________________________________
Int_t /* Local Track Search SdcIn w/Fiber */
LocalTrackSearchSdcInFiber(const std::vector<DCHC>& HC,
                           const DCPairPlaneInfo* PpInfo,
                           Int_t npp, DCLocalTC& TrackCont,
                           Int_t MinNumOfHits /*=6*/)
{
  std::vector<ClusterList> CandCont(npp);

  for(Int_t i=0; i<npp; ++i){

    Bool_t ppFlag    = PpInfo[i].pair;
    Bool_t honeycomb = PpInfo[i].honeycomb;
    Int_t  layer1    = PpInfo[i].id1;
    Int_t  layer2    = PpInfo[i].id2;

    if(ppFlag) {
      MakePairPlaneHitCluster(HC[layer1], HC[layer2],
                              PpInfo[i].CellSize, CandCont[i], honeycomb);
    }else if(layer1==layer2){
      MakeMWPCPairPlaneHitCluster(HC[layer1], CandCont[i]);
    }else{
      MakeUnPairPlaneHitCluster(HC[layer1], CandCont[i], honeycomb);
    }
  }

  IndexList nCombi(npp);
  for(Int_t i=0; i<npp; ++i) {
    Int_t n = CandCont[i].size();
    nCombi[i] = n>MaxNumOfCluster ? 0 : n;
  }

  Bool_t status = true;
  std::vector<IndexList> CombiIndex = MakeIndex(npp, nCombi, status);

  for(Int_t i=0, n=CombiIndex.size(); i<n; ++i){
    DCLocalTrack *track = MakeTrack(CandCont, CombiIndex[i]);
    if(!track) continue;
    if(true
       && track->GetNHitSFT() > 1
       && track->GetNHit()>=MinNumOfHits
       && track->DoFit()
       && track->GetChiSquare()<MaxChisquare
      ){
      TrackCont.push_back(track);
    }
    else
      delete track;
  }

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackCompSdcInFiber(), CandCont);
  return status? TrackCont.size() : -1;
}

// BC3&4, SDC1 VUX Tracking ___________________________________________
Int_t
LocalTrackSearchVUX(const std::vector<DCHC>& HC,
                    const DCPairPlaneInfo* PpInfo,
                    Int_t npp, DCLocalTC& TrackCont,
                    Int_t MinNumOfHits /*=6*/)
{
  DCLocalTC TrackContV;
  DCLocalTC TrackContX;
  DCLocalTC TrackContU;
  std::vector<ClusterList>   CandCont(npp);
  std::vector<ClusterList>   CandContV(npp);
  std::vector<ClusterList>   CandContX(npp);
  std::vector<ClusterList>   CandContU(npp);

  //  Int_t NumOfLayersDC_12 = 12 ;
  Int_t iV=0, iX=0, iU=0;
  Int_t nV=0, nX=0, nU=0;
  for(Int_t i=0; i<npp; ++i){
    Bool_t ppFlag    = PpInfo[i].pair;
    Bool_t honeycomb = PpInfo[i].honeycomb;
    Int_t  layer1    = PpInfo[i].id1;
    Int_t  layer2    = PpInfo[i].id2;
    Int_t  nh1       = HC[layer1].size();
    Int_t  nh2       = HC[layer2].size();
    Double_t TiltAngle = 0.;
    if(nh1>0) TiltAngle = HC[layer1][0]->GetTiltAngle();
    if(ppFlag && nh1==0 && nh2>0)
      TiltAngle = HC[layer2][0]->GetTiltAngle();
    if(ppFlag && nh1>0 && nh2>0){
      if(TiltAngle<0){
        MakePairPlaneHitClusterVUX(HC[layer1], HC[layer2],
                                   PpInfo[i].CellSize, CandContV[iV],
                                   honeycomb);
        ++iV; nV = nV+2;
      }
      if(TiltAngle==0){
        MakePairPlaneHitClusterVUX(HC[layer1], HC[layer2],
                                   PpInfo[i].CellSize, CandContX[iX],
                                   honeycomb);
        ++iX; nX = nX+2;
      }
      if(TiltAngle>0){
        MakePairPlaneHitClusterVUX(HC[layer1], HC[layer2],
                                   PpInfo[i].CellSize, CandContU[iU],
                                   honeycomb);
        ++iU; nU = nU+2;
      }
    }
    if(!ppFlag){
      if(TiltAngle<0){
        MakeUnPairPlaneHitCluster(HC[layer1], CandContV[iV], honeycomb);
        ++nV; ++iV;
      }
      if(TiltAngle==0){
        MakeUnPairPlaneHitCluster(HC[layer1], CandContX[iX], honeycomb);
        ++nX; ++iX;
      }
      if(TiltAngle>0){
        MakeUnPairPlaneHitCluster(HC[layer1], CandContU[iU], honeycomb);
        ++nU; ++iU;
      }
    }
    if(ppFlag && nh1==0 && nh2>0){
      if(TiltAngle<0){
        MakeUnPairPlaneHitCluster(HC[layer2], CandContV[iV], honeycomb);
        ++nV; ++iV;
      }
      if(TiltAngle==0){
        MakeUnPairPlaneHitCluster(HC[layer2], CandContX[iX], honeycomb);
        ++nX; ++iX;
      }
      if(TiltAngle>0){
        MakeUnPairPlaneHitCluster(HC[layer2], CandContU[iU], honeycomb);
        ++nU; ++iU;
      }
    }
  }

  IndexList nCombi(npp);
  IndexList nCombiV(npp);
  IndexList nCombiX(npp);
  IndexList nCombiU(npp);

  for(Int_t i=0; i<npp; ++i){
    nCombiV[i]=(CandContV[i]).size();
    nCombiX[i]=(CandContX[i]).size();
    nCombiU[i]=(CandContU[i]).size();
  }

#if 0
  DebugPrint(nCombiV, CandContV, FUNC_NAME+" V");
  DebugPrint(nCombiX, CandContX, FUNC_NAME+" X");
  DebugPrint(nCombiU, CandContU, FUNC_NAME+" U");
#endif

  Bool_t status[3] = {true, true, true};
  std::vector<IndexList> CombiIndexV = MakeIndex(npp, nCombiV, status[0]);
  Int_t nnCombiV=CombiIndexV.size();
  std::vector<IndexList> CombiIndexX = MakeIndex(npp, nCombiX, status[1]);
  Int_t nnCombiX=CombiIndexX.size();
  std::vector<IndexList> CombiIndexU = MakeIndex(npp, nCombiU, status[2]);
  Int_t nnCombiU=CombiIndexU.size();

  for(Int_t i=0; i<nnCombiV; ++i){
    DCLocalTrack *track = MakeTrack(CandContV, CombiIndexV[i]);
    if(!track) continue;
    if(track->GetNHit()>=3 && track->DoFitVXU() &&
       track->GetChiSquare()<MaxChisquareVXU){
      TrackContV.push_back(track);
    }
    else{
      delete track;
    }
  }

  for(Int_t i=0; i<nnCombiX; ++i){
    DCLocalTrack *track = MakeTrack(CandContX, CombiIndexX[i]);
    if(!track) continue;
    if(track->GetNHit()>=3 && track->DoFitVXU() &&
       track->GetChiSquare()<MaxChisquareVXU){
      TrackContX.push_back(track);
    }
    else{
      delete track;
    }
  }

  for(Int_t i=0; i<nnCombiU; ++i){
    DCLocalTrack *track = MakeTrack(CandContU, CombiIndexU[i]);
    if(!track) continue;
    if(track->GetNHit()>=3 && track->DoFitVXU() &&
       track->GetChiSquare()<MaxChisquareVXU){
      TrackContU.push_back(track);
    }
    else{
      delete track;
    }
  }

  // Clear Flags
  if(nV>3) ClearFlags(TrackContV);
  if(nX>3) ClearFlags(TrackContX);
  if(nU>3) ClearFlags(TrackContU);

  std::stable_sort(TrackContV.begin(), TrackContV.end(), DCLTrackComp1());
  std::stable_sort(TrackContX.begin(), TrackContX.end(), DCLTrackComp1());
  std::stable_sort(TrackContU.begin(), TrackContU.end(), DCLTrackComp1());

#if 0
  DebugPrint(TrackContV, FUNC_NAME+" V After Sorting.");
  DebugPrint(TrackContX, FUNC_NAME+" X After Sorting.");
  DebugPrint(TrackContU, FUNC_NAME+" U After Sorting.");
#endif

  // Delete Duplicated Tracks (cut chisqr>100 & flag)
  Double_t chiV = ChisquareCutVXU;
  Double_t chiX = ChisquareCutVXU;
  Double_t chiU = ChisquareCutVXU;

  DeleteDuplicatedTracks(TrackContV, chiV);
  DeleteDuplicatedTracks(TrackContX, chiX);
  DeleteDuplicatedTracks(TrackContU, chiU);
  CalcTracks(TrackContV);
  CalcTracks(TrackContX);
  CalcTracks(TrackContU);

#if 0
  DebugPrint(TrackContV, FUNC_NAME+" V After Delete.");
  DebugPrint(TrackContX, FUNC_NAME+" X After Delete.");
  DebugPrint(TrackContU, FUNC_NAME+" U After Delete.");
#endif

  Int_t nnV=1, nnX=1, nnU=1;
  Int_t nkV=0, nkX=0, nkU=0;
  Int_t nnVT=1, nnXT=1, nnUT=1;
  Int_t checkV=0, checkX=0, checkU=0;

  Int_t cV=TrackContV.size();
  if(chiV>1.5 || cV<5) checkV++;
  Int_t cX=TrackContX.size();
  if(chiX>1.5 || cX<5) checkX++;
  Int_t cU=TrackContU.size();
  if(chiU>1.5 || cU<5) checkU++;

  std::vector<IndexList> CombiIndexSV;
  std::vector<IndexList> CombiIndexSX;
  std::vector<IndexList> CombiIndexSU;

  {
    if((nV>=3) && (cV)){
      nnV = TrackContV.size();
      nnVT = TrackContV.size();
      ++nkV;
    }
    if((nV>0) && checkV){
      CombiIndexSV = MakeIndex_VXU(npp, 2, nCombiV);
      nnV = nnV +  CombiIndexSV.size();
    }
  }

  {
    if((nX>=3) && (cX)){
      nnX = TrackContX.size();
      nnXT = TrackContX.size();
      ++nkX;
    }
    if((nX>0) && checkX){
      CombiIndexSX = MakeIndex_VXU(npp, 2, nCombiX);
      nnX = nnX + CombiIndexSX.size();
    }
  }

  {
    if((nU>=3) && (cU)){
      nnU = TrackContU.size();
      nnUT = TrackContU.size();
      ++nkU;
    }
    if((nU>0) && checkU){
      CombiIndexSU = MakeIndex_VXU(npp, 2, nCombiU);
      nnU = nnU + CombiIndexSU.size();
    }
  }

  // Double_t DifVXU=0.0;
  // Double_t Av=0.0;
  // Double_t Ax=0.0;
  // Double_t Au=0.0;

  Double_t chiv, chix, chiu;

#if 0
  for(Int_t i=0; i<nnV; ++i){
    for(Int_t j=0; j<nnX; ++j){
      for(Int_t k=0; k<nnU; ++k){

        chiv=-1.0,chix=-1.0,chiu=-1.0;

        DCLocalTrack *track = new DCLocalTrack();

        /* V Plane  */
        if(nkV){
          DCLocalTrack *trackV=TrackContV[i];
          //Av=trackV->GetVXU_A();
          chiv=trackV->GetChiSquare();
          for(Int_t l=0; l<(trackV->GetNHit()); ++l){
            DCLTrackHit *hitpV=trackV->GetHit(l);
            if(hitpV){
              track->AddHit(hitpV) ;
            }
          }
          //delete trackV;
        }
        if((!nkV) && (nV>0)){
          DCLocalTrack *trackV = MakeTrack(CandContV, CombiIndexV[i]);
          for(Int_t l=0; l<(trackV->GetNHit()); ++l){
            DCLTrackHit *hitpV=trackV->GetHit(l);
            if(hitpV){
              track->AddHit(hitpV) ;
            }
          }
          delete trackV;
        }

        /* X Plane  */
        if(nkX){
          DCLocalTrack *trackX=TrackContX[j];
          //Ax=trackX->GetVXU_A();
          chix=trackX->GetChiSquare();
          for(Int_t l=0; l<(trackX->GetNHit()); ++l){
            DCLTrackHit *hitpX=trackX->GetHit(l);
            if(hitpX){
              track->AddHit(hitpX) ;
            }
          }
          //delete trackX;
        }
        if((!nkX) && (nX>0)){
          DCLocalTrack *trackX = MakeTrack(CandContX, CombiIndexX[j]);
          for(Int_t l=0; l<(trackX->GetNHit()); ++l){
            DCLTrackHit *hitpX=trackX->GetHit(l);
            if(hitpX){
              track->AddHit(hitpX) ;
            }
          }
          delete trackX;
        }

        /* U Plane  */
        if(nkU){
          DCLocalTrack *trackU=TrackContU[k];
          //Au=trackU->GetVXU_A();
          chiu=trackU->GetChiSquare();
          for(Int_t l=0; l<(trackU->GetNHit()); ++l){
            DCLTrackHit *hitpU=trackU->GetHit(l);
            if(hitpU){
              track->AddHit(hitpU) ;
            }
          }
          //delete trackU;
        }
        if((!nkU) && (nU>0)){
          DCLocalTrack *trackU = MakeTrack(CandContU, CombiIndexU[k]);
          for(Int_t l=0; l<(trackU->GetNHit()); ++l){
            DCLTrackHit *hitpU=trackU->GetHit(l);
            if(hitpU){
              track->AddHit(hitpU) ;
            }
          }
          delete trackU;
        }
        //track->SetAv(Av);
        //track->SetAx(Ax);
        //track->SetAu(Au);
        //DifVXU = track->GetDifVXU();

        track->SetChiv(chiv);
        track->SetChix(chix);
        track->SetChiu(chiu);

        if(!track) continue;
        if(track->GetNHit()>=MinNumOfHits && track->DoFit() &&
           track->GetChiSquare()<MaxChisquare){
          TrackCont.push_back(track);
        }
        else{
          delete track;
        }
      }
    }
  }
#endif

  for(Int_t i=-1; i<nnV; ++i){
    for(Int_t j=-1; j<nnX; ++j){
      for(Int_t k=-1; k<nnU; ++k){
        if(((i+j)==-2) || ((j+k)==-2) || ((k+i)==-2)) continue;

        chiv=-1.0,chix=-1.0,chiu=-1.0;

        DCLocalTrack *track = new DCLocalTrack();

        /* V Plane  */
        if(i>-1){
          if(nkV && i<nnVT){
            DCLocalTrack *trackV=TrackContV[i];
            // Av=trackV->GetVXU_A();
            chiv=trackV->GetChiSquare();

            for(Int_t l=0; l<(trackV->GetNHit()); ++l){
              DCLTrackHit *hitpV=trackV->GetHit(l);
              if(hitpV){
                track->AddHit(hitpV) ;
              }
            }
          }
          if((i>=nnVT) && (nV>0)){
            DCLocalTrack *trackV = MakeTrack(CandContV, CombiIndexSV[i-nnVT]);
            for(Int_t l=0; l<(trackV->GetNHit()); ++l){
              DCLTrackHit *hitpV=trackV->GetHit(l);
              if(hitpV){
                track->AddHit(hitpV) ;
              }
            }
            delete trackV ;
          }
          // Int_t NHitV = track->GetNHit();
        }

        /* X Plane  */
        if(j>-1){
          if(nkX && j<nnXT){
            DCLocalTrack *trackX=TrackContX[j];
            // Ax=trackX->GetVXU_A();
            chix=trackX->GetChiSquare();
            for(Int_t l=0; l<(trackX->GetNHit()); ++l){
              DCLTrackHit *hitpX=trackX->GetHit(l);
              if(hitpX){
                track->AddHit(hitpX) ;
              }
            }
          }
          if((j>=nnXT) && (nX>0)){
            DCLocalTrack *trackX = MakeTrack(CandContX, CombiIndexSX[j-nnXT]);
            for(Int_t l=0; l<(trackX->GetNHit()); ++l){
              DCLTrackHit *hitpX=trackX->GetHit(l);
              if(hitpX){
                track->AddHit(hitpX) ;
              }
            }
            delete trackX;
          }
          // Int_t NHitX = track->GetNHit();
        }
        /* U Plane  */
        if(k>-1){
          if(nkU && k<nnUT){
            DCLocalTrack *trackU=TrackContU[k];
            // Au=trackU->GetVXU_A();
            chiu=trackU->GetChiSquare();
            for(Int_t l=0; l<(trackU->GetNHit()); ++l){
              DCLTrackHit *hitpU=trackU->GetHit(l);
              if(hitpU){
                track->AddHit(hitpU) ;
              }
            }
          }
          if((k>=nnUT) && (nU>0)){
            DCLocalTrack *trackU = MakeTrack(CandContU, CombiIndexSU[k-nnUT]);
            for(Int_t l=0; l<(trackU->GetNHit()); ++l){
              DCLTrackHit *hitpU=trackU->GetHit(l);
              if(hitpU){
                track->AddHit(hitpU) ;
              }
            }
            delete trackU;
          }
          // Int_t NHitU = track->GetNHit();
        }

        //track->SetAv(Av);
        //track->SetAx(Ax);
        //track->SetAu(Au);
        //DifVXU = track->GetDifVXU();

        track->SetChiv(chiv);
        track->SetChix(chix);
        track->SetChiu(chiu);

        if(!track) continue;
        if(track->GetNHit()>=MinNumOfHits && track->DoFit() &&
           track->GetChiSquare()<MaxChisquare){//MaXChisquare
          TrackCont.push_back(track);
        }
        else{
          delete track;
        }
      }
    }
  }

  {

    for(Int_t i=0; i<Int_t(TrackContV.size()); ++i){
      DCLocalTrack *tp=TrackContV[i];
      delete tp;
      TrackContV.erase(TrackContV.begin()+i);
    }
    for(Int_t i=0; i<Int_t(TrackContX.size()); ++i){
      DCLocalTrack *tp=TrackContX[i];
      delete tp;
      TrackContX.erase(TrackContX.begin()+i);
    }
    for(Int_t i=0; i<Int_t(TrackContU.size()); ++i){
      DCLocalTrack *tp=TrackContU[i];
      delete tp;
      TrackContU.erase(TrackContU.begin()+i);
    }

  }

  ClearFlags(TrackCont);

  std::stable_sort(TrackCont.begin(), TrackCont.end(), DCLTrackComp1());

#if 1
  // Delete Tracks about  (Nhit1 > Nhit2+1) (Nhit1 > Nhit2  && chi1 < chi2)
  for(Int_t i=0; i<Int_t(TrackCont.size()); ++i){
    DCLocalTrack *tp=TrackCont[i];
    Int_t nh=tp->GetNHit();
    Double_t chi=tp->GetChiSquare();
    for(Int_t j=0; j<nh; ++j) tp->GetHit(j)->JoinTrack();

    for(Int_t i2=TrackCont.size()-1; i2>i; --i2){
      DCLocalTrack *tp2=TrackCont[i2];
      Int_t nh2=tp2->GetNHit(), flag=0;
      Double_t chi2=tp2->GetChiSquare();
      for(Int_t j=0; j<nh2; ++j)
        if(tp2->GetHit(j)->BelongToTrack()) ++flag;
      if((flag>=2) && ((nh==nh2) || ((nh>nh2) && (chi<chi2)))){
        //      if((flag) && ((nh>nh2+1) || ((nh==nh2) || (nh>nh2) && (chi<chi2)))){
        //if((nh>nh2) && (chi<chi2)){
        delete tp2;
        TrackCont.erase(TrackCont.begin()+i2);
      }
    }
  }
#endif

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackCompSdcInFiber(), CandCont);

  std::stable_sort(TrackCont.begin(), TrackCont.end(), DCLTrackComp());

  // Clear Flags
  {
    Int_t nbefore=TrackCont.size();
    for(Int_t i=0; i<nbefore; ++i){
      DCLocalTrack *tp=TrackCont[i];
      Int_t nh=tp->GetNHit();
      for(Int_t j=0; j<nh; ++j){
        tp->GetHit(j)->QuitTrack();
      }
    }
  }

  // Delete Duplicated Tracks

  for(Int_t i=0; i<Int_t(TrackCont.size()); ++i){
    DCLocalTrack *tp=TrackCont[i];
    Int_t nh=tp->GetNHit();
    for(Int_t j=0; j<nh; ++j) tp->GetHit(j)->JoinTrack();

    for(Int_t i2=TrackCont.size()-1; i2>i; --i2){
      DCLocalTrack *tp2=TrackCont[i2];
      Int_t nh2=tp2->GetNHit(), flag=0;
      for(Int_t j=0; j<nh2; ++j)
        if(tp2->GetHit(j)->BelongToTrack()) ++flag;
      if(flag){
        delete tp2;
        TrackCont.erase(TrackCont.begin()+i2);
      }
    }
  }

  // Clear Flags
  {
    Int_t nbefore=TrackCont.size();
    for(Int_t i=0; i<nbefore; ++i){
      DCLocalTrack *tp=TrackCont[i];
      Int_t nh=tp->GetNHit();
      for(Int_t j=0; j<nh; ++j){
        tp->GetHit(j)->QuitTrack();
      }
    }
  }

  {
    Int_t nn=TrackCont.size();
    for(Int_t i=0; i<nn; ++i){
      DCLocalTrack *tp=TrackCont[i];
      Int_t nh=tp->GetNHit();
      for(Int_t j=0; j<nh; ++j){
        Int_t lnum = tp->GetHit(j)->GetLayer();
        Double_t zz = gGeom.GetLocalZ(lnum);
        tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }

#if 0
  if(TrackCont.size()>0){
    Int_t nn=TrackCont.size();
    hddaq::cout << FUNC_NAME << ": After Deleting. #Tracks = " << nn << std::endl;
    for(Int_t i=0; i<nn; ++i){
      DCLocalTrack *track=TrackCont[i];
      hddaq::cout << std::setw(3) << i << " #Hits="
                  << std::setw(2) << track->GetNHit()
                  << " ChiSqr=" << track->GetChiSquare()
                  << std::endl;
      hddaq::cout << std::endl;
      for(Int_t j=0; j<(track->GetNHit()); ++j){
        DCLTrackHit *hit = track->GetHit(j);
        hddaq::cout << "layer = " << hit->GetLayer()+1 << " Res = " << hit->GetResidual() << std::endl ;
        hddaq::cout << "hitp = " << hit->GetLocalHitPos() << " calp = " << hit->GetLocalCalPos() << std::endl ;
        hddaq::cout << "X = " << hit->GetXcal() << " Y = " << hit->GetYcal() << std::endl ;
        hddaq::cout << std::endl;
      }
      hddaq::cout << "*********************************************" << std::endl;
    }
    hddaq::cout << std::endl;
  }
#endif

  del::ClearContainerAll(CandCont);
  del::ClearContainerAll(CandContV);
  del::ClearContainerAll(CandContX);
  del::ClearContainerAll(CandContU);

  Bool_t status_all = true;
  status_all = status_all && status[0];
  status_all = status_all && status[1];
  status_all = status_all && status[2];

  return status_all? TrackCont.size() : -1;
}

//_____________________________________________________________________________
Int_t
LocalTrackSearchBcOutSdcIn(const std::vector<DCHC>& BcHC,
                           const DCPairPlaneInfo *BcPpInfo,
                           const std::vector<DCHC>& SdcHC,
                           const DCPairPlaneInfo *SdcPpInfo,
                           Int_t BcNpp, Int_t SdcNpp,
                           DCLocalTC& TrackCont,
                           Int_t MinNumOfHits)
{
  const Int_t npp = BcNpp + SdcNpp;

  std::vector<ClusterList> CandCont(npp);

  for(Int_t i=0; i<BcNpp; ++i){
    Bool_t ppFlag=BcPpInfo[i].pair;
    Int_t layer1=BcPpInfo[i].id1, layer2=BcPpInfo[i].id2;
    if(ppFlag)
      MakePairPlaneHitCluster(BcHC[layer1], BcHC[layer2],
                              BcPpInfo[i].CellSize, CandCont[i]);
    else
      MakeUnPairPlaneHitCluster(BcHC[layer1], CandCont[i]);
  }

  for(Int_t i=0; i<SdcNpp; ++i){
    Bool_t ppFlag=SdcPpInfo[i].pair;
    Int_t layer1=SdcPpInfo[i].id1, layer2=SdcPpInfo[i].id2;
    if(ppFlag)
      MakePairPlaneHitCluster(SdcHC[layer1], SdcHC[layer2],
                              SdcPpInfo[i].CellSize, CandCont[i+BcNpp]);
    else
      MakeUnPairPlaneHitCluster(SdcHC[layer1], CandCont[i+BcNpp]);
  }

  IndexList nCombi(npp);
  for(Int_t i=0; i<npp; ++i){
    Int_t n = CandCont[i].size();
    nCombi[i] = n>MaxNumOfCluster ? 0 : n;
  }

#if 0
  DebugPrint(nCombi, CandCont, FUNC_NAME);
#endif

  Bool_t status = true;
  std::vector<IndexList> CombiIndex = MakeIndex(npp, nCombi, status);
  Int_t nnCombi = CombiIndex.size();

#if 0
  hddaq::cout << " ===> " << nnCombi << " combinations will be checked.."
              << std::endl;
  for(Int_t i=0; i<nnCombi; ++i){
    for(Int_t j=0;j<npp;j++) {
      hddaq::cout << CombiIndex[i][j] << " ";
    }
    hddaq::cout << std::endl;
  }
#endif

  for(Int_t i=0; i<nnCombi; ++i){
    DCLocalTrack *track = MakeTrack(CandCont, CombiIndex[i]);
    if(!track) continue;
    if(track->GetNHit()>=MinNumOfHits && track->DoFitBcSdc() &&
       track->GetChiSquare()<MaxChisquare)
      TrackCont.push_back(track);
    else
      delete track;
  }

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackComp(), CandCont);

  return TrackCont.size();
}

//_____________________________________________________________________________
Int_t
LocalTrackSearchCFT( const std::vector<CFTFiberCC>& HC,
		     Int_t NPlane,
		     CFTLocalTC& TrackCont,
		     Int_t MinNumOfHitsPhi,  Int_t MinNumOfHitsU)
{
  static const std::string funcname = "[LocalTrackSearchCFT]";

  Double_t cut_range = 250.; // center-250~center+250

  IndexList nCombi(NPlane);

  for( Int_t i=0; i<NPlane; ++i ){
    nCombi[i]=(HC[i]).size();

    if( nCombi[i]>MaxNumOfClusterCFT ) nCombi[i]=0;
  }

  Bool_t status = true;
  std::vector < IndexList >
    CombiIndex = MakeIndex( NPlane, nCombi, status );
  Int_t nnCombi=CombiIndex.size();

#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.."
	    << std::endl;
#endif
  if( nnCombi>MaxCombi ) return 0;

  for( Int_t i=0; i<nnCombi; ++i ){
    CFTLocalTrack *track = MakeTrack( HC, CombiIndex[i]);
    Int_t nMinus = 0;
    for (Int_t j=0; j<CombiIndex[i].size(); j++) {
      if ( CombiIndex[i][j] == -1)
	nMinus++;
    }
    /*
    if (nMinus==0) {
      std::cout << "( " ;
      for (int j=0; j<CombiIndex[i].size(); j++) {
	std::cout << CombiIndex[i][j] << ", " ;
      }
      std::cout << " ) " << std::endl;
    }
    */
    if( !track ) continue;

    if( track->GetNHit()>=MinNumOfHitsPhi&& track->CheckPhi_1st() && track->DoFitXY() &&
	track->GetChiSquareXY()<MaxChisquareCFT) {
#if 0
      std::cout << "NHit = " << track->GetNHit() << " Chisqr = " << track->GetChiSquareXY() << std::endl;

      std::cout << "******** i = " << i << " *******************" << std::endl;
#endif

      if (track->GetNHitU()>=MinNumOfHitsU &&
	  track->DoFitZTrack() &&
	  track->GetChiSquareZ()<MaxChisquareCFT  &&
	  fabs(track->GetVtxZ())<cut_range // vtx cut
	  ){

	//std::cout << "pass 1st Z Tracking : chi2 = " << track->GetChiSquareZ()
	//<< std::endl;

	track->SetCalculatedValue();

	track->SetCorPhi();

	// 2nd tracking
	track->DoFitXY_2nd();
	track->DoFitZTrack();
	Double_t chi2_0 = track->GetChiSquareZ();

	track->SetThetaCFT();

	track->DoFitZTrack();
	Double_t chi2_1 = track->GetChiSquareZ();
	if (chi2_1 > chi2_0 && chi2_1>1.5 && track->GetNHitU()==4) {
	  track->ResetThetaCFT();
	  track->DoFitZTrack();
	}

	Bool_t flag_chi2_uv = true;
	if (track->GetNHitU()==3 && track->GetChiSquareZ() > 30)
	  flag_chi2_uv = false;

	Bool_t flag_chi2_phi = true;
	if (track->GetNHit()==3 && track->GetChiSquareXY() > 30)
	  flag_chi2_phi = false;

	/*
	std::cout << "2nd tracking : chi2(XY) = " << track->GetChiSquareXY()
		  << ", nhitX = " << track->GetNHit()
		  << ", chi2(Z) = " << track->GetChiSquareZ()
		  << ", nhitZ = " << track->GetNHitU()
		  << std::endl;
	*/
	if (track->GetChiSquareXY() < MaxChi2CFT2nd_phi  &&
	    track->GetChiSquareZ() < MaxChi2CFT2nd_uv  &&
	    flag_chi2_uv &&
	    flag_chi2_phi &&
	    fabs(track->GetVtxZ())<cut_range /*&&
	    track->SetCalculatedValue() &&
	    track->GetChiSquareDE() < 10*/) {
	  TrackCont.push_back(track);
	} else {
	  delete track;
	}
      } else {
	delete track;
      }
    }
    else {
      delete track;
    }

  }

  // Clear Flags
  Int_t nbefore=TrackCont.size();
  for( Int_t i=0; i<nbefore; ++i ){
    CFTLocalTrack *tp=TrackCont[i];
    Int_t nh=tp->GetNHit();
    for( Int_t j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();

    Int_t nhU=tp->GetNHitU();
    for( Int_t j=0; j<nhU; ++j ) tp->GetHitU(j)->clearFlags();
  }

#if 0
  {
    Int_t nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = "
 	      << nn << std::endl;
    for( Int_t i=0; i<nn; ++i ){
      CFTLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits (XY)="
 		<< std::setw(2) << track->GetNHit()
 		<< " ChiSqr(XY)=" << track->GetChiSquareXY()
		<< std::setw(3) << i << " #Hits (Z)="
 		<< std::setw(2) << track->GetNHitU()
 		<< " ChiSqr(Z)=" << track->GetChiSquareXY()
 		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  //partial_sort( TrackCont.begin(), TrackCont.end(),
  //TrackCont.end(), CFTLTrackComp() );
  std::stable_sort( TrackCont.begin(), TrackCont.end(), CFTLTrackComp() );

#if 0
  {
    Int_t nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = "
 	      << nn << std::endl;

    for( Int_t i=0; i<nn; ++i ){
      CFTLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits(XY)="
		<< std::setw(2) << track->GetNHit()
		<< " ChiSqr(XY)=" << track->GetChiSquareXY()
		<< std::setw(3) << i << " #Hits(Z)="
		<< std::setw(2) << track->GetNHitU()
		<< " ChiSqr(Z)=" << track->GetChiSquareZ()
		<< " Total ChiSqr=" << track->GetChiSquareZ() + track->GetChiSquareXY()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif


  // Delete Duplicated Tracks
  for( Int_t i=0; i<Int_t(TrackCont.size()); ++i ){
    CFTLocalTrack *tp=TrackCont[i];

    Int_t nh=tp->GetNHit();
    for( Int_t j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

    Int_t nhU=tp->GetNHitU();
    for( Int_t j=0; j<nhU; ++j ) tp->GetHitU(j)->setFlags();

     for( Int_t i2=TrackCont.size()-1; i2>i; --i2 ){
       CFTLocalTrack *tp2=TrackCont[i2];
       Int_t nh2=tp2->GetNHit();
       Int_t flag=0;
       for( Int_t j=0; j<nh2; ++j ) {
	 //std::cout << "Phi " << j << " : " <<  tp2->GetHit(j)->showFlags() << std::endl;
	 if( tp2->GetHit(j)->showFlags() ) ++flag;
       }

       Int_t nhU2=tp2->GetNHitU();
       for( Int_t j=0; j<nhU2; ++j ) {
	 //std::cout << "U " << j << " : " <<  tp2->GetHitU(j)->showFlags() << std::endl;
	 if( tp2->GetHitU(j)->showFlags() ) ++flag;
       }
       if(flag){
	 //std::cout << "delete" << std::endl;
	 delete tp2;
	 TrackCont.erase(TrackCont.begin()+i2);
       } else {
	 //std::cout << "accept" << std::endl;
       }
     }
   }

  {
    Int_t nn=TrackCont.size();
    for(Int_t i=0; i<nn; ++i ){
      CFTLocalTrack *tp=TrackCont[i];
      tp->SetCalculatedValue();
    }
  }

#if 0
  //  if(TrackCont.size()>1){
  Int_t nn=TrackCont.size();
  std::cout << funcname << ": After Deleting. #Tracks = "
	    << nn << std::endl;

  for( Int_t i=0; i<nn; ++i ){
    CFTLocalTrack *track=TrackCont[i];
    std::cout << std::setw(3) << i << " #Hits(XY)="
	      << std::setw(2) << track->GetNHit()
	      << " ChiSqr(XY)=" << track->GetChiSquareXY()
	      << std::setw(3) << i << " #Hits(Z)="
	      << std::setw(2) << track->GetNHitU()
	      << " ChiSqr(Z)=" << track->GetChiSquareZ()
	      << " Total ChiSqr=" << track->GetChiSquareZ() + track->GetChiSquareXY()
	      << " Theta=" << track->GetThetaCFT()
	      << std::endl;
    std::cout << std::endl;

    /*
      for( int j=0; j<(track->GetNHit()); ++j){
      DCLTrackHit *hit = track->GetHit(j);
      std::cout << "layer = " << hit->GetLayer()+1 << " Res = " << hit->GetResidual() << std::endl ;
      std::cout << "hitp = " << hit->GetLocalHitPos() << " calp = " << hit->GetLocalCalPos() << std::endl ;
      std::cout << "X = " << hit->GetXcal() << " Y = " << hit->GetYcal() << std::endl ;
      std::cout << std::endl;
       }
    */
    std::cout << "*********************************************" << std::endl;
  }
  std::cout << std::endl;

  //  }
#endif

  Int_t ntCFT = TrackCont.size();
  for (Int_t i=0; i<ntCFT; i++) {
    CFTLocalTrack *track=TrackCont[i];
    track->CalcNormalizedDE();
  }

  return TrackCont.size();
}

//_____________________________________________________________________________
Int_t LocalTrackSearchCFT_16layer( CFTLocalTC& TrackCont,
				   CFTLocalTC& TrackCont_16layer)
{
  static const std::string funcname = "[LocalTrackSearchCFT_16layer]";

  CFTLocalTrack *track = new CFTLocalTrack();

  if(!track){
    std::cerr << funcname << ": new fail" << std::endl;
    return 0;
  }


  if (TrackCont.size() != 2)
    return 0;

  CFTLocalTrack *track1 = TrackCont[0];
  Int_t nhitPhi1 = track1->GetNHit();
  for (Int_t i=0; i<nhitPhi1; i++) {
    CFTFiberCluster *fcl = track1->GetHit(i);
    track->AddHit(fcl);
  }
  Int_t nhitUV1 = track1->GetNHitU();
  for (Int_t i=0; i<nhitUV1; i++) {
    CFTFiberCluster *fcl = track1->GetHitU(i);
    track->AddHitU(fcl);
  }

  CFTLocalTrack *track2 = TrackCont[1];
  Int_t nhitPhi2 = track2->GetNHit();
  for (Int_t i=0; i<nhitPhi2; i++) {
    CFTFiberCluster *fcl = track2->GetHit(i);
    track->AddHit(fcl);
  }
  Int_t nhitUV2 = track2->GetNHitU();
  for (Int_t i=0; i<nhitUV2; i++) {
    CFTFiberCluster *fcl = track2->GetHitU(i);
    track->AddHitU(fcl);
  }

  if (track->DoFitXY() && track->GetChiSquareXY() < MaxChisquareCFT) {
    //std::cout << "8layer : chi2_XY = "  << track->GetChiSquareXY() << std::endl;
    if (track->DoFitZTrack_16layer() && track->GetChiSquareZ()<MaxChisquareCFT ){
      //std::cout << "8layer : chi2_Z = "  << track->GetChiSquareZ() << std::endl;
      TrackCont_16layer.push_back(track);
    } else {
      //std::cout << "8layer : chi2_Z = "  << track->GetChiSquareZ() << std::endl;
      delete track;
    }
  } else {
    //std::cout << "8layer : chi2_XY = "  << track->GetChiSquareXY() << std::endl;
    delete track;
  }


  {
    Int_t nn=TrackCont_16layer.size();
    for(Int_t i=0; i<nn; ++i ){
      CFTLocalTrack *tp=TrackCont_16layer[i];
      tp->SetCalculatedValue_16layer();

      tp->SetCorPhi();

      // 2nd tracking
      tp->DoFitXY_2nd();
      tp->DoFitZTrack_16layer();
      tp->SetCalculatedValue_16layer();

    }
  }

#if 0
  //  if(TrackCont.size()>1){
  Int_t nn=TrackCont_16layer.size();
  std::cout << funcname << ": After Deleting. #Tracks = "
	       << nn << std::endl;

  for( Int_t i=0; i<nn; ++i ){
    CFTLocalTrack *track=TrackCont_8layer[i];
    std::cout << std::setw(3) << i << " #Hits(XY)="
	      << std::setw(2) << track->GetNHit()
	      << " ChiSqr(XY)=" << track->GetChiSquareXY()
	      << std::setw(3) << i << " #Hits(Z)="
	      << std::setw(2) << track->GetNHitU()
	      << " ChiSqr(Z)=" << track->GetChiSquareZ()
	      << " Total ChiSqr=" << track->GetChiSquare()
	      << std::endl;
    std::cout << std::endl;

    /*
      for( int j=0; j<(track->GetNHit()); ++j){
      DCLTrackHit *hit = track->GetHit(j);
      std::cout << "layer = " << hit->GetLayer()+1 << " Res = " << hit->GetResidual() << std::endl ;
      std::cout << "hitp = " << hit->GetLocalHitPos() << " calp = " << hit->GetLocalCalPos() << std::endl ;
      std::cout << "X = " << hit->GetXcal() << " Y = " << hit->GetYcal() << std::endl ;
      std::cout << std::endl;
       }
    */
    std::cout << "*********************************************" << std::endl;
  }
  std::cout << std::endl;

  //  }
#endif

  Int_t ntCFT = TrackCont_16layer.size();
  for (Int_t i=0; i<ntCFT; i++) {
    CFTLocalTrack *track=TrackCont_16layer[i];
    track->CalcNormalizedDE();
  }

  return TrackCont_16layer.size();
}


//_____________________________________________________________________________
inline Bool_t
IsDeletionTarget(const std::vector<std::pair<Int_t,Int_t>>& nh,
                 std::size_t NDelete, Int_t layer)
{
  if(NDelete==0) return false;

  for(auto itr=nh.begin(), end=nh.end(); itr!=end; ++itr){
    if(itr->first==0 && itr->second==layer)
      return true;
  }

  if(nh[0].second == layer)
    return true;

  return false;
}

//_____________________________________________________________________________
Int_t
LocalTrackSearchSdcInSdcOut(const std::vector<DCHC>& SdcInHC,
                            const DCPairPlaneInfo* SdcInPpInfo,
                            const std::vector<DCHC>& SdcOutHC,
                            const DCPairPlaneInfo* SdcOutPpInfo,
                            Int_t SdcInNpp, Int_t SdcOutNpp,
                            DCLocalTC& TrackCont,
                            Int_t MinNumOfHits)
{
  const Int_t npp = SdcInNpp+SdcOutNpp;

  std::vector<ClusterList> CandCont(npp);

  for(Int_t i=0; i<SdcInNpp; ++i){
    Bool_t ppFlag = SdcInPpInfo[i].pair;
    Int_t  layer1 = SdcInPpInfo[i].id1;
    Int_t  layer2 = SdcInPpInfo[i].id2;
    if(ppFlag){
      MakePairPlaneHitCluster(SdcInHC[layer1], SdcInHC[layer2],
                              SdcInPpInfo[i].CellSize, CandCont[i]);
    } else {
      MakeUnPairPlaneHitCluster(SdcInHC[layer1], CandCont[i]);
    }
  }

  for(Int_t i=0; i<SdcOutNpp; ++i){
    Bool_t ppFlag = SdcOutPpInfo[i].pair;
    Int_t  layer1 = SdcOutPpInfo[i].id1;
    Int_t  layer2 = SdcOutPpInfo[i].id2;
    if(ppFlag){
      MakePairPlaneHitCluster(SdcOutHC[layer1], SdcOutHC[layer2],
                              SdcOutPpInfo[i].CellSize, CandCont[i+SdcInNpp]);
    } else {
      MakeUnPairPlaneHitCluster(SdcOutHC[layer1], CandCont[i+SdcInNpp]);
    }
  }

  IndexList nCombi(npp);
  for(Int_t i=0; i<npp; ++i){
    Int_t n = CandCont[i].size();
    nCombi[i] = n>MaxNumOfCluster ? 0 : n;
  }

  Bool_t status = true;
  std::vector<IndexList> CombiIndex = MakeIndex(npp, nCombi, status);
  Int_t nnCombi = CombiIndex.size();

  for(Int_t i=0; i<nnCombi; ++i){
    DCLocalTrack *track = MakeTrack(CandCont, CombiIndex[i]);
    if(!track) continue;
    if(track->GetNHit()>=MinNumOfHits &&
       track->DoFit() &&
       track->GetChiSquare()<MaxChisquare)
      TrackCont.push_back(track);
    else
      delete track;
  }

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackComp(), CandCont);
  return status? TrackCont.size() : -1;
}

//For MWPC
//_____________________________________________________________________________
Int_t MWPCLocalTrackSearch(const std::vector<DCHC>& HC,
                           DCLocalTC& TrackCont)

{
  std::vector<ClusterList> CandCont(NumOfLayersBcIn);

  for(Int_t i=0; i<NumOfLayersBcIn; ++i){
    MakeMWPCPairPlaneHitCluster(HC[i], CandCont[i]);
  }

  //   IndexList nCombi(NumOfLayersBcIn);
  //   for(Int_t i=0; i<NumOfLayersBcIn; ++i){
  //     nCombi[i]=(CandCont[i]).size();

  //     // If #Cluster>MaxNumerOfCluster,  error return

  //     if(nCombi[i]>MaxNumOfCluster){
  //       hddaq::cout << FUNC_NAME << " too many clusters " << FUNC_NAME
  //                   << "  layer = " << i << " : " << nCombi[i] << std::endl;
  //       del::ClearContainerAll(CandCont);
  //       return 0;
  //     }
  //   }

  //   DebugPrint(nCombi, FUNC_NAME);

#if 0
  DebugPrint(nCombi, CandCont, FUNC_NAME);
#endif

  const Int_t MinNumOfHitsBcIn   = 6;
  TrackMaker trackMaker(CandCont, MinNumOfHitsBcIn, MaxCombi, MaxChisquare);
  trackMaker.MakeTracks(TrackCont);

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackComp(), CandCont);
  return TrackCont.size();
}

//_____________________________________________________________________________
Int_t
MWPCLocalTrackSearch(const std::vector<std::vector<DCHC>>& hcList,
                     DCLocalTC& trackCont)
{
  for(auto itr=hcList.begin(), end=hcList.end(); itr!=end; ++itr){
    const std::vector<DCHC>& l = *itr;
    DCLocalTC tc;
    MWPCLocalTrackSearch(l, tc);
    trackCont.insert(trackCont.end(), tc.begin(), tc.end());
    // hddaq::cout << " tc " << tc.size()
    //             << " : " << trackCont.size() << std::endl;
  }

  ClearFlags(trackCont);
  DeleteDuplicatedTracks(trackCont);
  CalcTracks(trackCont);
  return trackCont.size();
}
}
