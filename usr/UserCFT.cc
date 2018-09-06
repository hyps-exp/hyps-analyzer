/**
 *  file: UserCFT.cc
 *  date: 2018.07.20
 *  based on UserBFT.cc
 */

#include <cmath>
#include <iostream>
#include <sstream>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"

#include "Kinematics.hh"
#include "EventDisplayCFT.hh"
#include "CFTParticle.hh"

#define HodoCut 0 // with BH1/BH2
#define TimeCut 0 // in cluster analysis

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

namespace
{
  using namespace root;
  const std::string& class_name("EventCFT");
  RMAnalyzer&         gRM   = RMAnalyzer::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const bool& FlagEvDisp = ConfMan::Get<bool>("EVDISP_CFT");  
}

//______________________________________________________________________________
VEvent::VEvent( void )
{
}

//______________________________________________________________________________
VEvent::~VEvent( void )
{
}

//______________________________________________________________________________
class EventCFT : public VEvent
{
private:
  RawData      *rawData;
  DCAnalyzer   *DCAna;
  HodoAnalyzer *hodoAna;

public:
        EventCFT( void );
       ~EventCFT( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventCFT::EventCFT( void )
  : VEvent(),
    rawData(0),
    DCAna( new DCAnalyzer ),
    hodoAna( new HodoAnalyzer )
{
}

//______________________________________________________________________________
EventCFT::~EventCFT( void )
{
  if ( hodoAna ){
    delete hodoAna;
    hodoAna = NULL;
  }
  if ( DCAna ){
    delete DCAna;
    DCAna   = NULL;
  }
  if ( rawData ){
    delete rawData;
    rawData = NULL;
  }
}

//______________________________________________________________________________
struct Event
{
  int evnum;

  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  // Fiber Hit
  
  int    nhits;
  int    unhits;
  int    dnhits;

  // CFT
  int ntCFT;
  int ncl;
  double phi[MaxDepth];
  double theta[MaxDepth];

  int seg[NumOfPlaneCFT][MaxDepth];

  double dphi[NumOfPlaneCFT][MaxDepth];
  double phi_ini[NumOfPlaneCFT][MaxDepth];
  double phi_track[NumOfPlaneCFT][MaxDepth];

  double dz[NumOfPlaneCFT][MaxDepth];
  double z_ini[NumOfPlaneCFT][MaxDepth];
  double z_track[NumOfPlaneCFT][MaxDepth];
  ThreeVector Pos[MaxDepth];
  ThreeVector Dir[MaxDepth];
  double vtx_x[MaxDepth], vtx_y[MaxDepth], vtx_z[MaxDepth];
  double vtxAB_x, vtxAB_y, vtxAB_z;
  
  //double pede[NumOfPlaneCFT][60];         // N_group_cable
  //double adc_raw_h[NumOfPlaneCFT][60][16];// N_group_cable
  //double adc_raw_l[NumOfPlaneCFT][60][16];// N_group_cable
  double adc_raw_h[NumOfPlaneCFT][NumOfSegCFT_PHI4];
  double adc_raw_l[NumOfPlaneCFT][NumOfSegCFT_PHI4];

  double adc[NumOfPlaneCFT][MaxDepth], adc_max[NumOfPlaneCFT][MaxDepth];
  double mip[NumOfPlaneCFT][MaxDepth], mip_max[NumOfPlaneCFT][MaxDepth];
  double dE[NumOfPlaneCFT][MaxDepth] , dE_max[NumOfPlaneCFT][MaxDepth] ;
  double TotaldE[MaxDepth],    TotaldE_max[MaxDepth];
  double TotaldEphi[MaxDepth], TotaldEphi_max[MaxDepth];
  double TotaldEuv[MaxDepth],  TotaldEuv_max[MaxDepth];

  // for cosmic ray tracking
  int seg16[NumOfPlaneCFT][2];

  double dphi16[NumOfPlaneCFT][2];
  double phi16_ini[NumOfPlaneCFT][2];
  double phi16_track[NumOfPlaneCFT][2];
  double phi16__track[NumOfPlaneCFT][2];

  double dz16[NumOfPlaneCFT][2];
  double z16_ini[NumOfPlaneCFT][2];
  double z16_track[NumOfPlaneCFT][2];

  // BGO
  int nhBGO;
  int segBGO[NumOfSegBGO];
  int segBGOt[NumOfSegBGO];// matched to track
  double adcbgo[NumOfSegBGO];
  double tdcbgo[NumOfSegBGO];

};

//______________________________________________________________________________
namespace root
{
  Event event;
  TH1   *h[MaxHist];
  TTree *tree;
}

//______________________________________________________________________________
bool
EventCFT::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
bool
EventCFT::ProcessingNormal( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"]");

  bool disp_flag = false;

  rawData = new RawData;
  rawData->DecodeHits();

  gRM.Decode();

  int evnum = gRM.EventNumber();
  event.evnum = evnum;

  //**************************************************************************
  //******************RawData

  // BGO
  {
    const HodoRHitContainer &cont = rawData->GetBGORawHC();
    int nhit = cont.size();
    for(int i = 0; i<nhit; ++i){	
      HodoRawHit *hit = cont.at(i);
      int seg = hit->SegmentId(); 
      int NhitT = hit->GetSizeTdcUp();
      int NhitA = hit->GetSizeAdcUp();	  

      for(int m = 0; m<NhitT; ++m){
	int tdc = hit->GetTdcUp(m);	
	HF2 (100, seg, tdc);
	//std::cout << "BGO raw tdc : seg=" << seg << ", tdc=" << tdc << ", m=" << m << std::endl;	
      }
      for(int m = 0; m<NhitA; ++m){
	int integral = hit->GetAdcUp();	
	HF2 (101, seg, integral);
	//if(NhitT>0){HF2 (111, seg, integral);}
	//std::cout << "BGO raw integral : seg=" << seg << ", integral=" << integral << std::endl;	
      }      
     
    }
  }
  hodoAna->DecodeBGOHits(rawData);
  int nhBGO = hodoAna->GetNHitsBGO();
  event.nhBGO = nhBGO;
  for(int i=0; i<nhBGO; ++i){
    Hodo1Hit *hit = hodoAna->GetHitBGO(i);
    int seg = hit->SegmentId();
    int nh = hit->GetNumOfHit();
    double adc = hit->GetAUp();

    event.segBGO[i] = seg;
    event.adcbgo[seg] = adc;

    for(int m=0; m<nh; ++m){
      double cmt = hit->CMeanTime(m);
      event.tdcbgo[seg] = cmt;
      HF2 (110, seg, cmt);
      if(-50<cmt&&cmt<50){
	HF2 (111, seg, adc);
	
	if (FlagEvDisp && adc>0) {
	  std::cout << "BGO decode : seg=" << seg << ", integral="
		    << adc << ", ct=" << cmt << std::endl;	
	  const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	  evDisp.ShowHitBGO(seg, adc);
	}    
      }
      
    }
  }

  // PiID counter
  {
    const HodoRHitContainer &cont = rawData->GetPiIDRawHC();
    int nhit = cont.size();

    for(int i = 0; i<nhit; ++i){	
      HodoRawHit *hit = cont.at(i);
      int seg = hit->SegmentId(); 
      int NhitTl = hit->GetSizeTdcUp();
      int NhitTt = hit->GetSizeTdcTUp();
	
      for(int m = 0; m<NhitTl; ++m){
	int tdc = hit->GetTdcUp(m);	
	HF2 (200, seg, tdc);
      }
      for(int m = 0; m<NhitTt; ++m){
	int tdc = hit->GetTdcTUp(m);	
	HF2 (201, seg, tdc);
      }      
    }
  }
  hodoAna->DecodePiIDHits(rawData);
  int nhit = hodoAna->GetNHitsPiID();          
  for(int i = 0; i<nhit; ++i){
    const FiberHit* hit = hodoAna->GetHitPiID(i);
    int mhit = hit->GetNumOfHit();
    int seg = hit->PairId();
    for(int m = 0; m<mhit; ++m){	
      double ctime  = hit->GetCTime(m);
      HF2 (202, seg, ctime);

      if (FlagEvDisp && ctime>-50 && ctime<50) {
	std::cout << "PiID : seg=" << seg << ", ct=" << ctime << std::endl;	
	const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	evDisp.ShowHitPiID(seg);
      }    

    }
  }


  // CFT   
  {
    for(int p = 0; p<NumOfPlaneCFT; ++p){
      int sumnhit = 0;
      int layer = p;

      double Maxpe = -10.;
      double Maxpeup = -10.;
      double Maxpedown = -10.;      
      int Maxpe_seg = -10;
      int Maxpe_segup = -10;
      int Maxpe_segdown = -10;
      double NofP = 0.;
      double newPeak = 0.;
          
      const HodoRHitContainer &cont = rawData->GetCFTRawHC(p);
      int nhit = cont.size();
      HF1 (1000*(p+1)+1, nhit);
      //event.nHitCFT[p] = nhit;
      
      for(int i = 0; i<nhit; ++i){	
	HodoRawHit *hit = cont.at(i);
	int seg = hit->SegmentId(); 
	int NhitT = hit->GetSizeTdcUp();
	//int NhitT_tr = hit->GetSizeTdcDown();
	int NhitT_tr = hit->GetSizeTdcTUp();
	int NhitAH = hit->GetSizeAdcUp();	  
	int NhitAL = hit->GetSizeAdcDown();
	
	//TDC
	int width = -999, bufTtr = -999;	      
	if(NhitT==NhitT_tr){
	  for(int m = 0; m<NhitT; ++m){	    
	    int bufT = hit->GetTdcUp(m);
	    int bufTtr = hit->GetTdcTUp(m);	      
	    HF2 (1000*(p+1)+100, seg, bufT);  //TDC Nhits 
	    HF2 (1000*(p+1)+101, seg, bufTtr);
	    
	    width = bufT - bufTtr;	  
	    HF2 (1000*(p+1)+104, seg, width);
	    //if(width>60){HF2 (1000*(p+1) +122, seg, bufT);}
	    //else{HF2 (1000*(p+1) +123, seg, bufT);}	    

	  }
	}else{
	  for(int m = 0; m<NhitT; ++m){	    
	    int bufT = hit->GetTdcUp(m);
	    HF2 (1000*(p+1)+100, seg, bufT);  //TDC Nhits 
	  }	  
	  for(int m = 0; m<NhitT_tr; ++m){
	    int bufTtr = hit->GetTdcTUp(m);	      
	    HF2 (1000*(p+1)+101, seg, bufTtr);
	  }	    
	}
		
	//ADC Hi
	for(int m = 0; m<NhitAH; ++m){
	  int bufAH = hit->GetAdcUp();	
	  HF2 (1000*(p+1)+200, seg, bufAH);
	}
	
	//ADC Low
	for(int m = 0; m<NhitAL; ++m){	    
	  int bufAL = hit->GetAdcDown();	
	  HF2 (1000*(p+1)+201, seg, bufAL);
	}
      }
    }
  }
  

  // FiberHitCFT    
  hodoAna->DecodeCFTHits(rawData);
  for(int p = 0; p<NumOfPlaneCFT; ++p){

    int nhit = hodoAna->GetNHitsCFT(p);          
    //if(nhit>20){ disp_flag = true;}
    int nhit_t = 0;
    for(int i = 0; i<nhit; ++i){
      const FiberHit* hit = hodoAna->GetHitCFT(p, i);
      int mhit = hit->GetNumOfHit();
      int seg_id = hit->PairId();
      double adcHi  = hit->GetAdcHi();
      double adcLow = hit->GetAdcLow();  
      double MIPLow = hit->GetMIPLow();  
      double dELow  = hit->GetdELow();  
      double CFT_r       = hit->GetPositionR();
      double CFT_phi   = hit->GetPositionPhi();	  
      /*
      int ig = seg_id / 16;
      int is = seg_id % 16;
      event.adc_raw_h[p][ig][is] = adcHi;
      event.adc_raw_l[p][ig][is] = adcLow;
      */
      event.adc_raw_h[p][seg_id] = adcHi;
      event.adc_raw_l[p][seg_id] = adcLow;

      if(mhit>0){
	HF1(1000*(p+1) +20, seg_id);
	HF2 (1000*(p+1) +202, seg_id, adcHi);
	HF2 (1000*(p+1) +203, seg_id, adcLow);
	HF2 (1000*(p+1) +205, seg_id, MIPLow);
	HF2 (1000*(p+1) +207, seg_id, dELow);	
      }

      bool fl_m = false;
      for(int m = 0; m<mhit; ++m){

	double leading = hit->GetLeading(m);	
	//double trailing = hit->GetTrailing(m);	
	double ctime  = hit->GetCTime(m);
	double time   = hit->GetTime(m);
	double width  = hit->GetWidth(m);	

	HF2 (1000*(p+1) +102, seg_id, leading);
	HF2 (1000*(p+1) +103, seg_id, ctime);
	HF2 (1000*(p+1) +144, seg_id, width);
	
	if(width>60){
	  HF2(1000*(p+1) +122, seg_id, leading);
	  HF2(1000*(p+1) +243, seg_id, adcLow);
	}else{HF2 (1000*(p+1) +123, seg_id, leading);}
	

	if(-30 < ctime && ctime < 30){
	  if(fl_m==false){
	    fl_m=true;
	    nhit_t++;
	    HF1(1000*(p+1)  +21, seg_id);
	    HF2( 3,CFT_r*cos(CFT_phi*Deg2Rad),CFT_r*sin(CFT_phi*Deg2Rad));
	    HF2(1000*(p+1) +222, seg_id, adcHi);
	    HF2(1000*(p+1) +233, seg_id, adcLow);
	    HF2(1000*(p+1) +255, seg_id, MIPLow);
	    HF2(1000*(p+1) +277, seg_id, dELow);	

	    if (FlagEvDisp) {
	      const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	      //evDisp.ShowHitFiber(p, seg_id, 25);
	      evDisp.ShowHitFiber(p, seg_id, adcLow);
	      /*
	      std::cout << "Fiber Hit : layer=" << p << ", seg=" << seg_id
			<< ", adcLow=" << adcLow << ", tdc=" << ctime 
			<< ", MIP=" << MIPLow << ", dE=" << dELow << ", width=" << width
			<< std::endl;
	      */
	    }		      
	    
	  }
	}

      }// mhit      
    }//nhit
    HF1 (1000*(p+1)+11, nhit_t);      
  }  

#if 1
  // Fiber Cluster
  for(int p = 0; p<NumOfPlaneCFT; ++p){
    hodoAna->TimeCutCFT(p, -30, 30); // CATCH@J-PARC  
    //hodoAna->AdcCutCFT(p, 0, 4000); // CATCH@J-PARC  
    hodoAna->AdcCutCFT(p, 50, 4000); // CATCH@J-PARC  for proton
    //hodoAna->WidthCutCFT(p, 60, 300); // pp scattering
    hodoAna->WidthCutCFT(p, 30, 300); // cosmic ray

    int ncl = hodoAna->GetNClustersCFT(p);
    HF1 (1000*(p+1)+12, ncl);      
    
    for(int i=0; i<ncl; ++i){
      FiberCluster *cl = hodoAna->GetClusterCFT(p,i);
      if(!cl) continue;
      double size  = cl->ClusterSize();
      double seg = cl->MeanSeg();
      double ctime = cl->CMeanTime();
      //double width = -cl->minWidth();
      double width = cl->Width();
      double sumADC= cl->SumAdcLow();
      double sumMIP= cl->SumMIPLow();
      double sumdE = cl->SumdELow();
      double r     = cl->MeanPositionR();
      double phi   = cl->MeanPositionPhi();
      HF2(4, r*cos(phi*Deg2Rad),r*sin(phi*Deg2Rad));      

      if (FlagEvDisp) {
	/*
	std::cout << "Cluster : layer=" << p << ", seg=" << seg
		  << ", size=" << size 
		  << ", sumADC=" << sumADC
		  << ", sumMIP=" << sumMIP
		  << ", sumdE=" << sumdE
		  << ", minWidth=" << width
		  << std::endl;
	*/
      }

    }
    
  }
#endif


  // CFT tracking
#if 1
  DCAna->DecodeCFTHits( rawData );
  DCAna->TrackSearchCFT();

  int ntCFT=DCAna->GetNtracksCFT();// vtx limit ver.
  event.ntCFT = ntCFT;
  if(ntCFT>=3)disp_flag = true;	
  DCLocalTrack *tpp[2];
  for( int i=0; i<ntCFT; ++i ){

    DCLocalTrack *tp=DCAna->GetTrackCFT(i);
    if(i==0){tpp[0]=tp;}
    else if(i==1){tpp[1]=tp;}

    int nh   = tp->GetNHit();
    int nhUV = tp->GetNHitUV();
    double chisqrXY=tp->GetChiSquareXY();
    double chisqrXYZ=tp->GetChiSquareZ();
    double vtx_z =tp->GetVtxZ();
    double theta =tp->GetThetaCFT();
    int xyFlag = tp->GetCFTxyFlag();
    int zFlag  = tp->GetCFTzFlag() ;
       
    ThreeVector Pos0 = tp->GetPos0();
    ThreeVector Dir = tp->GetDir();
    double A=(Dir.x()*Dir.x()+Dir.y()*Dir.y());

    // aka 
    double D=(Dir.x()*Dir.x()+Dir.y()*Dir.y()+Dir.z()*Dir.z());

    double phi = -999.;
    if(Dir.x()>=0 && Dir.y()>=0){
      phi = acos(Dir.x()/sqrt(A))*Rad2Deg;
    }//0~90
    else if (Dir.x()<0 && Dir.y()>=0){
      phi = acos(Dir.x()/sqrt(A))*Rad2Deg;
    }//90~180
    else if (Dir.x()<0 && Dir.y()<0){
      phi = 360. - acos(Dir.x()/sqrt(A))*Rad2Deg;
    }//180~270
    else if (Dir.x()>=0 && Dir.y()<0){
      phi = 360. - acos(Dir.x()/sqrt(A))*Rad2Deg; ;
    }//270~360
    else{}

    HF1(5, vtx_z); 
    HF1(6, theta);
    event.theta[i] = theta;
    event.phi[i] = phi;

    event.Pos[i] = Pos0;
    event.Dir[i] = Dir;

    // vertex
    ThreeVector  bPos(0., 0., 0.);
    ThreeVector  bdir(0., 0., 1.);
    double dist = 1000.;
    if(Pos0.x()>-500.){
      ThreeVector vtx = Kinematics::VertexPoint3D(bPos, Pos0, bdir, Dir, dist);
      event.vtx_x[i] = vtx.x();
      event.vtx_y[i] = vtx.y();
      event.vtx_z[i] = vtx.z();
    }

    event.TotaldE[i]   = tp->GetTotalSumdE()   ;
    event.TotaldEphi[i]= tp->GetTotalSumdEphi();
    event.TotaldEuv[i] = tp->GetTotalSumdEuv ();
    event.TotaldE_max[i]   = tp->GetTotalMaxdE()   ;
    event.TotaldEphi_max[i]= tp->GetTotalMaxdEphi();
    event.TotaldEuv_max[i] = tp->GetTotalMaxdEuv ();

    // straight layer
    for(int ip=0; ip<nh; ip++){
      DCLTrackHit *hit = tp->GetHit(ip);
      int layer = hit->GetLayer();
      int seg = (int)hit->GetMeanSeg();

      double phi_ini   = tp->GetPhiIni(layer);      
      double phi_track   = tp->GetPhiTrack(layer);      
      double z_track = tp->GetZTrack(layer);
      double dphi  = tp->GetdPhi(layer);
      event.seg[layer][i] = seg;
      event.phi_ini[layer][i] = phi_ini;
      event.phi_track[layer][i] = phi_track;
      event.dphi[layer][i] = dphi;
      event.z_track[layer][i] = z_track;
      HF2(1000*(layer+1)+300, z_track, dphi);

      event.adc[layer][i] = tp->GetSumAdc(layer);
      event.mip[layer][i] = tp->GetSumMIP(layer);
      event.dE[layer][i]  = tp->GetSumdE (layer);
      event.adc_max[layer][i] = tp->GetMaxAdc(layer);      
      event.mip_max[layer][i] = tp->GetMaxMIP(layer);
      event.dE_max[layer][i]  = tp->GetMaxdE (layer);

      if (FlagEvDisp) {
	std::cout << "track#" << i << ", layer=" << layer << ", seg=" << seg
		  << ", ini_phi=" << phi_ini << std::endl;
      }

    }

    // spiral layer
    double xmin, xmax;
    double ymin, ymax;
    for(int ip=0; ip<nhUV; ip++){
      DCLTrackHit *hit = tp->GetHitUV(ip);
      int layer = hit->GetLayer();
      int seg = (int)hit->GetMeanSeg();
      
      double phi_track   = tp->GetPhiTrack(layer);      
      double z_track = tp->GetZTrack(layer);
      double z_ini   = tp->GetZIni(layer);      
      double dz    = tp->GetdZ(layer);

      event.seg[layer][i] = seg;
      event.phi_track[layer][i] = phi_track;
      event.z_ini[layer][i] = z_ini;
      event.z_track[layer][i] = z_track;
      event.dz[layer][i] = dz;
      HF2(1000*(layer+1)+310, phi_track, dz);

      event.adc[layer][i] = tp->GetSumAdc(layer);
      event.mip[layer][i] = tp->GetSumMIP(layer);
      event.dE[layer][i]  = tp->GetSumdE (layer);
      event.adc_max[layer][i] = tp->GetMaxAdc(layer);      
      event.mip_max[layer][i] = tp->GetMaxMIP(layer);
      event.dE_max[layer][i]  = tp->GetMaxdE (layer);

      if (FlagEvDisp) {
	const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	double rr = 49+10*(layer/2);
	double xz = rr*cos(phi_track*math::Deg2Rad());
	double yz = rr*sin(phi_track*math::Deg2Rad());
	evDisp.ShowHitPosZX(z_ini, xz, 25);
	evDisp.ShowHitPosZY(z_ini, yz, 25);
	if(ip==0){xmin=xz;ymin=yz;}
	else if(ip==nhUV-1){xmax=xz;ymax=yz;}

	std::cout << "track#" << i << ", layer=" << layer << ", seg=" << seg
		  << ", phi=" << phi_track << ", z_ini=" << z_ini << std::endl;      	
      }

    }    

    CFTParticle * CFTPart = new CFTParticle(tp, rawData);  
    int segBGOt = CFTPart->GetTrackBGOSeg();// BGO track segment
    event.segBGOt[i] = segBGOt;

    // 2 track ver.
    if (FlagEvDisp) {
      if(chisqrXY<150){
	double Axy = tp->GetAxy(); double Bxy = tp->GetBxy();
	double Az  = tp->GetAz() ; double Bz  = tp->GetBz();

	double x1=-120,x2=120;
	double y1=-120,y2=120;
	//double x1=0,x2=0;
	std::cout << "phi = " << event.phi[i] << std::endl;
	if( (phi>=0&&phi<=45)||(phi>=135&&phi<=225)||(phi>=315&&phi<=360) ){
	  // y = a*x + b
	  if( (phi>=0&&phi<=45)||(phi>=315&&phi<=360) ){ // +
	    x1=0;x2=120;
	  }else{ x1=0;x2=-120; }// -
	  y1 = x1*Axy+Bxy; y2 = x2*Axy+Bxy;
	}else{
	  // x = a*y + b
	  if( (phi>45&&phi<135) ){ // +
	    y1=0;y2=120;
	  }else{y1=0;y2=-120; }// -
	  x1 = y1*Axy+Bxy; x2 = y2*Axy+Bxy;
	}

	const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	evDisp.DrawTrackInXYPlane(x1, y1, x2, y2);
	printf("(%f,%f),(%f,%f)\n",x1,y1,x2,y2);      
	
	double zxDir = Dir.z()/Dir.x();
	double zyDir = Dir.z()/Dir.y();
	
	if(x2>0){x1=-20; x2=80.;}else{x1=20; x2=-80;} 	
	if(y2>0){y1=-20; y2=80.;}else{y1=20; y2=-80;} 
	
	double z1x = Pos0.z() + zxDir*(x1-Pos0.x());  
	double z2x = Pos0.z() + zxDir*(x2-Pos0.x());  
	double z1y = Pos0.z() + zyDir*(y1-Pos0.y());  
	double z2y = Pos0.z() + zyDir*(y2-Pos0.y());  	
	evDisp.DrawTrackInZXPlane(z1x, x1, z2x, x2);
	evDisp.DrawTrackInZYPlane(z1y, y1, z2y, y2);

	z1x = Pos0.z() + zxDir*(xmin-Pos0.x());  
	z2x = Pos0.z() + zxDir*(xmax-Pos0.x());  
	z1y = Pos0.z() + zyDir*(ymin-Pos0.y());  
	z2y = Pos0.z() + zyDir*(ymax-Pos0.y());  	
	evDisp.DrawTrackInZXPlane_(z1x, xmin, z2x, xmax);
	evDisp.DrawTrackInZYPlane_(z1y, ymin, z2y, ymax);
	
	std::cout << " Pos0 =(" << Pos0.x()
		  << ", " << Pos0.y() 
		  << ", " << Pos0.z() << ")"  << std::endl;
	std::cout << " Dir =(" << Dir.x() 
		  << ", " << Dir.y() 
		  << ", " << Dir.z() << ")"  << std::endl;
	std::cout << " theta = " << event.theta[i] 
		  << ", phi = " << event.phi[i] << std::endl;
	std::cout << "final: Axy = " << Axy
		  << ",Bxy = " << Bxy << ", phi = " << phi
		  << std::endl;  	
      }
    }

    
  }
  if(ntCFT>1){
    // vertex of 2 tracks
    ThreeVector PosA = event.Pos[0];
    ThreeVector PosB = event.Pos[1];
    ThreeVector DirA = event.Dir[0];
    ThreeVector DirB = event.Dir[1];
    double distAB = -999.;

    ThreeVector vtxAB = Kinematics::VertexPoint3D(PosA, PosB, DirA, DirB, distAB);
    event.vtxAB_x = vtxAB.x();
    event.vtxAB_y = vtxAB.y();
    event.vtxAB_z = vtxAB.z();
  }
#endif 


#if 0
  if(ntCFT>1){
    // conbine to 16 layers tracking    
    double d_phi = fabs(event.phi[0]-event.phi[1]);
    if(d_phi>0){
      DCAna->DecodeCFT16Hits(rawData,tpp[0],0); // 1st track
      DCAna->DecodeCFT16Hits(rawData,tpp[1],1); // 2nd track
      
      DCAna->TrackSearchCFT16();
      int ntCFT16=DCAna->GetNtracksCFT16();
      //if(ntCFT16>=1)disp_flag = true;	
      for( int it=0; it<ntCFT16; ++it ){

	DCLocalTrack *tp=DCAna->GetTrackCFT16(it);
	int nh   = tp->GetNHit();
	int nhUV = tp->GetNHitUV();
	
	// straight layer
	int i=0;
	for(int ip=0; ip<nh; ip++){
	  DCLTrackHit *hit = tp->GetHit(ip);
	  int layer16 = hit->GetLayer();
	  int seg = (int)hit->GetMeanSeg();	  
	  double phi_ini   = tp->GetPhiIni(layer16);      
	  double phi_track = tp->GetPhiTrack(layer16);      
	  double z_track   = tp->GetZTrack(layer16);
	  double dphi      = tp->GetdPhi(layer16);
	  
	  int layer = layer16;
	  if(layer<8){i=0;}
	  else if(layer>=8){layer-=8;i=1;}
	  event.seg16[layer][i]      = seg;
	  event.phi16_ini[layer][i]   = phi_ini;
	  event.phi16_track[layer][i] = phi_track;
	  event.dphi16[layer][i]      = dphi;
	  event.z16_track[layer][i]   = z_track;
	  
	  if (FlagEvDisp) {
	    std::cout << "16 layer tracking : i=" << ip << ", layer=" << layer16 << ", seg=" << seg
		      << ", ini_phi=" << phi_ini << std::endl;	    
	  }

	}	
	// spiral layer
	for(int ip=0; ip<nhUV; ip++){
	  DCLTrackHit *hit = tp->GetHitUV(ip);
	  int layer = hit->GetLayer();
	  int seg = (int)hit->GetMeanSeg();	  
	  double phi_track = tp->GetPhiTrack(layer);      
	  double z_track   = tp->GetZTrack(layer);
	  double z_ini     = tp->GetZIni(layer);      
	  double dz        = tp->GetdZ(layer);

	  if(layer<8){i=0;}
	  else if(layer>=8){layer-=8;i=1;}

	  event.seg16[layer][i] = seg;
	  event.phi16_track[layer][i] = phi_track;
	  if(phi_track<180){
	    event.phi16__track[layer][i] = phi_track;
	  }else{
	    event.phi16__track[layer][i] = phi_track+360;
	  }
	  event.z16_ini[layer][i] = z_ini;
	  event.z16_track[layer][i] = z_track;
	  event.dz16[layer][i] = dz;

	  if (FlagEvDisp) {
	    std::cout << "16 layer tracking layer=" << layer << ", seg=" << seg
		      << ", phi=" << phi_track << ", z_ini=" << z_ini << std::endl;	    
	  }

	}	     

	ThreeVector Pos0 = tp->GetPos0();
	ThreeVector Dir = tp->GetDir();
	double A=(Dir.x()*Dir.x()+Dir.y()*Dir.y());
	double D=(Dir.x()*Dir.x()+Dir.y()*Dir.y()+Dir.z()*Dir.z());

	// Event Display	
	if (FlagEvDisp) {
	  hddaq::cout << "Event Display 16, flag = " << FlagEvDisp << std::endl;
	  double Axy = tp->GetAxy(); double Bxy = tp->GetBxy();
	  double Az  = tp->GetAz() ; double Bz  = tp->GetBz();
	  int xyFlag = tp->GetCFTxyFlag();

	  double x1=-120,x2=120;
	  double y1=-120,y2=120;
	  if(xyFlag==0){
	    //y=ax+b
	    y1 = x1*Axy+Bxy; y2 = x2*Axy+Bxy;
	  }else{
	    //x=ay+b
	    x1 = y1*Axy+Bxy; x2 = y2*Axy+Bxy;
	  }
	  
	  const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	  evDisp.DrawTrackInXYPlane_(x1, y1, x2, y2);
	  printf("(%f,%f),(%f,%f)\n",x1,y1,x2,y2);      
	  std::cout << "final: Axy = " << Axy
		    << ",Bxy = " << Bxy //<< ", phi = " << phi
		    << std::endl;  		  

	  double zxDir = Dir.z()/Dir.x();
	  double zyDir = Dir.z()/Dir.y();	  
	  x1=-80; x2=80.;	
	  y1=-80; y2=80.;	  
	  double z1x = Pos0.z() + zxDir*(x1-Pos0.x());  
	  double z2x = Pos0.z() + zxDir*(x2-Pos0.x());  
	  double z1y = Pos0.z() + zyDir*(y1-Pos0.y());  
	  double z2y = Pos0.z() + zyDir*(y2-Pos0.y());  	
	  evDisp.DrawTrackInZXPlane__(z1x, x1, z2x, x2);
	  evDisp.DrawTrackInZYPlane__(z1y, y1, z2y, y2);
	}

      }
      
    }

  }
#endif


  // Trigger Flag
  {
    const HodoRHitContainer &cont=rawData->GetTrigRawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg = hit->SegmentId()+1;
      int tdc = hit->GetTdc1();
      if( tdc ){
	event.trigpat[i]      = seg;
	event.trigflag[seg-1] = tdc;
      }
    }
  }

#if 1
  if (FlagEvDisp) {
    const EvDispCFT & evDisp = EvDispCFT::GetInstance();    
    if(disp_flag ){
      evDisp.UpdateCanvas();
      evDisp.get_command();
    }
    evDisp.EndOfEvent();
  }    
#endif
  
  return true;
}

//______________________________________________________________________________
bool
EventCFT::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
void
EventCFT::InitializeEvent( void )
{
  event.evnum  = 0;
  event.nhits  = 0;
  event.unhits = 0;
  event.dnhits = 0;
  event.ncl    = 0;

  for( int it=0; it<NumOfSegTrig; ++it ){
    event.trigpat[it]  = -1;
    event.trigflag[it] = -1;
  }


  event.ntCFT  = -1;  
  for(int i = 0; i<MaxDepth; ++i){
    event.phi[i]  = -999.;
    event.theta[i]  = -999.;
    event.vtx_x[i]  = -999.;
    event.vtx_y[i]  = -999.;
    event.vtx_z[i]  = -999.;
    for(int j = 0; j<3; ++j){
      event.Pos[i][j]  = -999.;
      event.Dir[i][j]  = -999.;
    }
    
    for( int p=0; p<NumOfPlaneCFT; ++p ){      
      event.seg[p][i]       = -999;
      event.dphi[p][i]       = -999.;
      event.phi_ini[p][i]    = -999.;
      event.phi_track[p][i]  = -999.;
      event.dz[p][i]       = -999.;
      event.z_ini[p][i]    = -999.;
      event.z_track[p][i]  = -999.;

      event.adc[p][i] = -999.; event.adc_max[p][i] = -999.;
      event.mip[p][i] = -999.; event.mip_max[p][i] = -999.;
      event.dE[p][i]  = -999.; event.dE_max[p][i]  = -999.;
      /*
      for(int j = 0; j<60; ++j){// group
	event.pede[p][j] = -999.;
	for(int k = 0; k<16; ++k){// seg in the group
	  event.adc_raw_h[p][j][k] = -999.;
	  event.adc_raw_l[p][j][k] = -999.;
	}
      }
      */
      for(int j = 0; j<NumOfSegCFT_PHI4; ++j){
	event.adc_raw_h[p][j] = -999.;
	event.adc_raw_l[p][j] = -999.;
      }
    }
    event.TotaldE[i]    = -999.; event.TotaldE_max[i]    = -999.;
    event.TotaldEphi[i] = -999.; event.TotaldEphi_max[i] = -999.;
    event.TotaldEuv[i]  = -999.; event.TotaldEuv_max[i]  = -999.;
  }
  event.vtxAB_x = -999.;
  event.vtxAB_y = -999.; 
  event.vtxAB_z = -999.;

  for( int p=0; p<NumOfPlaneCFT; ++p ){      
    for(int i = 0; i<2; ++i){
      event.seg16[p][i]       = -999;
      event.dphi16[p][i]       = -999.;
      event.phi16_ini[p][i]    = -999.;
      event.phi16_track[p][i]  = -999.;
      event.phi16__track[p][i]  = -999.;
      event.dz16[p][i]       = -999.;
      event.z16_ini[p][i]    = -999.;
      event.z16_track[p][i]  = -999.;
    }
  }

  event.nhBGO  = -1; 
  for( int i=0; i<NumOfSegBGO; ++i ){      
    event.segBGO[i] = -1;
    event.segBGOt[i] = -1;
    event.adcbgo[i] = -999; 
    event.tdcbgo[i] = -999;
  }

}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventCFT;
}

//______________________________________________________________________________
namespace
{
  const int    NbinTdc = 1000;
  const double MinTdc  =    0.;
  const double MaxTdc  = 1000.;

  const int    NbinTot =  136;
  const double MinTot  =   -8.;
  const double MaxTot  =  128.;

  const int    NbinTime = 1000;
  const double MinTime  = -500.;
  const double MaxTime  =  500.;
}
//______________________________________________________________________________
bool
ConfMan:: InitializeHistograms( void )
{
  HB1( 1, "Status", 20, 0., 20. );
  HB1( 2, "analys flag", 20, 0., 20. );

  //BGO
  HB2( 100, "BGO TDC", 25,0,25, 2000,0,2000 );
  HB2( 110, "BGO Corrected TDC", 25,0,25, 200,-100,100 );
  HB2( 101, "BGO Integral", 25,0,25, 5000,0,100000 );
  HB2( 111, "BGO Integral w/ TDC", 25,0,25, 5000,0,100000 );

  //PiID
  HB2( 200, "PiID TDC Leading" , 33,0,33, 1024,0,1024 );
  HB2( 201, "PiID TDC Trailing", 33,0,33, 1024,0,1024 );
  HB2( 202, "PiID ctime" , 33,0,33, 200,-100,100 );

  //CFT
  HB2( 3, "Fiber Position", 400,-100,100,400,-100,100 );
  HB2( 4, "Fiber Mean Position", 400,-100,100,400,-100,100 );
  HB1( 5, "vertex z (CFT tracking)", 1000,-500,500 );
  HB1( 6, "theta (CFT tracking)", 180, 0 ,180 );
  for(int i=0; i<NumOfPlaneCFT; i++){
    std::ostringstream title[30];
   
    if(i%2 == 0){// spiral layer
      int layer = (int)i/2 +1;
      title[0] << "CFT UV"<< layer << " : N hit" ;
      title[1] << "CFT UV"<< layer << " : N hit w/ Tdc cut";
      title[2] << "CFT UV"<< layer << " : NCluster hit"    ;
      title[3] << "CFT UV"<< layer << " : Hit Pattern"     ;
      title[4] << "CFT UV"<< layer << " : Hit Pattern w/ Tdc cut";
      // TDC
      title[5] << "CFT UV"<< layer << " : Tdc(Leading) vs seg" ;  //100
      title[6] << "CFT UV"<< layer << " : Tdc(Trailing) vs seg";  //101
      // TDC Fiber Hit
      title[7] << " CFT UV"<< layer << " : Time(Leading) vs seg"; //102  
      title[8] << " CFT UV"<< layer << " : Time(Leading) vs seg (width>60)"; //122 
      title[25]<< " CFT UV"<< layer << " : Time(Leading) vs seg (width<=60)";//123  
      title[9] << " CFT UV"<< layer << " : CTime vs seg"; //103 
      title[10]<< " CFT UV"<< layer << " : width vs seg"; //104
      title[24]<< " CFT UV"<< layer << " : width vs seg (Fiber Hit)"; //144
      // ADC
      title[11]<< "CFT UV"<< layer << " : Adc(High) vs seg";//200
      title[12]<< "CFT UV"<< layer << " : Adc(Low)  vs seg";//201     
      // ADC Fiber Hit
      title[13] << "CFT UV"<< layer << " : Adc(High)-pedestal vs seg"            ;//202
      title[14] << "CFT UV"<< layer << " : Adc(Low)-pedestal vs seg"             ;//203 
      title[15] << "CFT UV"<< layer << " : Adc(High)-pedestal vs seg w/ Tdc cut" ;//222 
      title[16] << "CFT UV"<< layer << " : Adc(Low)-pedestal vs seg w/ Tdc cut"  ;//233 
      title[26] << "CFT UV"<< layer << " : Adc(Low)-pedestal vs seg w/ width cut";//243 
      title[17] << "CFT UV"<< layer << " : MIP(High) vs seg"                     ;//204
      title[18] << "CFT UV"<< layer << " : MIP(Low) vs seg"                      ;//205
      title[19] << "CFT UV"<< layer << " : MIP(Low) vs seg w/ Tdc cut"           ;//255
      title[20] << "CFT UV"<< layer << " : dE(Low) vs seg"                       ;//207
      title[21] << "CFT UV"<< layer << " : dE(Low) vs seg w/ Tdc cut"            ;//277
      // tracking
      title[22]<< "CFT UV"<< layer << " : dphi vs z";
      title[23]<< "CFT UV"<< layer << " : dz vs phi";
    }else if(i%2 == 1){// straight layer
      int layer = (int)i/2 +1;
      title[0] << "CFT Phi"<< layer << " : N hit" ;
      title[1] << "CFT Phi"<< layer << " : N hit w/ Tdc cut";
      title[2] << "CFT Phi"<< layer << " : NCluster hit" ;
      title[3] << "CFT Phi"<< layer << " : Hit Pattern"  ;
      title[4] << "CFT Phi"<< layer << " : Hit Pattern w/ Tdc cut";
      // TDC
      title[5] << "CFT Phi"<< layer << " : Tdc(Leading) vs seg" ;     //100
      title[6] << "CFT Phi"<< layer << " : Tdc(Trailing) vs seg";     //101
      // TDC Fiber Hit
      title[7] << " CFT Phi"<< layer << " : Time(Leading) vs seg" ;  //102
      title[8] << " CFT Phi"<< layer << " : Time(Leading) vs seg (width>60)";//122  
      title[25]<< " CFT Phi"<< layer << " : Time(Leading) vs seg (width<=60)";//123  
      title[9] << " CFT Phi"<< layer << " : CTime vs seg";  
      title[10]<< " CFT Phi"<< layer << " : width vs seg"; //104
      title[24]<< " CFT Phi"<< layer << " : width vs seg (Fiber Hit)"; //144
      // ADC
      title[11]<< "CFT Phi"<< layer << " : Adc(High) vs seg";//200
      title[12]<< "CFT Phi"<< layer << " : Adc(Low) vs seg" ;//201
      // ADC Fiber Hit
      title[13] << "CFT Phi"<< layer << " : Adc(High)-pedestal vs seg"           ; //202
      title[14] << "CFT Phi"<< layer << " : Adc(Low)-pedestal vs seg"            ; //203
      title[15] << "CFT Phi"<< layer << " : Adc(High)-pedestal vs seg w/ Tdc cut"; //222
      title[16] << "CFT Phi"<< layer << " : Adc(Low)-pedestal vs seg w/ Tdc cut" ; //233
      title[26] << "CFT Phi"<< layer << " : Adc(Low)-pedestal vs seg w/ width cut";//243 
      title[17] << "CFT Phi"<< layer << " : MIP(High) vs seg"                    ; //204
      title[18] << "CFT Phi"<< layer << " : MIP(Low) vs seg"                     ; //205
      title[19] << "CFT Phi"<< layer << " : MIP(Low) vs seg w/ Tdc cut"          ; //255
      title[20] << "CFT Phi"<< layer << " : dE(Low) vs seg"                      ; //207
      title[21] << "CFT Phi"<< layer << " : dE(Low) vs seg w/ Tdc cut"           ; //277
      // tracking
      title[22]<< "CFT Phi"<< layer << " : dphi vs z";
      title[23]<< "CFT Phi"<< layer << " : dz vs phi";     
    }

    HB1( 1000*(i+1)+1,  title[0].str().c_str(),50,0,50);
    HB1( 1000*(i+1)+11, title[1].str().c_str(),50,0,50);
    HB1( 1000*(i+1)+12, title[2].str().c_str(),50,0,50);
    HB1( 1000*(i+1)+20, title[3].str().c_str(), NumOfSegCFT[i]+20, -10, NumOfSegCFT[i]+10);
    HB1( 1000*(i+1)+21, title[4].str().c_str(), NumOfSegCFT[i]+20, -10, NumOfSegCFT[i]+10);

    //TDC
    HB2( 1000*(i+1)+100, title[5].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 1024,0,1024);
    HB2( 1000*(i+1)+101, title[6].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 1024,0,1024);    
    // Fiber Hit
    HB2( 1000*(i+1)+102, title[7].str().c_str(),  NumOfSegCFT[i], 0, NumOfSegCFT[i],1000,0,1000);
    HB2( 1000*(i+1)+122, title[8].str().c_str(),  NumOfSegCFT[i], 0, NumOfSegCFT[i],1000,0,1000);
    HB2( 1000*(i+1)+123, title[25].str().c_str(),  NumOfSegCFT[i], 0, NumOfSegCFT[i],1000,0,1000);
    HB2( 1000*(i+1)+103, title[9].str().c_str(),  NumOfSegCFT[i], 0, NumOfSegCFT[i],1000,-500,500);
    HB2( 1000*(i+1)+104, title[10].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i],150,0,150);
    HB2( 1000*(i+1)+144, title[24].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i],150,0,150);

    //ADC
    HB2( 1000*(i+1)+200, title[11].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 4096,0,4096);
    HB2( 1000*(i+1)+201, title[12].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 4096,0,4096);
    
    // Fiber Hit
    HB2( 1000*(i+1)+202, title[13].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 4010,-10,4000);
    HB2( 1000*(i+1)+203, title[14].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 1010,-10,1000);
    HB2( 1000*(i+1)+222, title[15].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 4010,-10,4000);
    HB2( 1000*(i+1)+233, title[16].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 1010,-10,1000);
    HB2( 1000*(i+1)+243, title[26].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 1010,-10,1000);
    HB2( 1000*(i+1)+204, title[17].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 440,-2,20);
    HB2( 1000*(i+1)+205, title[18].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 440,-2,20);
    HB2( 1000*(i+1)+255, title[19].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 440,-2,20);

    HB2( 1000*(i+1)+207, title[20].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 220,-10,100);
    HB2( 1000*(i+1)+277, title[21].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 220,-10,100);
    
    // tracking
    HB2( 1000*(i+1)+300, title[22].str().c_str(), 200,0,400,200,-5,5);
    HB2( 1000*(i+1)+310, title[23].str().c_str(), 180,0,360,200,-10,10);
    
  }

  HB1( 10000, "nhvalid", 10,0,10);

  //Tree
  HBTree( "tree","tree of Counter" );
  tree->Branch("evnum",    &event.evnum,    "evnum/I");
  tree->Branch("trigpat",   event.trigpat,  Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",  event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));
  tree->Branch("nhits",     &event.nhits,        "nhits/I");
  
  tree->Branch("ntCFT",     &event.ntCFT,    "ntCFT/I");
  tree->Branch("theta",     event.theta,    "theta[ntCFT]/D");
  tree->Branch("phi",       event.phi,      "phi[ntCFT]/D");
  tree->Branch("vtx_x",     event.vtx_x,    "vtx_x[ntCFT]/D");
  tree->Branch("vtx_y",     event.vtx_y,    "vtx_y[ntCFT]/D");
  tree->Branch("vtx_z",     event.vtx_z,    "vtx_z[ntCFT]/D");
  tree->Branch("vtxAB_x",   &event.vtxAB_x,  "vtxAB_x/D");
  tree->Branch("vtxAB_y",   &event.vtxAB_y,  "vtxAB_y/D");
  tree->Branch("vtxAB_z",   &event.vtxAB_z,  "vtxAB_z/D");
  
  tree->Branch("seg",      event.seg,     Form("seg[%d][%d]/I", NumOfPlaneCFT, MaxDepth ) );

  tree->Branch("dphi",     event.dphi,     Form("dphi[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("phi_ini",  event.phi_ini,  Form("phi_ini[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("phi_track",event.phi_track,Form("phi_track[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );

  tree->Branch("dz",      event.dz,      Form("dz[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("z_ini",   event.z_ini,   Form("z_ini[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("z_track", event.z_track, Form("z_track[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );

  tree->Branch("adc",     event.adc,      Form("adc[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("mip",     event.mip,      Form("mip[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("dE",      event.dE,       Form("dE[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("adc_max", event.adc_max,  Form("adc_max[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("mip_max", event.mip_max,  Form("mip_max[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("dE_max",  event.dE_max,   Form("dE_max[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );

  tree->Branch("TotaldE",    event.TotaldE,    "totaldE/D");
  tree->Branch("TotaldEphi", event.TotaldEphi, "totaldEphi/D");
  tree->Branch("TotaldEuv",  event.TotaldEuv,  "totaldEuv/D");
  tree->Branch("TotaldE_max",    event.TotaldE_max,    "totaldE_max/D");
  tree->Branch("TotaldEphi_max", event.TotaldEphi_max, "totaldEphi_max/D");
  tree->Branch("TotaldEuv_max",  event.TotaldEuv_max,  "totaldEuv_max/D");

  // pedestal correction
  //tree->Branch("pede",       event.pede,       Form("adc[%d][%d]/D", NumOfPlaneCFT, 60 ) );
  tree->Branch("adc_raw_h",  event.adc_raw_h,  Form("adc_raw_h[%d][%d]/D", NumOfPlaneCFT, NumOfSegCFT_PHI4 ) );
  tree->Branch("adc_raw_l",  event.adc_raw_l,  Form("adc_raw_l[%d][%d]/D", NumOfPlaneCFT, NumOfSegCFT_PHI4 ) );

  // for 16 layer tracking
  tree->Branch("seg16",      event.seg16,      Form("seg16[%d][%d]/I", NumOfPlaneCFT, 2 ) );

  tree->Branch("dphi16",     event.dphi16,     Form("dphi16[%d][%d]/D", NumOfPlaneCFT, 2 ) );
  tree->Branch("phi16_ini",  event.phi16_ini,  Form("phi16_ini[%d][%d]/D", NumOfPlaneCFT, 2 ) );
  tree->Branch("phi16_track",event.phi16_track,Form("phi16_track[%d][%d]/D", NumOfPlaneCFT, 2 ) );
  tree->Branch("phi16__track",event.phi16__track,Form("phi16__track[%d][%d]/D", NumOfPlaneCFT, 2 ) );

  tree->Branch("dz16",      event.dz16,      Form("dz16[%d][%d]/D", NumOfPlaneCFT, 2 ) );
  tree->Branch("z16_ini",   event.z16_ini,   Form("z16_ini[%d][%d]/D", NumOfPlaneCFT, 2 ) );
  tree->Branch("z16_track", event.z16_track, Form("z16_track[%d][%d]/D", NumOfPlaneCFT, 2 ) );

  // BGO
  tree->Branch("nhBGO",     &event.nhBGO,    "nhBGO/I");
  tree->Branch("segBGO",    event.segBGO,    "segBGO[nhBGO]/I");
  tree->Branch("segBGOt",   event.segBGOt,   "segBGOt[nhBGO]/I");
  tree->Branch("adcBGO",    event.adcbgo,    "adcbgo[24]/D");
  tree->Branch("tdcBGO",    event.tdcbgo,    "adcbgo[24]/D");


  HPrint();
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")    &&
      InitializeParameter<HodoParamMan>("HDPRM") &&
      InitializeParameter<HodoPHCMan>("HDPHC")   &&
      InitializeParameter<UserParamMan>("USER")  );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
