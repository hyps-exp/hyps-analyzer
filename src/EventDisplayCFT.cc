/*
  EvDispCFT.cc

  2012/1/24
*/

//Unit is mm.

#include "EventDisplayCFT.hh"
#include "DCGeomMan.hh"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <sstream>
#include <string>
#include <algorithm>

#include "TemplateLib.hh"

#include "DCLocalTrack.hh"
#include "DCLTrackHit.hh"
//#include "BFTCluster.hh"
//#include "BFTHit.hh"

const int MaxChar = 200;

const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

EvDispCFT *EvDispCFT::evDisp_ = 0;
TApplication *EvDispCFT::theApp =0;

const int NumColor = 10;
const int Color[NumColor] = {kBlue, kRed, kCyan,  kGreen, kPink, kOrange, kYellow, kMagenta, kViolet, kSpring};

EvDispCFT::EvDispCFT(void)
{

}

EvDispCFT::~EvDispCFT(void)
{

}

EvDispCFT & EvDispCFT::GetInstance( void )
{
  if( !evDisp_ ){
    evDisp_ = new EvDispCFT();
  }
  if( !theApp ){
    theApp=new TApplication( "App", 0, 0 );
  }
  return *evDisp_;
}

void EvDispCFT::Initialize(int RunNum)
{
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  tc_ = new TCanvas("canvas","CFT Event Display",1600,800);

  tp_[0] = new TPad("pad0","PHI",0.001,0.001,0.499,0.999,10);
  tp_[0]->SetFrameFillColor(kGray);
  tp_[0]->Draw();

  tp_[5] = new TPad("pad5","",0.501,0.501,0.999,0.999,10);
  tp_[5]->Draw();

  tp_[6] = new TPad("pad6","",0.501,0.001,0.999,0.499,10);
  tp_[6]->Draw();

  tp_[0]->cd();
  gPad->SetGridx();  
  gPad->SetGridy();  

  hbase_ = new TH2F("hbase","Event Display XY plane", 120, -120, 120, 120, -120, 120);
  hbase_->SetMaximum(30);
  hbase_->SetMinimum(-1);
  hbase_->Draw();


  tp_[5]->cd();
  gPad->SetGridx();  
  gPad->SetGridy();  
  hbaseZX_ = new TH2F("hbaseZX", "Event Display ZX plane", 125, -100, 400, 100, -100, 100);
  hbaseZX_->SetMaximum(30);
  hbaseZX_->SetMinimum(-1);
  hbaseZX_->Draw();

  tp_[6]->cd();
  gPad->SetGridx();  
  gPad->SetGridy();  
  hbaseZY_ = new TH2F("hbaseZY", "Event Display ZY plane", 125, -100, 400, 100, -100, 100);
  hbaseZY_->SetMaximum(30);
  hbaseZY_->SetMinimum(-1);
  hbaseZY_->Draw();

  ConstructCFT();

  tc_->cd();
  tc_->Update();

}


void EvDispCFT::ConstructCFT(void)
{

  static const std::string funcname = "EvDispCFT::ConstructCFT";  

  tp_[0]->cd(); 

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  Tgt_Arc_ = new TArc(0, 0, 40./2);
  Tgt_Arc_->SetLineColor(kCyan);  
  Tgt_Arc_->SetFillStyle(0);
  Tgt_Arc_->Draw("same");

  CFRP_Arc_ = new TArc(0, 0, 80./2);
  CFRP_Arc_->SetLineColor(kBlack);  
  CFRP_Arc_->SetFillStyle(0);
  CFRP_Arc_->SetLineWidth(2);
  CFRP_Arc_->Draw("same");

  // UV1
  for (int seg=0; seg<NumOfSegCFT_UV1; seg++) {
    double x ,y;
    FiberPosPhi(0, seg, &x, &y); //UV1
    printf("UV1: seg=%d, x=%f, y=%f \n",seg,x,y);
    UV1_Arc_[seg] = new TArc(x, y, 0.75/2);
    UV1_Arc_[seg]->SetLineColor(kGray);

    UV1_Arc_[seg]->SetFillStyle(0);
    UV1_Arc_[seg]->Draw("same");
  }

  // UV2
  for (int seg=0; seg<NumOfSegCFT_UV2; seg++) {
    double x ,y;
    FiberPosPhi(2, seg, &x, &y); //UV2
    printf("UV2: seg=%d, x=%f, y=%f \n",seg,x,y);
    UV2_Arc_[seg] = new TArc(x, y, 0.75/2);
    UV2_Arc_[seg]->SetLineColor(kGray);

    UV2_Arc_[seg]->SetFillStyle(0);
    UV2_Arc_[seg]->Draw("same");
  }

  // UV3
  for (int seg=0; seg<NumOfSegCFT_UV3; seg++) {
    double x ,y;
    FiberPosPhi(4, seg, &x, &y); //UV3
    printf("UV3: seg=%d, x=%f, y=%f \n",seg,x,y);
    UV3_Arc_[seg] = new TArc(x, y, 0.75/2);
    UV3_Arc_[seg]->SetLineColor(kGray);

    UV3_Arc_[seg]->SetFillStyle(0);
    UV3_Arc_[seg]->Draw("same");
  }

  // UV4
  for (int seg=0; seg<NumOfSegCFT_UV4; seg++) {
    double x ,y;
    FiberPosPhi(6, seg, &x, &y); //UV4
    printf("UV4: seg=%d, x=%f, y=%f \n",seg,x,y);
    UV4_Arc_[seg] = new TArc(x, y, 0.75/2);
    UV4_Arc_[seg]->SetLineColor(kGray);

    UV4_Arc_[seg]->SetFillStyle(0);
    UV4_Arc_[seg]->Draw("same");
  }

  // PHI1
  for (int seg=0; seg<NumOfSegCFT_PHI1; seg++) {
    double x, y;
    double z=200;
    FiberPosPhi(1, seg, &x, &y);//PHI1
    printf("PHI1: seg=%d, x=%f, y=%f \n",seg,x,y);
    Phi1_Arc_[seg] = new TArc(x, y, 0.75/2);
    Phi1_Arc_[seg]->SetFillStyle(0);
    Phi1_Arc_[seg]->Draw("same");
  }

  // PHI2
  for (int seg=0; seg<NumOfSegCFT_PHI2; seg++) {
    double x, y;
    double z=200;
    FiberPosPhi(3, seg, &x, &y);
    printf("PHI2: seg=%d, x=%f, y=%f \n",seg,x,y);
    Phi2_Arc_[seg] = new TArc(x, y, 0.75/2);
    Phi2_Arc_[seg]->SetFillStyle(0);
    Phi2_Arc_[seg]->Draw("same");
  }

  // PHI3
  for (int seg=0; seg<NumOfSegCFT_PHI3; seg++) {
    double x, y;
    double z=200;
    FiberPosPhi(5, seg, &x, &y);//PHI3
    printf("PHI3: seg=%d, x=%f, y=%f \n",seg,x,y);
    Phi3_Arc_[seg] = new TArc(x, y, 0.75/2);
    Phi3_Arc_[seg]->SetFillStyle(0);
    Phi3_Arc_[seg]->Draw("same");
  }

  // PHI4
  for (int seg=0; seg<NumOfSegCFT_PHI4; seg++) {
    double x, y;
    double z=200;
    FiberPosPhi(7, seg, &x, &y);//PHI4
    printf("PHI4: seg=%d, x=%f, y=%f \n",seg,x,y);
    Phi4_Arc_[seg] = new TArc(x, y, 0.75/2);
    Phi4_Arc_[seg]->SetFillStyle(0);
    Phi4_Arc_[seg]->Draw("same");
  }


  // BGO
  double D2R = acos(-1.)/180.;
  for (int seg=0; seg<NumOfSegBGO; seg++) {

    double x1,x2,x3,x4,y1,y2,y3,y4;
    // box
#if 1 // new
    if(seg==0 || seg==1){      
      x1 = 100.0; x2 = 100.0 + 25.;
      if(seg==0){
	y1 = -30.0; y2 = 0.; 
      }else if(seg==1){
	y1 = 30.0; y2 = 0.; 
      }      
    }else if(seg==6 || seg==7){
      y1 = 100.0; y2 = 100.0 + 25.;
      if(seg==6){
	x1 = 30.0; x2 = 0.; 
      }else if(seg==7){
	x1 = -30.0; x2 = 0.; 
      }
    }else if(seg==12 || seg==13){
      x1 = -100.0; x2 = -100.0 - 25.;
      if(seg==12){
	y1 = 30.0; y2 = 0.; 
      }else if(seg==13){
	y1 = -30.0; y2 = 0.; 
      }
    }else if(seg==18 || seg==19){
      y1 = -100.0; y2 = -100.0 - 25.;
      if(seg==18){
	x1 = -30.0; x2 = 0.; 
      }else if(seg==19){
	x1 = 30.0; x2 = 0.; 
      }
    }else if(seg==3 || seg==4){ // Line
      double angle = 45.;
      x1 = 100.*cos(angle*D2R);
      y1 = 100.*sin(angle*D2R);
      x2 = x1 + 25.*cos(45.*D2R);
      y2 = y1 + 25.*cos(45.*D2R);      
      if(seg==4){
	x3 = x1 - 30.*cos(45.*D2R);
	y3 = y1 + 30.*cos(45.*D2R);
	x4 = x2 - 30.*cos(45.*D2R);
	y4 = y2 + 30.*cos(45.*D2R);	
      }else if(seg==3){
	x3 = x1 + 30.*cos(45.*D2R);
	y3 = y1 - 30.*cos(45.*D2R);
	x4 = x2 + 30.*cos(45.*D2R);
	y4 = y2 - 30.*cos(45.*D2R);	
      }
    }else if(seg==9 || seg==10){ // Line
      double angle = 135.;
      x1 = 100.*cos(angle*D2R);
      y1 = 100.*sin(angle*D2R);
      x2 = x1 - 25.*cos(45.*D2R);
      y2 = y1 + 25.*cos(45.*D2R);      
      if(seg==9){
	x3 = x1 + 30.*cos(45.*D2R);
	y3 = y1 + 30.*cos(45.*D2R);
	x4 = x2 + 30.*cos(45.*D2R);
	y4 = y2 + 30.*cos(45.*D2R);	
      }else if(seg==10){
	x3 = x1 - 30.*cos(45.*D2R);
	y3 = y1 - 30.*cos(45.*D2R);
	x4 = x2 - 30.*cos(45.*D2R);
	y4 = y2 - 30.*cos(45.*D2R);	
      }
    }else if(seg==15 || seg==16){ // Line
      double angle = -135.;
      x1 = 100.*cos(angle*D2R);
      y1 = 100.*sin(angle*D2R);
      x2 = x1 - 25.*cos(45.*D2R);
      y2 = y1 - 25.*cos(45.*D2R);      
      if(seg==16){
	x3 = x1 + 30.*cos(45.*D2R);
	y3 = y1 - 30.*cos(45.*D2R);
	x4 = x2 + 30.*cos(45.*D2R);
	y4 = y2 - 30.*cos(45.*D2R);	
      }else if(seg==15){
	x3 = x1 - 30.*cos(45.*D2R);
	y3 = y1 + 30.*cos(45.*D2R);
	x4 = x2 - 30.*cos(45.*D2R);
	y4 = y2 + 30.*cos(45.*D2R);	
      }
    }else if(seg==21 || seg==22){ // Line
      double angle = -45.;
      x1 = 100.*cos(angle*D2R);
      y1 = 100.*sin(angle*D2R);
      x2 = x1 + 25.*cos(45.*D2R);
      y2 = y1 - 25.*cos(45.*D2R);      
      if(seg==21){
	x3 = x1 - 30.*cos(45.*D2R);
	y3 = y1 - 30.*cos(45.*D2R);
	x4 = x2 - 30.*cos(45.*D2R);
	y4 = y2 - 30.*cos(45.*D2R);	
      }else if(seg==22){
	x3 = x1 + 30.*cos(45.*D2R);
	y3 = y1 + 30.*cos(45.*D2R);
	x4 = x2 + 30.*cos(45.*D2R);
	y4 = y2 + 30.*cos(45.*D2R);	
      }
    }else if(seg%3==2){ // Line
      double n=(seg+1)/3;
      double angle = +22.5+45.*(n-1);
      printf("BGO seg=%d, angle=%f\n", seg,angle);
      double xc = (120.+25./2.)*cos(angle*D2R);
      double yc = (120.+25./2.)*sin(angle*D2R);
      x1 = xc +(-25./2)*cos(angle*D2R) - (-30./2)*sin(angle*D2R);
      y1 = yc +(-25./2)*sin(angle*D2R) + (-30./2)*cos(angle*D2R);
      x2 = xc +(-25./2)*cos(angle*D2R) - (30./2)*sin(angle*D2R);
      y2 = yc +(-25./2)*sin(angle*D2R) + (30./2)*cos(angle*D2R);
      x3 = xc +(25./2)*cos(angle*D2R)  - (-30./2)*sin(angle*D2R);
      y3 = yc +(25./2)*sin(angle*D2R)  + (-30./2)*cos(angle*D2R);
      x4 = xc +(25./2)*cos(angle*D2R) - (30./2)*sin(angle*D2R);
      y4 = yc +(25./2)*sin(angle*D2R) + (30./2)*cos(angle*D2R);
    }
#endif

    if(seg%6==3 || seg%6==4 || seg%3==2){
      BGO_Line_[seg][0] = new TLine(x1, y1, x2, y2);
      BGO_Line_[seg][1] = new TLine(x1, y1, x3, y3);
      BGO_Line_[seg][2] = new TLine(x2, y2, x4, y4);
      BGO_Line_[seg][3] = new TLine(x3, y3, x4, y4);
      printf("BGO seg=%d, (%f,%f),(%f,%f),(%f,%f),(%f,%f)\n", seg,x1,y1,x2,y2,x3,y3,x4,y4);
      for(int i=0;i<4;i++){
	BGO_Line_[seg][i]->SetLineColor(kBlack);		
	BGO_Line_[seg][i]->Draw("same");
      }

    }else if(seg==0||seg==1){
    }else{
      BGO_Box_[seg] = new TBox(x1, y1, x2, y2);
      printf("BGO seg=%d, (x1,y1)=(%f,%f) , (x2,y2)=(%f,%f)\n", seg,x1,y1,x2,y2);
      BGO_Box_[seg]->SetLineColor(kBlack);
      
      BGO_Box_[seg]->SetFillStyle(0);
      BGO_Box_[seg]->Draw("same");
    }

  }

  tp_[5]->cd();
  Tgt_Box_ = new TBox(-50, -20, 300, 20);
  Tgt_Box_->SetLineColor(kCyan);  
  Tgt_Box_->SetFillStyle(0);
  Tgt_Box_->Draw("same");
  TgtHead_Arc_ = new TArc(300, 0, 40./2);
  TgtHead_Arc_->SetLineColor(kCyan);  
  TgtHead_Arc_->SetFillStyle(0);
  TgtHead_Arc_->Draw("same");

  CFRP_Box_ = new TBox(-100, -40, 390, 40);
  CFRP_Box_->SetLineColor(kBlack);  
  CFRP_Box_->SetFillStyle(0);
  CFRP_Box_->SetLineWidth(2);
  CFRP_Box_->Draw("same");
      
  tp_[6]->cd();
  Tgt_Box_->Draw("same");
  TgtHead_Arc_->Draw("same");
  CFRP_Box_->Draw("same");

}


void EvDispCFT::ShowHitFiber(int layer, int segment, double pe) const
{
  static const std::string funcname = "EvDispCFT::ShowHitFiber";
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  //printf("layer %d, seg %d, pe %f\n", layer, segment, pe);
  static int hit_unum=0;

  double x, y;
  if(layer%2==1){ //PHI
    FiberPosPhi(layer, segment, &x, &y);
  }else{   // UV
    FiberPosPhi(layer, segment, &x, &y);
  }
  //printf("plane= %d  (x,y)=(%f,%f) ,pe=%f\n",layer,x,y,pe);
  
  hbase_->Fill(x, y, pe);
  //hp_[layer]->Fill(segment, pe);  
}

void EvDispCFT::ShowHitPosZX(double Z, double X, double pe) const
{
  static const std::string funcname = "EvDispCFT::ShowHitFiberZX";
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  //printf("layer %d, seg %d, pe %f\n", layer, segment, pe);
  static int hit_unum=0;
  //hp_[5]->Fill(Z,X);
  hbaseZX_->Fill(Z,X,pe);
}

void EvDispCFT::ShowHitPosZY(double Z, double Y, double pe) const
{
  static const std::string funcname = "EvDispCFT::ShowHitFiberZY";
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  //printf("layer %d, seg %d, pe %f\n", layer, segment, pe);
  static int hit_unum=0;
  //hp_[6]->Fill(Z,Y);
  hbaseZY_->Fill(Z,Y,pe);      
}


void EvDispCFT::ShowHitBGO(int segment, int ADC) const
{
  static const std::string funcname = "EvDispCFT::ShowHitBGO";
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  //  printf("BGo seg=%d, ADC= %d\n", segment, ADC);
  static int hit_unum=0;

  //hbase_->Fill(x, y, ADC);
  if(segment%6==3||segment%6==4 || segment%3==2){
    for(int i=0;i<4;i++){
      BGO_Line_[segment][i]->SetLineColor(kRed);
    }
  }else if(segment==0||segment==1){
  }else{
    BGO_Box_[segment]->SetLineColor(kRed);
  } 
}

void EvDispCFT::ShowHitPos(double X, double Y, double pe) const
{
  static const std::string funcname = "EvDispCFT::ShowHitFiber";
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  //printf("layer %d, seg %d, pe %f\n", layer, segment, pe);

  static int hit_unum=0;
  hp_[0]->Fill(X,Y);
  hp_[2]->Fill((int)X, pe);     
}

void EvDispCFT::DrawTrackInXYPlane(double x0, double y0, double x1, double y1) const
{

  TLine *l = new TLine(x0, y0, x1, y1);
  l->SetLineColor(kYellow);
  //l->SetLineStyle(9);
  TrackXYCont.push_back(l);

}

void EvDispCFT::DrawTrackInZXPlane(double z0, double x0, double z1, double x1) const
{

  TLine *l = new TLine(z0, x0, z1, x1);
  //l->SetLineColor(kRed);
  l->SetLineColor(kYellow);
  //l->SetLineStyle(9);
  l->SetLineWidth(1);
  TrackZXCont.push_back(l);

}
void EvDispCFT::DrawTrackInZYPlane(double z0, double y0, double z1, double y1) const
{

  TLine *l = new TLine(z0, y0, z1, y1);
  //l->SetLineColor(kRed);
  l->SetLineColor(kYellow);
  //l->SetLineStyle(9);
  l->SetLineWidth(1);
  TrackZYCont.push_back(l);

}

void EvDispCFT::DrawTrackInXYPlane_(double x0, double y0, double x1, double y1) const
{

  TLine *l = new TLine(x0, y0, x1, y1);
  l->SetLineColor(kRed);
  TrackXYCont.push_back(l);

}

void EvDispCFT::DrawTrackInZXPlane_(double z0, double x0, double z1, double x1) const
{

  TLine *l = new TLine(z0, x0, z1, x1);
  l->SetLineColor(kBlack);
  l->SetLineWidth(2);
  TrackZXCont.push_back(l);

}
void EvDispCFT::DrawTrackInZYPlane_(double z0, double y0, double z1, double y1) const
{

  TLine *l = new TLine(z0, y0, z1, y1);
  l->SetLineColor(kBlack);
  l->SetLineWidth(2);
  TrackZYCont.push_back(l);

}

void EvDispCFT::DrawTrackInZXPlane__(double z0, double x0, double z1, double x1) const
{

  TLine *l = new TLine(z0, x0, z1, x1);
  l->SetLineColor(kRed);
  l->SetLineWidth(1);
  TrackZXCont.push_back(l);

}
void EvDispCFT::DrawTrackInZYPlane__(double z0, double y0, double z1, double y1) const
{

  TLine *l = new TLine(z0, y0, z1, y1);
  l->SetLineColor(kRed);
  l->SetLineWidth(1);
  TrackZYCont.push_back(l);

}


void EvDispCFT::UpdateCanvas() const
{

  tp_[0]->cd();
  //hbase_->Draw("colz");
  hbase_->Draw("col");
  //hbase_->Draw("");
  Tgt_Arc_->Draw("same");
  CFRP_Arc_->Draw("same");


  for (int seg=0; seg<NumOfSegCFT_UV1; seg++) {
    UV1_Arc_[seg]->Draw("same");
  }
  for (int seg=0; seg<NumOfSegCFT_PHI1; seg++) {
    Phi1_Arc_[seg]->Draw("same");
  }
  for (int seg=0; seg<NumOfSegCFT_UV2; seg++) {
    UV2_Arc_[seg]->Draw("same");
  }
  for (int seg=0; seg<NumOfSegCFT_PHI2; seg++) {
    Phi2_Arc_[seg]->Draw("same");
  }
  for (int seg=0; seg<NumOfSegCFT_UV3; seg++) {
    UV3_Arc_[seg]->Draw("same");
  }
  for (int seg=0; seg<NumOfSegCFT_PHI3; seg++) {
    Phi3_Arc_[seg]->Draw("same");
  }
  for (int seg=0; seg<NumOfSegCFT_UV4; seg++) {
    UV4_Arc_[seg]->Draw("same");
  }
  for (int seg=0; seg<NumOfSegCFT_PHI4; seg++) {
    Phi4_Arc_[seg]->Draw("same");
  }

  for (int seg=0; seg<NumOfSegBGO; seg++) {
    if(seg%6==3||seg%6==4||seg%3==2){
      for(int i=0;i<4;i++){
	BGO_Line_[seg][i]->Draw("same");
      }      
    }else if(seg==0||seg==1){
    }else{
      BGO_Box_[seg]->Draw("same");
    }
  }
  
  
  int nt = TrackXYCont.size();
  for (int i=0; i<nt; i++) {
    TLine *l = TrackXYCont[i];
    l->Draw("same");
  }

  tp_[5]->cd();
  hbaseZX_->Draw("col");
  Tgt_Box_->Draw("same");
  TgtHead_Arc_->Draw("same");
  CFRP_Box_->Draw("same");
  int ntZX = TrackZXCont.size();
  for (int i=0; i<ntZX; i++) {
    TLine *l = TrackZXCont[i];
    l->Draw("same");
  }

  tp_[6]->cd();
  hbaseZY_->Draw("col");
  Tgt_Box_->Draw("same");
  TgtHead_Arc_->Draw("same");
  CFRP_Box_->Draw("same");
  int ntZY = TrackZYCont.size();
  for (int i=0; i<ntZY; i++) {
    TLine *l = TrackZYCont[i];
    l->Draw("same");
  }


  tc_->cd();
  tc_->Update();
  tc_->Modified();


}

void EvDispCFT::EndOfEvent() const
{  

  //BGO
  for (int seg=0; seg<NumOfSegBGO; seg++) {
    if(seg%6==3||seg%6==4||seg%3==2){
      for(int i=0;i<4;i++){
	BGO_Line_[seg][i]->SetLineColor(kBlack);
      }
    }else if(seg==0||seg==1){
    }else{
      BGO_Box_[seg]->SetLineColor(kBlack);
    }
  }
  
  hbase_->Reset("ICES");
  hbaseZX_->Reset("ICES");
  hbaseZY_->Reset("ICES");
  
  std::for_each(TrackXYCont.begin(), TrackXYCont.end(), DeleteObject());
  TrackXYCont.clear();
  
  std::for_each(TrackZXCont.begin(), TrackZXCont.end(), DeleteObject());
  TrackZXCont.clear();
  
  std::for_each(TrackZYCont.begin(), TrackZYCont.end(), DeleteObject());
  TrackZYCont.clear();

}

void EvDispCFT::FiberPosPhi(int layer, int seg, double *x, double *y) const
{
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum=301+layer;

  double R     = geomMan.CalcCFTPositionR(lnum, seg);
  double Phi = geomMan.CalcCFTPositionPhi(lnum, seg) ;

  *x = R * cos(Phi*Deg2Rad);
  *y = R * sin(Phi*Deg2Rad);
  /*
  std::cout << "FiberPos layer" << layer << ", seg" << seg 
	    << "(" << R * cos(Theta*Deg2Rad) 
	    << "," << R * sin(Theta*Deg2Rad)
	    << ")" << std::endl;
  */
}

void EvDispCFT::FiberPosU(int layer, int seg, double z, double *x, double *y) const
{
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum=layer;
  double phi = -(360./NumOfSegCFT_UV2/*SegNumU*/)*(double)seg;;

  double r   = geomMan.GetLocalZ(lnum);
  if (seg%2 == 0)
    r -= 0.4755/2;
  else
    r += 0.4755/2;
  
  *x = r * cos(phi*Deg2Rad);
  *y = r * sin(phi*Deg2Rad);

}


void EvDispCFT::get_command(void) const
{
  tc_->Update();

  char ch;
  char data[100];
  static int stat=0;
  static int Nevent=0;
  static int ev=0;
  
  if (stat == 1 && Nevent > 0 && ev<Nevent) {

    ev++;
    return;
  } 
  if (ev==Nevent) {
    stat=0;
    ev=0;
  }

  if (stat == 0) {
    printf("q|n|p>");

    /* get command */
    scanf("%c",&ch);
    if (ch!='\n')
      while(getchar() != '\n');

    switch (ch) {
    case 'q': exit(0);
    case 'n':
      stat = 1;
      do {
	printf("event#>");
	scanf("%s",data);
      } while ((Nevent=atoi(data))<=0);
      std::cout << "Continue " << Nevent << "event" << std::endl;
      break;
    case 'p':
      theApp->Run(kTRUE);
      break;
    }
  }
}
