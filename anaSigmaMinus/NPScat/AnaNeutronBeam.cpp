#define DstPiKAna_SM_cxx
#include "DstPiKAna_SM.h"
#include <TH2.h>
#include <TH1.h>
#include <TF2.h>
#include <TFile.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLorentzVector.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string> 
#include "math.h"

bool getTargetFlag( TVector3 pos );
const int NumOfMom = 30;
const int NumOfMomHist = 3;
const  double MinMom[NumOfMomHist] = { 0.35, 0.45, 0.55 };
const  double MaxMom[NumOfMomHist] = { 0.45, 0.45, 0.65 };
#include "NPScatParam.txt"

// mass is MeV unit //
const double PionMass=0.1396;
const double KaonMass=0.4937;
const double ProtonMass=0.9383;
const double SigmaMass=1.1974;

const double TgtR = 20;//mm
const double TgtZ = 150;//mm

int main(){
  std::ifstream reading_file;
  std::string reading_line_buffer;
  reading_file.open("/home/had/kunpei/work/E40/analyzer2/ProductionRunList.txt");
  int RunList[1000];
  int run_list_number;
  int i=0;
  while( !reading_file.eof() ){
	std::getline(reading_file, reading_line_buffer); 
	if( reading_line_buffer[0]=='#' || reading_line_buffer.empty() ){ continue;
	}else{
	  std::istringstream line_separater(reading_line_buffer); 
	  RunList[i] = -1;
	  if( line_separater>>run_list_number ){
		RunList[i]=run_list_number;
		i++;
	  }
	}
  }
  RunList[i]=-9999;

  TChain *chain = new TChain("pik");
  for( int j=0; RunList[j]!=-9999; j++){
	TString rootfile = Form("/home/had/kunpei/work/E40/analyzer2/dst_rootfile/run0%d_DstPiKAna_SM_Sigma.root",RunList[j]);
	std::ifstream ifs(rootfile);
	if(ifs.is_open()){
	  chain->Add(rootfile);
	  std::cout<<"run0"<<RunList[j]<<"_DstPiKAna_SM_Sigma.root is reading..."<<std::endl;
	}
  }

  DstPiKAna_SM obj( chain );
  obj.Loop();

  return 0;
}

void DstPiKAna_SM::Loop(){
  if (fChain == 0) return;
  fChain->SetBranchStatus("*",0);  
  fChain->SetBranchStatus("NBeamMom",1);
  fChain->SetBranchStatus("NBeamMomx",1);
  fChain->SetBranchStatus("NBeamMomy",1);
  fChain->SetBranchStatus("NBeamMomz",1);
  fChain->SetBranchStatus("NBeamLength",1);
  fChain->SetBranchStatus("vtx_NBeam",1);
  fChain->SetBranchStatus("vty_NBeam",1);
  fChain->SetBranchStatus("vtz_NBeam",1);
  fChain->SetBranchStatus("vdist_NBeam",1);
  fChain->SetBranchStatus("cdistNBeam",1);
  fChain->SetBranchStatus("pitrno_NBeam",1);
  fChain->SetBranchStatus("pikno_NBeam",1);
  fChain->SetBranchStatus("vtz_KURAMA",1);//npik

  TH1D *h1 = new TH1D("h1","NeutronBeam Momentum",250,0,1.0);
  TH1D *h2 = new TH1D("h2","NeutronBeam Length(all)",320,0,320);
  TH1D *h22 = new TH1D("h22","NeutronBeam Length(all) vertex all",320,0,320);
  TH1D *h23 = new TH1D("h23","NeutronBeam Length(all) Length!=0",320,0,320);
  TH1D *h24 = new TH1D("h24","NeutronBeam Length(all) Up",500,0,50);
  TH1D *h3 = new TH1D("h3","NeutronBeam Length(0.35<p<0.45)",320,0,320);
  TH1D *h4 = new TH1D("h4","NeutronBeam Length(0.45<p<0.55)",320,0,320);
  TH1D *h5 = new TH1D("h5","NeutronBeam Length(0.55<p<0.65)",320,0,320);
  TH2D *h6 = new TH2D("h6","NeutronBeam Length vs SigmaDecayPosZ",320,0,320,200,-200,200);
  TH2D *vtxy = new TH2D("vtxy","vtx_SigmaDecay vs vty_SigmaDecay",250,-50,50,250,-50,50);
  TH1D *vtz = new TH1D("vtz","vtz_SigmaDecay",250,-500,500);
  TH1D *dz = new TH1D("dz","vtz_SigmaDecay - vtz_KURAMA",250,-100,100);
  TH1D *cdist = new TH1D("cdist","vtz_SigmaDecay",250,0,100);

  double FlightLen[NumOfMom];
  for( int i=0; i<NumOfMom; i++ ) FlightLen[i]=0;

  int nev = fChain->GetEntries();
  std::cout<<"[Total Entry:"<<nev<<"]"<<std::endl;

  for(int evenum =0; evenum<nev; ++evenum ){
  //for(int evenum =0; evenum<20000; ++evenum ){
	if( evenum%1000000==0) std::cout<<"["<<(long int)100*evenum/nev<<"/100]"<<std::endl;

	NBeamMom=-9999;     
	NBeamMomx=-9999;     
	NBeamMomy=-9999;     
	NBeamMomz=-9999;     
	NBeamLength=-9999;     
	vtx_NBeam=-9999;     
	vty_NBeam=-9999;     
	vtz_NBeam=-9999;     
	vdist_NBeam=-9999;     
	cdistNBeam=-9999;     
	pitrno_NBeam=-9999;     
	pikno_NBeam=-9999;   
	for( int i=0; i<20; i++ ){
	  vtz_KURAMA[i]=-9999;
	}
     	
    fChain->GetEntry(evenum);

	cdist->Fill(cdistNBeam);
	dz->Fill(vtz_NBeam-vtz_KURAMA[pikno_NBeam]);
	if( cdistNBeam>Limit_Real_cdistDecayPi || vtz_NBeam-vtz_KURAMA[pikno_NBeam]<Limit_Real_dz1 ) continue;

	h22->Fill(NBeamLength);
	if( NBeamLength!=0 ) h23->Fill(NBeamLength);

	vtxy->Fill(vtx_NBeam,vty_NBeam);
	vtz->Fill(vtz_NBeam);

	TVector3 pos( vtx_NBeam, vty_NBeam, vtz_NBeam );
	if( !getTargetFlag( pos ) ) continue;

	h1->Fill(NBeamMom);
	h2->Fill(NBeamLength);
	h24->Fill(NBeamLength);
	h6->Fill(NBeamLength,vtz_NBeam);
	if( MinMom[0]<=NBeamMom&&NBeamMom<MaxMom[0] ) h3->Fill(NBeamLength);
	if( MinMom[2]<=NBeamMom&&NBeamMom<MaxMom[1] ) h4->Fill(NBeamLength);
	if( MinMom[2]<=NBeamMom&&NBeamMom<MaxMom[2] ) h5->Fill(NBeamLength);

	for ( int i=0; i<NumOfMom; i++ ) {
	  double low_range = MinMom[0]+0.01*i;
	  double up_range  = low_range+0.01;
	  if ( low_range < NBeamMom&&NBeamMom < up_range ) {
		FlightLen[i]+=NBeamLength/10.;
		break;
	  }
	}

  }//event loop

  FILE *fp = fopen("NeutronBeam.txt","w");
  if(!fp) fprintf(stderr, "cannot open output file\n");
  
  fprintf(fp, "double FlightLen[NumOfMom] = { \n");
  for( int i=0; i<NumOfMom; i++ ){
	  fprintf( fp, " %.3lf", FlightLen[i] );
	  if(i != NumOfMom-1) fprintf(fp, ", \n");
  }
  fprintf(fp, "};\n");
    
  fclose(fp);

  ///write rootfile////
  TFile *fout = new TFile("ForNeutronBeam.root", "recreate");
  h1->Write();
  h2->Write();
  h22->Write();
  h23->Write();
  h24->Write();
  h3->Write();
  h4->Write();
  h5->Write();
  h6->Write();
  vtxy->Write();
  vtz->Write();
  dz->Write();
  cdist->Write();
  fout->Close();
}
//////////////////////////////////////////////////////////////////////////
bool getTargetFlag( TVector3 pos ){
  const double LH2TgtR = 40;//mm
  const double Target_Length = 300;//mm

  bool flagTgt = false;
  double r = sqrt((pos.X())*(pos.X())+(pos.Y())*(pos.Y()));
  if ( r <= LH2TgtR/2.&& abs(pos.Z()) <= Target_Length/2. ) {
    flagTgt = true;
  }
  return flagTgt;	
}
//////////////////////////////////////////////////////////////////////////
