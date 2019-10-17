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
double GetFlightLength( TVector3 v1, TVector3 v2, TVector3 v3 );
const double Deg2Rad = 3.14159265358/180;
const double Rad2Deg = 180/3.14159265358;

const int NbinCos=20;
const int NbindE=100;
const int NumOfMom = 20;
const int NumOfMomHist = 3;
const double EffThr = 0.05;

const  double MinMom[NumOfMomHist] = { 0.45, 0.55, 0.65 };
const  double MaxMom[NumOfMomHist] = { 0.55, 0.65, 0.85 };

#include "LambdaNConvParam.txt"
#include "SigmaBeam.txt"//FlightLen[NumOfMom];
#include "LambdaNConvAcceptance.txt"//Acceptance[NumOfMom][NbinCos],Acceptance_err[NumOfMom][NbinCos];
#include "LambdaNConvCutEfficiency.txt"//Efficiency[2][NumOfMomHist][NbinCos],Efficiency_err[2][NumOfMomHist][NbinCos];

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
  fChain->SetBranchStatus("DeltaE_SigmaPScat2npi",1);
  fChain->SetBranchStatus("DeltaE_NPScat",1);
  fChain->SetBranchStatus("DeltaP_LambdaNConv",1);
  fChain->SetBranchStatus("thetaCM_LamdaNConv",1);
  fChain->SetBranchStatus("cdistLambdaNConv",1);
  fChain->SetBranchStatus("cdistLambdaDecay",1);
  fChain->SetBranchStatus("vtz_LambdaNConv",1);
  fChain->SetBranchStatus("vtz_LambdaDecay",1);
  fChain->SetBranchStatus("pikno_LambdaNConv",1);
  fChain->SetBranchStatus("ntOther",1);
  fChain->SetBranchStatus("ntProton",1);
  fChain->SetBranchStatus("MissMassCorrDE",1);//npik 
  fChain->SetBranchStatus("chisqrKuramapik",1);//npik 
  fChain->SetBranchStatus("pKuramapik",1);//npik 
  fChain->SetBranchStatus("qKuramapik",1);//npik 
  fChain->SetBranchStatus("ccm2",1);//npik 
  fChain->SetBranchStatus("vtx_KURAMA",1);//npik 
  fChain->SetBranchStatus("vty_KURAMA",1);//npik 
  fChain->SetBranchStatus("vtz_KURAMA",1);//npik 
  fChain->SetBranchStatus("MissMomcal",1);//npik  

  ///prepare hist////
  TH1D *hdE_LN[NumOfMomHist][NbinCos], *hdsdE_LN[NumOfMomHist][NbinCos];
  TH1D *hdE_bg[NumOfMomHist][NbinCos], *hdsdE_bg[NumOfMomHist][NbinCos];
  for (int j=0; j<NumOfMomHist; j++) {
    for (int i=0; i<NbinCos; i++) {
	  hdsdE_LN[j][i]  = new TH1D( Form("hdsdE%d_%d", j, i), Form("#LambdaN d#sigma/dEd#OmegadE %d", i), NbindE, -0.5, 0.5);
	  hdE_LN[j][i]  = new TH1D( Form("hdE%d_%d", j, i), Form("#LambdaN dE %d", i), NbindE, -0.5, 0.5);
      hdsdE_bg[j][i]  = new TH1D( Form("hdsdE%d_bg_%d", j, i), Form("d#sigma/dEd#OmegadE", i), NbindE, -0.5, 0.5);
      hdE_bg[j][i]  = new TH1D( Form("hdE%d_bg_%d", j, i), Form("dE %d", i), NbindE, -0.5, 0.5);
    }
  }

  TH1D *h0_LN[NumOfMomHist], *h1_LN[NumOfMomHist], *h2_LN[NumOfMomHist];
  TH1D *h0_bg[NumOfMomHist],*h1_bg[NumOfMomHist], *h2_bg[NumOfMomHist];
  for (int i=0; i<NumOfMomHist; i++) {
	h0_LN[i]  = new TH1D( Form("hdE%d",i), Form("#Delta E (#LambdaN conversion) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbindE, -0.5, 0.5);
	h1_LN[i]  = new TH1D( Form("hdsdE%d",i), Form("d#sigma/dE (#LambdaN conversion) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbindE, -0.5, 0.5);
	h2_LN[i]  = new TH1D( Form("hdsdO%d",i), Form("d#sigma/d#Omega (#LambdaN conversion) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbinCos, -1., 1.);
	h2_LN[i]->GetXaxis()->SetTitle("cos(#theta_{CM})");
	h2_LN[i]->GetYaxis()->SetTitle("mb/sr");

	h0_bg[i]  = new TH1D( Form("hdE%d_bg",i), Form("#Delta E (back ground) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbindE, -0.5, 0.5);
	h1_bg[i]  = new TH1D( Form("hdsdE%d_bg",i), Form("d#sigma/dE (back ground) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbindE, -0.5, 0.5);
	h2_bg[i]  = new TH1D( Form("hdsdO%d_bg",i), Form("d#sigma/d#Omega (back ground) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbinCos, -1., 1.);
	h2_bg[i]->GetXaxis()->SetTitle("cos(#theta_{CM})");
	h2_bg[i]->GetYaxis()->SetTitle("mb/sr");
  }

  double dsigma1_LN[NumOfMomHist][NbinCos][NbindE], stat_err_dsigma1_LN[NumOfMomHist][NbinCos][NbindE];
  double dsigma1_bg[NumOfMomHist][NbinCos][NbindE], stat_err_dsigma1_bg[NumOfMomHist][NbinCos][NbindE];
  double dsigma2_LN[NumOfMomHist][NbinCos], stat_err_dsigma2_LN[NumOfMomHist][NbinCos];
  double dsigma2_bg[NumOfMomHist][NbinCos], stat_err_dsigma2_bg[NumOfMomHist][NbinCos];
  double dsigma3_LN[NumOfMomHist][NbindE], stat_err_dsigma3_LN[NumOfMomHist][NbindE];
  double dsigma3_bg[NumOfMomHist][NbindE], stat_err_dsigma3_bg[NumOfMomHist][NbindE];
  double dsigma4_LN[NumOfMomHist], stat_err_dsigma4_LN[NumOfMomHist];
  double dsigma4_bg[NumOfMomHist], stat_err_dsigma4_bg[NumOfMomHist];
  for (int i=0; i<NumOfMomHist; i++) {
    for (int j=0; j<NbinCos; j++) {
      for (int k=0; k<NbindE; k++) {
        dsigma1_LN[i][j][k] = 0.;
        stat_err_dsigma1_LN[i][j][k] = 0.;
        dsigma1_bg[i][j][k] = 0.;
        stat_err_dsigma1_bg[i][j][k] = 0.;
		if( j==0 ){
        dsigma3_LN[i][k] = 0.;
        stat_err_dsigma3_LN[i][k] = 0.;
        dsigma3_bg[i][k] = 0.;
        stat_err_dsigma3_bg[i][k] = 0.;
		}
      }
	  dsigma2_LN[i][j] = 0.;
	  stat_err_dsigma2_LN[i][j] = 0.;
	  dsigma2_bg[i][j] = 0.;
	  stat_err_dsigma2_bg[i][j] = 0.;
    }
	dsigma4_LN[i] = 0.;
	stat_err_dsigma4_LN[i] = 0.;
	dsigma4_bg[i] = 0.;
	stat_err_dsigma4_bg[i] = 0.;
  }
  ///end prepare hist////

  double FL[NumOfMomHist][NbinCos];
  double FL2[NumOfMomHist];
  for( int i=0; i<NumOfMomHist; i++ ){
    for( int j=0; j<NbinCos; j++ ){
      FL[i][j] = 0.;
	}
	FL2[i]=0;
  }	
  for( int i=0; i<NumOfMomHist; i++ ){
    for( int j=0; j<NumOfMom; j++ ){
	  double mom = MinMom[0]+0.02*j+0.01;
      if ( MinMom[i] <mom&&mom< MaxMom[i] ){
        for( int k=0; k<NbinCos; k++ ){
		  if( Acceptance[j][k]>EffThr ) FL[i][k] += FlightLen[j];
        }
		FL2[i] += FlightLen[j];
      }
    }
  }
  double dcos = h2_LN[0]->GetBinWidth(1);
  double dOmega = 2.*3.14159265358*dcos;
  double coeff = 10000./0.424; // For Liq H2

  int nev = fChain->GetEntries();
  std::cout<<"[Total Entry:"<<nev<<"]"<<std::endl;

  for(int evenum =0; evenum<nev; ++evenum ){
	//for(int evenum =0; evenum<20000; ++evenum ){
	if( evenum%1000000==0) std::cout<<"["<<(long int)100*evenum/nev<<"/100]"<<std::endl;

	DeltaE_SigmaPScat2npi=-9999;     
	thetaCM_LamdaNConv=-9999;     
	cdistLambdaNConv=-9999;     
	cdistLambdaDecay=-9999;     
	vtz_LambdaNConv=-9999;     
	vtz_LambdaDecay=-9999;     
	pikno_LambdaNConv=-1;
	DeltaE_NPScat=-9999;
	DeltaP_LambdaNConv=-9999;
	ntOther=-1;      
	ntProton=-1;   
	for( int i=0; i<20; i++ ){
	  chisqrKuramapik[i]=-9999.;
	  pKuramapik[i]=-9999.;
	  qKuramapik[i]=-9999.;
	  MissMassCorrDE[i]=-9999.;
	  ccm2[i]=-9999.;
	  MissMomcal[i]=-9999.;	  
	  vtx_KURAMA[i]=-9999;
	  vty_KURAMA[i]=-9999;
	  vtz_KURAMA[i]=-9999;
	}
     	
    fChain->GetEntry(evenum);

	TVector3 Vertex( vtx_KURAMA[pikno_LambdaNConv], vty_KURAMA[pikno_LambdaNConv], vtz_KURAMA[pikno_LambdaNConv] );
	bool target_flag=false;	
	if( getTargetFlag( Vertex ) ) target_flag=true;

	if( chisqrKuramapik[pikno_LambdaNConv]<50 
		&& pKuramapik[pikno_LambdaNConv]<0.89 
		&& qKuramapik[pikno_LambdaNConv]>0 
		&& 0.15<ccm2[pikno_LambdaNConv]&&ccm2[pikno_LambdaNConv]<0.35 
		&& 1.16<MissMassCorrDE[pikno_LambdaNConv]&&MissMassCorrDE[pikno_LambdaNConv]<1.25 
		&& target_flag ){

	  bool dE_flag=false;
	  if( abs(DeltaP_LambdaNConv)<0.5 ) dE_flag=true;
	  bool cdist1_flag=false;
	  if( cdistLambdaNConv<Limit_Real_cdistLambdaNConv ) cdist1_flag=true;
	  bool dz1_flag=false;
	  if( vtz_LambdaNConv-vtz_KURAMA[pikno_LambdaNConv]>Limit_Real_dz1 ) dz1_flag=true;
	  bool cdist2_flag=false;
	  if( cdistLambdaDecay<Limit_Real_cdistLambdaDecay ) cdist2_flag=true;
	  bool dz2_flag=false;
	  if( vtz_LambdaDecay-vtz_LambdaNConv>Limit_Real_dz2 ) dz2_flag=true;
	  bool NPCut_flag=false;
	  if( abs(DeltaE_NPScat-PeakValue)>Limit_Real_NPCut ) NPCut_flag=true;
	  bool SigmaPCut_flag=false;
	  if( abs(DeltaE_SigmaPScat2npi-PeakValue)>Limit_Real_SigmaPCut ) SigmaPCut_flag=true;
	  bool BG_flag=false;
	  if( 20<cdistLambdaNConv&&cdistLambdaNConv<40 ) BG_flag=true;
	
	  if( !dE_flag || ntProton==0 || ntOther==0 )	continue;      

      int indexMom = -1;
      for ( int i=0; i<NumOfMom; i++ ) {
		double low_range = MinMom[0]+0.02*i;
		double up_range  = low_range+0.02;
        if ( low_range < MissMomcal[pikno_LambdaNConv]&&MissMomcal[pikno_LambdaNConv] < up_range ) {
          indexMom = i;
          break;
        }
      }
      if ( !( 0 <= indexMom&&indexMom < NumOfMom ) ) continue;

      int indexMomHist = -1;
      for ( int i=0; i<NumOfMomHist; i++ ) {
        if ( MinMom[i] < MissMomcal[pikno_LambdaNConv]&&MissMomcal[pikno_LambdaNConv] < MaxMom[i] ) {
          indexMomHist = i;
          break;
        }
      }
      if ( !( 0 <= indexMomHist&&indexMomHist < NumOfMomHist ) ) continue;

	  double thetaCM = thetaCM_LamdaNConv;
	  double dE = DeltaP_LambdaNConv;
      int cosbin = h2_LN[0]->FindBin(cos(thetaCM*Deg2Rad));
      int dEbin  = hdsdE_LN[0][0]->FindBin(dE);
	  if( !( 0<cosbin&&cosbin<=NbinCos && 0<dEbin&&dEbin<NbindE ) ) continue;

      double eff = Acceptance[indexMom][cosbin-1];
      double cuteff1 = Efficiency[0][indexMomHist][cosbin-1];
      double cuteff2 = Efficiency[1][indexMomHist][cosbin-1];

      double Value1 = (1./(eff*cuteff1))/FL[indexMomHist][cosbin-1]*coeff*(1./dOmega);
      double Value2 = (1./(eff*cuteff2))/FL[indexMomHist][cosbin-1]*coeff*(1./dOmega);
      double Value3 = (1./(eff*cuteff1))/FL[indexMomHist][cosbin-1]*coeff;
      double Value4 = (1./(eff*cuteff2))/FL[indexMomHist][cosbin-1]*coeff;
	  if( FL[indexMomHist][cosbin-1]<1 || eff<EffThr || cuteff1<EffThr || cuteff2<EffThr ) continue;

      if( cdist1_flag && dz1_flag && cdist2_flag && dz2_flag && NPCut_flag ){
        dsigma1_LN[indexMomHist][cosbin-1][dEbin-1] += Value1;
        stat_err_dsigma1_LN[indexMomHist][cosbin-1][dEbin-1] += Value1*Value1;
        dsigma3_LN[indexMomHist][dEbin-1] += Value3;
        stat_err_dsigma3_LN[indexMomHist][dEbin-1] += Value3*Value3;

        hdsdE_LN[indexMomHist][cosbin-1]->Fill(dE, Value1);
        hdE_LN[indexMomHist][cosbin-1]->Fill(dE);
		h0_LN[indexMomHist]->Fill(dE);		
		h1_LN[indexMomHist]->Fill(dE,Value3);	
		if( fabs(dE)<Limit_Real_Cutdp ){
		  dsigma2_LN[indexMomHist][cosbin-1] += Value2;
		  stat_err_dsigma2_LN[indexMomHist][cosbin-1] += Value2*Value2;
		  dsigma4_LN[indexMomHist] += Value4;
		  stat_err_dsigma4_LN[indexMomHist] += Value4*Value4;
		  
		  h2_LN[indexMomHist]->Fill(cos(thetaCM*Deg2Rad), Value2);
		}
      }

      if( BG_flag ){
        dsigma1_bg[indexMomHist][cosbin-1][dEbin-1] += Value1;
        stat_err_dsigma1_bg[indexMomHist][cosbin-1][dEbin-1] += Value1*Value1;
        dsigma3_bg[indexMomHist][dEbin-1] += Value3;
        stat_err_dsigma3_bg[indexMomHist][dEbin-1] += Value3*Value3;
        
        hdsdE_bg[indexMomHist][cosbin-1]->Fill(dE, Value1);
        hdE_bg[indexMomHist][cosbin-1]->Fill(dE);
		h0_bg[indexMomHist]->Fill(dE);
		if( fabs(dE)<Limit_Real_Cutdp ){
		  dsigma2_bg[indexMomHist][cosbin-1] += Value2;
		  stat_err_dsigma2_bg[indexMomHist][cosbin-1] += Value2*Value2;
		  dsigma4_bg[indexMomHist] += Value4;
		  stat_err_dsigma4_bg[indexMomHist] += Value4*Value4;
		  		  
		  h2_bg[indexMomHist]->Fill(cos(thetaCM*Deg2Rad), Value2);
		}
	  }

	}//Beam
  }//event loop

  for(int i=0; i<NumOfMomHist; i++){
    for(int j=0; j<NbinCos; j++){
      for(int k=0; k<NbindE; k++){
        stat_err_dsigma1_LN[i][j][k] = sqrt(stat_err_dsigma1_LN[i][j][k]);
		hdsdE_LN[i][j]->SetBinError(k+1, stat_err_dsigma1_LN[i][j][k]);
		stat_err_dsigma1_bg[i][j][k] = sqrt(stat_err_dsigma1_bg[i][j][k]);
        hdsdE_bg[i][j]->SetBinError(k+1, stat_err_dsigma1_bg[i][j][k]);
		if( j==0 ){
		  stat_err_dsigma3_LN[i][k] = sqrt(stat_err_dsigma3_LN[i][k]);
		  h1_LN[i]->SetBinError(k+1, stat_err_dsigma3_LN[i][k]);
		  stat_err_dsigma3_bg[i][k] = sqrt(stat_err_dsigma3_bg[i][k]);
		  h1_bg[i]->SetBinError(k+1, stat_err_dsigma3_bg[i][k]);
		}
      }
	  stat_err_dsigma2_LN[i][j] = sqrt(stat_err_dsigma2_LN[i][j]);
	  h2_LN[i]->SetBinError(j+1, stat_err_dsigma2_LN[i][j]);
	  stat_err_dsigma2_bg[i][j] = sqrt(stat_err_dsigma2_bg[i][j]);
	  h2_bg[i]->SetBinError(j+1, stat_err_dsigma2_bg[i][j]);
    }
	stat_err_dsigma4_LN[i] = sqrt(stat_err_dsigma4_LN[i]);
	stat_err_dsigma4_bg[i] = sqrt(stat_err_dsigma4_bg[i]);
  }

  ///write rootfile////
  double X[NumOfMomHist], X_err[NumOfMomHist];
  for( int i=0; i<NumOfMomHist; i++ ){
	X[i] = (MaxMom[i]+MinMom[i])/2;
	X_err[i] = (MaxMom[i]-MinMom[i])/2;
  }
  TGraphErrors *h3_LN = new TGraphErrors( NumOfMomHist, X, dsigma4_LN, X_err, stat_err_dsigma4_LN );
  h3_LN->GetXaxis()->SetTitle("Momentum (GeV/c)");
  h3_LN->GetYaxis()->SetTitle("Cross Section (mb)");
  h3_LN->SetName("h3");
  TGraphErrors *h3_bg = new TGraphErrors( NumOfMomHist, X, dsigma4_bg, X_err, stat_err_dsigma4_bg );
  h3_bg->GetXaxis()->SetTitle("Momentum (GeV/c)");
  h3_bg->GetYaxis()->SetTitle("Cross Section (mb)");
  h3_bg->SetName("h3_bg");

  TFile *fout = new TFile("ForLambdaNConvCrossSection.root", "recreate");
  for(int i=0; i<NumOfMomHist; i++){
    for(int j=0; j<NbinCos; j++){
	  hdE_LN[i][j]->Write();
	  hdsdE_LN[i][j]->Write();
	  hdE_bg[i][j]->Write();
	  hdsdE_bg[i][j]->Write();
	}
	h0_LN[i]->Write();
	h1_LN[i]->Write();
	h2_LN[i]->Write();
	h0_bg[i]->Write();
	h1_bg[i]->Write();
	h2_bg[i]->Write();
  }
  h3_LN->Write();
  h3_bg->Write();
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
double GetFlightLength( TVector3 v1, TVector3 v2, TVector3 v3 ){
  TVector3 VertexLength1 = v2-v1;
  TVector3 VertexLength2 = v3-v2;
  double flightlength = VertexLength1.Mag()+VertexLength2.Mag();
  return flightlength;
}
