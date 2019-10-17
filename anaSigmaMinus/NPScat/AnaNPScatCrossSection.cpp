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
const double Deg2Rad = 3.14159265358/180;
const double Rad2Deg = 180/3.14159265358;

const int NbinCos=20;
const int NbindE=100;
const int NumOfMom = 30;
const int NumOfMomHist = 3;
const double EffThr = 0.05;

const  double MinMom[NumOfMomHist] = { 0.35, 0.45, 0.55 };
const  double MaxMom[NumOfMomHist] = { 0.45, 0.55, 0.65 };

#include "NPScatParam.txt"
#include "NeutronBeam.txt"//FlightLen[NumOfMom];
#include "NPScatAcceptance.txt"//Acceptance[NumOfMom][NbinCos],Acceptance_err[NumOfMom][NbinCos];
#include "NPScatCutEfficiency.txt"//Efficiency[2][NumOfMomHist][NbinCos],Efficiency_err[2][NumOfMomHist][NbinCos];

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
  fChain->SetBranchStatus("DeltaE_NPScat",1);
  fChain->SetBranchStatus("thetaCM_NPScat",1);
  fChain->SetBranchStatus("DecayNeutronMom",1);
  fChain->SetBranchStatus("cdistDecay2np",1);
  fChain->SetBranchStatus("cdistNPScat",1);
  fChain->SetBranchStatus("vtz_KURAMA",1);//npik
  fChain->SetBranchStatus("vtz_Decay2np",1);
  fChain->SetBranchStatus("vtz_NPScat",1);
  fChain->SetBranchStatus("pikno_NPScat",1);
  fChain->SetBranchStatus("ntOther",1);
  fChain->SetBranchStatus("ntProton",1);

  ///prepare hist////
  TH1D *hdE_NP[NumOfMomHist][NbinCos], *hdsdE_NP[NumOfMomHist][NbinCos];
  TH1D *hdE_bg[NumOfMomHist][NbinCos], *hdsdE_bg[NumOfMomHist][NbinCos];
  for (int j=0; j<NumOfMomHist; j++) {
    for (int i=0; i<NbinCos; i++) {
	  hdsdE_NP[j][i]  = new TH1D( Form("hdsdE%d_%d", j, i), Form("NP d#sigma/d#OmegadE %d", i), NbindE, -100., 100.);
	  hdE_NP[j][i]  = new TH1D( Form("hdE%d_%d", j, i), Form("NP dE %d", i), NbindE, -100., 100.);
      hdsdE_bg[j][i]  = new TH1D( Form("hdsdE%d_bg_%d", j, i), Form("d#sigma/d#OmegadE %d", i), NbindE, -100., 100.);
      hdE_bg[j][i]  = new TH1D( Form("hdE%d_bg_%d", j, i), Form("dE %d", i), NbindE, -100., 100.);
    }
  }

  TH1D *h0_NP[NumOfMomHist], *h1_NP[NumOfMomHist], *h2_NP[NumOfMomHist];
  TH1D *h0_bg[NumOfMomHist], *h1_bg[NumOfMomHist], *h2_bg[NumOfMomHist];
  for (int i=0; i<NumOfMomHist; i++) {
	h0_NP[i]  = new TH1D( Form("hdE%d",i), Form("#Delta E (NP scattering) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbindE, -100., 100.);
	h1_NP[i]  = new TH1D( Form("hdsdE%d",i), Form("d#sigma/dE (NP scattering) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbindE, -100., 100.);
	h2_NP[i]  = new TH1D( Form("hdsdO%d",i), Form("d#sigma/d#Omega (NP scattering) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbinCos, -1., 1.);
	h2_NP[i]->GetXaxis()->SetTitle("cos(#theta_{CM})");
	h2_NP[i]->GetYaxis()->SetTitle("mb/sr");

	h0_bg[i]  = new TH1D( Form("hdE%d_bg",i), Form("#Delta E (back ground) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbindE, -100., 100.);
	h1_bg[i]  = new TH1D( Form("hdsdE%d_bg",i), Form("d#sigma/dE (back ground) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbindE, -100., 100.);
	h2_bg[i]  = new TH1D( Form("hdsdO%d_bg",i), Form("d#sigma/d#Omega (back ground) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbinCos, -1., 1.);
	h2_bg[i]->GetXaxis()->SetTitle("cos(#theta_{CM})");
	h2_bg[i]->GetYaxis()->SetTitle("mb/sr");
  }
  
  double dsigma1_NP[NumOfMomHist][NbinCos][NbindE], stat_err_dsigma1_NP[NumOfMomHist][NbinCos][NbindE];
  double dsigma1_bg[NumOfMomHist][NbinCos][NbindE], stat_err_dsigma1_bg[NumOfMomHist][NbinCos][NbindE];
  double dsigma2_NP[NumOfMomHist][NbinCos], stat_err_dsigma2_NP[NumOfMomHist][NbinCos];
  double dsigma2_bg[NumOfMomHist][NbinCos], stat_err_dsigma2_bg[NumOfMomHist][NbinCos];
  double dsigma3_NP[NumOfMomHist][NbindE], stat_err_dsigma3_NP[NumOfMomHist][NbindE];
  double dsigma3_bg[NumOfMomHist][NbindE], stat_err_dsigma3_bg[NumOfMomHist][NbindE];
  double dsigma4_NP[NumOfMomHist], stat_err_dsigma4_NP[NumOfMomHist];
  double dsigma4_bg[NumOfMomHist], stat_err_dsigma4_bg[NumOfMomHist];
  for (int i=0; i<NumOfMomHist; i++) {
    for (int j=0; j<NbinCos; j++) {
      for (int k=0; k<NbindE; k++) {
        dsigma1_NP[i][j][k] = 0.;
        stat_err_dsigma1_NP[i][j][k] = 0.;
        dsigma1_bg[i][j][k] = 0.;
        stat_err_dsigma1_bg[i][j][k] = 0.;
		if( j==0 ){
		  dsigma3_NP[i][k] = 0.;
		  stat_err_dsigma3_NP[i][k] = 0.;
		  dsigma3_bg[i][j] = 0.;
		  stat_err_dsigma3_bg[i][k] = 0.;
		}
      }
	  dsigma2_NP[i][j] = 0.;
	  stat_err_dsigma2_NP[i][j] = 0.;
	  dsigma2_bg[i][j] = 0.;
	  stat_err_dsigma2_bg[i][j] = 0.;
    }
	dsigma4_NP[i] = 0.;
	stat_err_dsigma4_NP[i] = 0.;
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
	FL2[i] = 0.;
  }	
  for( int i=0; i<NumOfMomHist; i++ ){
    for( int j=0; j<NumOfMom; j++ ){
	  double mom = MinMom[0]+0.01*j+0.005;
      if ( MinMom[i] <mom&&mom< MaxMom[i] ){
        for( int k=0; k<NbinCos; k++ ){
          if( Acceptance[j][k]>EffThr ) FL[i][k] += FlightLen[j];
		}
		FL2[i] += FlightLen[j];
      }
    }
  }
  double dcos = h2_NP[0]->GetBinWidth(1);
  double dOmega = 2.*3.14159265358*dcos;
  double coeff = 10000./0.424; // For Liq H2

  int nev = fChain->GetEntries();
  std::cout<<"[Total Entry:"<<nev<<"]"<<std::endl;

  for(int evenum =0; evenum<nev; ++evenum ){
	//for(int evenum =0; evenum<20000; ++evenum ){
	if( evenum%100000==0) std::cout<<"["<<(long int)100*evenum/nev<<"/100]"<<std::endl;

	DeltaE_NPScat=-9999;     
	thetaCM_NPScat=-9999;     
	DecayNeutronMom=-9999;     
	cdistDecay2np=-9999;     
	cdistNPScat=-9999;     
	vtz_Decay2np=-9999;     
	vtz_NPScat=-9999;     
	pikno_NPScat=-1;   
	ntOther=-1;      
	ntProton=-1;   
	for( int i=0; i<20; i++ ){
	  vtz_KURAMA[i]=-9999;
	}
     	
    fChain->GetEntry(evenum);
	bool cdist1_flag=false;
	if( cdistDecay2np<Limit_Real_cdistDecayPi ) cdist1_flag=true;
	bool dz1_flag=false;
	if( vtz_Decay2np-vtz_KURAMA[pikno_NPScat]>Limit_Real_dz1 ) dz1_flag=true;
	bool target_flag=false;
	TVector3 DecayPiVert( vtx_Decay2np, vty_Decay2np, vtz_Decay2np );
	if( getTargetFlag( DecayPiVert ) ) target_flag=true;

    if( ntOther>0 && cdist1_flag && dz1_flag && target_flag ){//Neutron Beam required
      
      bool dE_flag=false;
	  if( abs(DeltaE_NPScat)<100 ) dE_flag=true;
      bool cdist2_flag=false;
	  if( cdistNPScat<Limit_Real_cdistNpScat ) cdist2_flag=true;
	  bool dz2_flag=false;
	  if( vtz_NPScat-vtz_Decay2np>Limit_Real_dz2 ) dz2_flag=true;
      bool BG_flag=false;
	  if( 20<cdistNPScat&&cdistNPScat<40 ) BG_flag=true;

	  if( !dE_flag || ntProton==0 )	continue;      
      
      int indexMom = -1;
      for ( int i=0; i<NumOfMom; i++ ) {
		double low_range = MinMom[0]+0.01*i;
		double up_range  = low_range+0.01;
        if ( low_range < DecayNeutronMom&&DecayNeutronMom < up_range ) {
          indexMom = i;
          break;
        }
      }
      if ( !( 0 <= indexMom&&indexMom < NumOfMom ) ) continue;

      int indexMomHist = -1;
      for ( int i=0; i<NumOfMomHist; i++ ) {
        if ( MinMom[i] < DecayNeutronMom&&DecayNeutronMom < MaxMom[i] ) {
          indexMomHist = i;
          break;
        }
      }
      if ( !( 0 <= indexMomHist&&indexMomHist < NumOfMomHist ) ) continue;

	  double thetaCM = thetaCM_NPScat;
	  double dE = DeltaE_NPScat;
      int cosbin = h2_NP[0]->FindBin(cos(thetaCM*Deg2Rad));
      int dEbin  = hdsdE_NP[0][0]->FindBin(dE);
	  if( !( 0<cosbin&&cosbin<=NbinCos && 0<dEbin&&dEbin<NbindE ) ) continue;

      double eff = Acceptance[indexMom][cosbin-1];
      double cuteff1 = Efficiency[0][indexMomHist][cosbin-1];
      double cuteff2 = Efficiency[1][indexMomHist][cosbin-1];

	  if( FL[indexMomHist][cosbin-1]<1 || eff<EffThr || cuteff1<EffThr || cuteff2<EffThr ) continue;
      double Value1 = (1./(eff*cuteff1))/FL[indexMomHist][cosbin-1]*coeff*(1./dOmega);
      double Value2 = (1./(eff*cuteff2))/FL[indexMomHist][cosbin-1]*coeff*(1./dOmega);
      double Value3 = (1./(eff*cuteff1))/FL2[indexMomHist]*coeff;
      double Value4 = (1./(eff*cuteff2))/FL2[indexMomHist]*coeff;
      if( cdist2_flag && dz2_flag ){
        dsigma1_NP[indexMomHist][cosbin-1][dEbin-1] += Value1;
        stat_err_dsigma1_NP[indexMomHist][cosbin-1][dEbin-1] += Value1*Value1;
        dsigma3_NP[indexMomHist][dEbin-1] += Value3;
        stat_err_dsigma3_NP[indexMomHist][dEbin-1] += Value3*Value3;
        
        hdsdE_NP[indexMomHist][cosbin-1]->Fill(dE, Value1);
        hdE_NP[indexMomHist][cosbin-1]->Fill(dE);
		h0_NP[indexMomHist]->Fill(dE);
		h1_NP[indexMomHist]->Fill(dE, Value3);
		if( abs(dE-PeakValue)<Limit_Real_CutdE ){
		  dsigma2_NP[indexMomHist][cosbin-1] += Value2;
		  stat_err_dsigma2_NP[indexMomHist][cosbin-1] += Value2*Value2;
		  dsigma4_NP[indexMomHist] += Value4;
		  stat_err_dsigma4_NP[indexMomHist] += Value4*Value4;
		  
		  h2_NP[indexMomHist]->Fill(cos(thetaCM*Deg2Rad), Value2);
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
		h1_bg[indexMomHist]->Fill(dE, Value3);
		if( abs(dE-PeakValue)<Limit_Real_CutdE ){
		  dsigma2_bg[indexMomHist][cosbin-1] += Value2;
		  stat_err_dsigma2_bg[indexMomHist][cosbin-1] += Value2*Value2;
		  dsigma4_bg[indexMomHist] += Value4;
		  stat_err_dsigma4_bg[indexMomHist] += Value4*Value4;
		  
		  h2_bg[indexMomHist]->Fill(cos(thetaCM*Deg2Rad), Value2);
		}
	  }

	}//NeutronBeam
  }//event loop

  for(int i=0; i<NumOfMomHist; i++){
    for(int j=0; j<NbinCos; j++){
      for(int k=0; k<NbindE; k++){
        stat_err_dsigma1_NP[i][j][k] = sqrt(stat_err_dsigma1_NP[i][j][k]);
		hdsdE_NP[i][j]->SetBinError(k+1, stat_err_dsigma1_NP[i][j][k]);
		stat_err_dsigma1_bg[i][j][k] = sqrt(stat_err_dsigma1_bg[i][j][k]);
        hdsdE_bg[i][j]->SetBinError(k+1, stat_err_dsigma1_bg[i][j][k]);
		if( j==0 ){
		  stat_err_dsigma3_NP[i][k] = sqrt(stat_err_dsigma3_NP[i][k]);
		  h1_NP[i]->SetBinError(k+1, stat_err_dsigma3_NP[i][k]);
		  stat_err_dsigma3_bg[i][k] = sqrt(stat_err_dsigma3_bg[i][k]);
		  h1_bg[i]->SetBinError(k+1, stat_err_dsigma3_bg[i][k]);
		}
	  }
	  stat_err_dsigma2_NP[i][j] = sqrt(stat_err_dsigma2_NP[i][j]);
	  h2_NP[i]->SetBinError(j+1, stat_err_dsigma2_NP[i][j]);
	  stat_err_dsigma2_bg[i][j] = sqrt(stat_err_dsigma2_bg[i][j]);
	  h2_bg[i]->SetBinError(j+1, stat_err_dsigma2_bg[i][j]);
    }
	stat_err_dsigma4_NP[i] = sqrt(stat_err_dsigma4_NP[i]);
	stat_err_dsigma4_bg[i] = sqrt(stat_err_dsigma4_bg[i]);
  }

  ///write rootfile////
  double X[NumOfMomHist], X_err[NumOfMomHist];
  for( int i=0; i<NumOfMomHist; i++ ){
	X[i] = (MaxMom[i]+MinMom[i])/2;
	X_err[i] = (MaxMom[i]-MinMom[i])/2;
  }
  TGraphErrors *h3_NP = new TGraphErrors( NumOfMomHist, X, dsigma4_NP, X_err, stat_err_dsigma4_NP );
  h3_NP->GetXaxis()->SetTitle("Momentum (GeV/c)");
  h3_NP->GetYaxis()->SetTitle("Cross Section (mb)");
  h3_NP->SetName("h3");
  TGraphErrors *h3_bg = new TGraphErrors( NumOfMomHist, X, dsigma4_bg, X_err, stat_err_dsigma4_bg );
  h3_bg->GetXaxis()->SetTitle("Momentum (GeV/c)");
  h3_bg->GetYaxis()->SetTitle("Cross Section (mb)");
  h3_bg->SetName("h3_bg");

  TFile *fout = new TFile("ForNPScatCrossSection.root", "recreate");
  for(int i=0; i<NumOfMomHist; i++){
    for(int j=0; j<NbinCos; j++){
	  hdE_NP[i][j]->Write();
	  hdsdE_NP[i][j]->Write();
	  hdE_bg[i][j]->Write();
	  hdsdE_bg[i][j]->Write();
	}
	h0_NP[i]->Write();
	h1_NP[i]->Write();
	h2_NP[i]->Write();
	h0_bg[i]->Write();
	h1_bg[i]->Write();
	h2_bg[i]->Write();
  }
  h3_NP->Write();
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
