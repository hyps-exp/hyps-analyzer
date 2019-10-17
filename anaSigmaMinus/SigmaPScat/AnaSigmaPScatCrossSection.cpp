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

#include "SigmaPScatParam.txt"
#include "SigmaBeam.txt"//FlightLen[NumOfMom];
#include "SigmaPScatAcceptance.txt"//Acceptance[NumOfMom][NbinCos],Acceptance_err[NumOfMom][NbinCos];
#include "SigmaPScatCutEfficiency.txt"//Efficiency[2][NumOfMomHist][NbinCos],Efficiency_err[2][NumOfMomHist][NbinCos];

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
  fChain->SetBranchStatus("thetaCM_SigmaPScat2npi",1);
  fChain->SetBranchStatus("cdistSigmaPScat",1);
  fChain->SetBranchStatus("cdistScatSigmaDecay2npi",1);
  fChain->SetBranchStatus("vtx_SigmaPScat",1);
  fChain->SetBranchStatus("vty_SigmaPScat",1);
  fChain->SetBranchStatus("vtz_SigmaPScat",1);
  fChain->SetBranchStatus("vtx_ScatSigmaDecay2npi",1);
  fChain->SetBranchStatus("vty_ScatSigmaDecay2npi",1);
  fChain->SetBranchStatus("vtz_ScatSigmaDecay2npi",1);
  fChain->SetBranchStatus("pikno_SigmaPScat2npi",1);
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
  TH1D *hdE_SigmaP[NumOfMomHist][NbinCos], *hdsdE_SigmaP[NumOfMomHist][NbinCos];
  TH1D *hdE_bg[NumOfMomHist][NbinCos], *hdsdE_bg[NumOfMomHist][NbinCos];
  for (int j=0; j<NumOfMomHist; j++) {
    for (int i=0; i<NbinCos; i++) {
	  hdsdE_SigmaP[j][i]  = new TH1D( Form("hdsdE%d_%d", j, i), Form("#SigmaP d#sigma/dEd#OmegadE %d", i), NbindE, -100., 100.);
	  hdE_SigmaP[j][i]  = new TH1D( Form("hdE%d_%d", j, i), Form("#SigmaP dE %d", i), NbindE, -100., 100.);
      hdsdE_bg[j][i]  = new TH1D( Form("hdsdE%d_bg_%d", j, i), Form("d#sigma/dEd#OmegadE", i), NbindE, -100., 100.);
      hdE_bg[j][i]  = new TH1D( Form("hdE%d_bg_%d", j, i), Form("dE %d", i), NbindE, -100., 100.);
    }
  }

  TH1D *h0_SigmaP[NumOfMomHist], *h1_SigmaP[NumOfMomHist], *h2_SigmaP[NumOfMomHist];
  TH1D *h0_bg[NumOfMomHist],*h1_bg[NumOfMomHist], *h2_bg[NumOfMomHist];
  for (int i=0; i<NumOfMomHist; i++) {
	h0_SigmaP[i]  = new TH1D( Form("hdE%d",i), Form("#Delta E (#SigmaP scattering) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbindE, -100., 100.);
	h1_SigmaP[i]  = new TH1D( Form("hdsdE%d",i), Form("d#sigma/dE (#SigmaP scattering) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbindE, -100., 100.);
	h2_SigmaP[i]  = new TH1D( Form("hdsdO%d",i), Form("d#sigma/d#Omega (#SigmaP scattering) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbinCos, -1., 1.);
	h2_SigmaP[i]->GetXaxis()->SetTitle("cos(#theta_{CM})");
	h2_SigmaP[i]->GetYaxis()->SetTitle("mb/sr");

	h0_bg[i]  = new TH1D( Form("hdE%d_bg",i), Form("#Delta E (back ground) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbindE, -100., 100.);
	h1_bg[i]  = new TH1D( Form("hdsdE%d_bg",i), Form("d#sigma/dE (back ground) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbindE, -100., 100.);
	h2_bg[i]  = new TH1D( Form("hdsdO%d_bg",i), Form("d#sigma/d#Omega (back ground) (%.2lf<p<%.2lf)", MinMom[i], MaxMom[i] ), NbinCos, -1., 1.);
	h2_bg[i]->GetXaxis()->SetTitle("cos(#theta_{CM})");
	h2_bg[i]->GetYaxis()->SetTitle("mb/sr");
  }

  double dsigma1_SigmaP[NumOfMomHist][NbinCos][NbindE], stat_err_dsigma1_SigmaP[NumOfMomHist][NbinCos][NbindE];
  double dsigma1_bg[NumOfMomHist][NbinCos][NbindE], stat_err_dsigma1_bg[NumOfMomHist][NbinCos][NbindE];
  double dsigma2_SigmaP[NumOfMomHist][NbinCos], stat_err_dsigma2_SigmaP[NumOfMomHist][NbinCos];
  double dsigma2_bg[NumOfMomHist][NbinCos], stat_err_dsigma2_bg[NumOfMomHist][NbinCos];
  double dsigma3_SigmaP[NumOfMomHist][NbindE], stat_err_dsigma3_SigmaP[NumOfMomHist][NbindE];
  double dsigma3_bg[NumOfMomHist][NbindE], stat_err_dsigma3_bg[NumOfMomHist][NbindE];
  double dsigma4_SigmaP[NumOfMomHist], stat_err_dsigma4_SigmaP[NumOfMomHist];
  double dsigma4_bg[NumOfMomHist], stat_err_dsigma4_bg[NumOfMomHist];
  for (int i=0; i<NumOfMomHist; i++) {
    for (int j=0; j<NbinCos; j++) {
      for (int k=0; k<NbindE; k++) {
        dsigma1_SigmaP[i][j][k] = 0.;
        stat_err_dsigma1_SigmaP[i][j][k] = 0.;
        dsigma1_bg[i][j][k] = 0.;
        stat_err_dsigma1_bg[i][j][k] = 0.;
		if( j==0 ){
        dsigma3_SigmaP[i][k] = 0.;
        stat_err_dsigma3_SigmaP[i][k] = 0.;
        dsigma3_bg[i][k] = 0.;
        stat_err_dsigma3_bg[i][k] = 0.;
		}
      }
	  dsigma2_SigmaP[i][j] = 0.;
	  stat_err_dsigma2_SigmaP[i][j] = 0.;
	  dsigma2_bg[i][j] = 0.;
	  stat_err_dsigma2_bg[i][j] = 0.;
    }
	dsigma4_SigmaP[i] = 0.;
	stat_err_dsigma4_SigmaP[i] = 0.;
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
  double dcos = h2_SigmaP[0]->GetBinWidth(1);
  double dOmega = 2.*3.14159265358*dcos;
  double coeff = 10000./0.424; // For Liq H2

  int nev = fChain->GetEntries();
  std::cout<<"[Total Entry:"<<nev<<"]"<<std::endl;

  for(int evenum =0; evenum<nev; ++evenum ){
	//for(int evenum =0; evenum<20000; ++evenum ){
	if( evenum%1000000==0) std::cout<<"["<<(long int)100*evenum/nev<<"/100]"<<std::endl;

	DeltaE_SigmaPScat2npi=-9999;     
	thetaCM_SigmaPScat2npi=-9999;     
	cdistSigmaPScat=-9999;     
	cdistScatSigmaDecay2npi=-9999;     
	vtx_SigmaPScat=-9999;     
	vty_SigmaPScat=-9999;     
	vtz_SigmaPScat=-9999;     
	vtx_ScatSigmaDecay2npi=-9999;     
	vty_ScatSigmaDecay2npi=-9999;     
	vtz_ScatSigmaDecay2npi=-9999;     
	pikno_SigmaPScat2npi=-1;
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

	TVector3 Vertex( vtx_KURAMA[pikno_SigmaPScat2npi], vty_KURAMA[pikno_SigmaPScat2npi], vtz_KURAMA[pikno_SigmaPScat2npi] );
	TVector3 VertexScat( vtx_SigmaPScat, vty_SigmaPScat, vtz_SigmaPScat );
	TVector3 VertexDecay( vtx_ScatSigmaDecay2npi, vty_ScatSigmaDecay2npi, vtz_ScatSigmaDecay2npi );
	bool target_flag=false;	
	if( getTargetFlag( Vertex ) ) target_flag=true;

	if( chisqrKuramapik[pikno_SigmaPScat2npi]<50 
		&& pKuramapik[pikno_SigmaPScat2npi]<0.89 
		&& qKuramapik[pikno_SigmaPScat2npi]>0 
		&& 0.15<ccm2[pikno_SigmaPScat2npi]&&ccm2[pikno_SigmaPScat2npi]<0.35 
		&& 1.16<MissMassCorrDE[pikno_SigmaPScat2npi]&&MissMassCorrDE[pikno_SigmaPScat2npi]<1.25 
		&& target_flag ){

	  bool dE_flag=false;
	  if( abs(DeltaE_SigmaPScat2npi)<100 ) dE_flag=true;
	  bool cdist1_flag=false;
	  if( cdistSigmaPScat<Limit_Real_cdistDecayPi ) cdist1_flag=true;
	  bool dz1_flag=false;
	  if( vtz_SigmaPScat-vtz_KURAMA[pikno_SigmaPScat2npi]>Limit_Real_dz1 ) dz1_flag=true;
	  bool cdist2_flag=false;
	  if( cdistScatSigmaDecay2npi<Limit_Real_cdistScat ) cdist2_flag=true;
	  bool dz2_flag=false;
	  if( vtz_ScatSigmaDecay2npi-vtz_SigmaPScat>Limit_Real_dz2 ) dz2_flag=true;
	  bool flength_flag=false;
	  if( GetFlightLength( Vertex, VertexScat, VertexDecay )<Limit_Real_flength ) flength_flag=true;
	  bool NPCut_flag=false;
	  if( abs(DeltaE_NPScat-PeakValue)>Limit_Real_NPCut ) NPCut_flag=true;
	  bool LambdaNCut_flag=false;
	  if( DeltaP_LambdaNConv<Limit_Real_LambdaNCut ) LambdaNCut_flag=true;
	  bool BG_flag=false;
	  if( 20<cdistSigmaPScat&&cdistSigmaPScat<40 ) BG_flag=true;
	
	  if( !dE_flag || ntProton==0 || ntOther==0 )	continue;      

      int indexMom = -1;
      for ( int i=0; i<NumOfMom; i++ ) {
		double low_range = MinMom[0]+0.02*i;
		double up_range  = low_range+0.02;
        if ( low_range < MissMomcal[pikno_SigmaPScat2npi]&&MissMomcal[pikno_SigmaPScat2npi] < up_range ) {
          indexMom = i;
          break;
        }
      }
      if ( !( 0 <= indexMom&&indexMom < NumOfMom ) ) continue;

      int indexMomHist = -1;
      for ( int i=0; i<NumOfMomHist; i++ ) {
        if ( MinMom[i] < MissMomcal[pikno_SigmaPScat2npi]&&MissMomcal[pikno_SigmaPScat2npi] < MaxMom[i] ) {
          indexMomHist = i;
          break;
        }
      }
      if ( !( 0 <= indexMomHist&&indexMomHist < NumOfMomHist ) ) continue;

	  double thetaCM = thetaCM_SigmaPScat2npi;
	  double dE = DeltaE_SigmaPScat2npi;
      int cosbin = h2_SigmaP[0]->FindBin(cos(thetaCM*Deg2Rad));
      int dEbin  = hdsdE_SigmaP[0][0]->FindBin(dE);
	  if( !( 0<cosbin&&cosbin<=NbinCos && 0<dEbin&&dEbin<NbindE ) ) continue;

      double eff = Acceptance[indexMom][cosbin-1];
      double cuteff1 = Efficiency[0][indexMomHist][cosbin-1];
      double cuteff2 = Efficiency[1][indexMomHist][cosbin-1];

      double Value1 = (1./(eff*cuteff1))/FL[indexMomHist][cosbin-1]*coeff*(1./dOmega);
      double Value2 = (1./(eff*cuteff2))/FL[indexMomHist][cosbin-1]*coeff*(1./dOmega);
      double Value3 = (1./(eff*cuteff1))/FL[indexMomHist][cosbin-1]*coeff;
      double Value4 = (1./(eff*cuteff2))/FL[indexMomHist][cosbin-1]*coeff;
	  if( FL[indexMomHist][cosbin-1]<1 || eff<EffThr || cuteff1<EffThr || cuteff2<EffThr ) continue;

      if( cdist1_flag && dz1_flag && cdist2_flag && dz2_flag && flength_flag && NPCut_flag ){
        dsigma1_SigmaP[indexMomHist][cosbin-1][dEbin-1] += Value1;
        stat_err_dsigma1_SigmaP[indexMomHist][cosbin-1][dEbin-1] += Value1*Value1;
        dsigma3_SigmaP[indexMomHist][dEbin-1] += Value3;
        stat_err_dsigma3_SigmaP[indexMomHist][dEbin-1] += Value3*Value3;
        
        hdsdE_SigmaP[indexMomHist][cosbin-1]->Fill(dE, Value1);
        hdE_SigmaP[indexMomHist][cosbin-1]->Fill(dE);
		h0_SigmaP[indexMomHist]->Fill(dE);		
		h1_SigmaP[indexMomHist]->Fill(dE,Value3);		
		if( abs(dE-PeakValue)<Limit_Real_CutdE ){
		  dsigma2_SigmaP[indexMomHist][cosbin-1] += Value2;
		  stat_err_dsigma2_SigmaP[indexMomHist][cosbin-1] += Value2*Value2;
		  dsigma4_SigmaP[indexMomHist] += Value4;
		  stat_err_dsigma4_SigmaP[indexMomHist] += Value4*Value4;
		  
		  h2_SigmaP[indexMomHist]->Fill(cos(thetaCM*Deg2Rad), Value2);
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
		if( abs(dE-PeakValue)<Limit_Real_CutdE ){
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
        stat_err_dsigma1_SigmaP[i][j][k] = sqrt(stat_err_dsigma1_SigmaP[i][j][k]);
		hdsdE_SigmaP[i][j]->SetBinError(k+1, stat_err_dsigma1_SigmaP[i][j][k]);
		stat_err_dsigma1_bg[i][j][k] = sqrt(stat_err_dsigma1_bg[i][j][k]);
        hdsdE_bg[i][j]->SetBinError(k+1, stat_err_dsigma1_bg[i][j][k]);
		if( j==0 ){
		  stat_err_dsigma3_SigmaP[i][k] = sqrt(stat_err_dsigma3_SigmaP[i][k]);
		  h1_SigmaP[i]->SetBinError(k+1, stat_err_dsigma3_SigmaP[i][k]);
		  stat_err_dsigma3_bg[i][k] = sqrt(stat_err_dsigma3_bg[i][k]);
		  h1_bg[i]->SetBinError(k+1, stat_err_dsigma3_bg[i][k]);
		}
	  }
	  stat_err_dsigma2_SigmaP[i][j] = sqrt(stat_err_dsigma2_SigmaP[i][j]);
	  h2_SigmaP[i]->SetBinError(j+1, stat_err_dsigma2_SigmaP[i][j]);
	  stat_err_dsigma2_bg[i][j] = sqrt(stat_err_dsigma2_bg[i][j]);
	  h2_bg[i]->SetBinError(j+1, stat_err_dsigma2_bg[i][j]);
    }
	stat_err_dsigma4_SigmaP[i] = sqrt(stat_err_dsigma4_SigmaP[i]);
	stat_err_dsigma4_bg[i] = sqrt(stat_err_dsigma4_bg[i]);
  }

  ///write rootfile////
  double X[NumOfMomHist], X_err[NumOfMomHist];
  for( int i=0; i<NumOfMomHist; i++ ){
	X[i] = (MaxMom[i]+MinMom[i])/2;
	X_err[i] = (MaxMom[i]-MinMom[i])/2;
  }
  TGraphErrors *h3_SigmaP = new TGraphErrors( NumOfMomHist, X, dsigma4_SigmaP, X_err, stat_err_dsigma4_SigmaP );
  h3_SigmaP->GetXaxis()->SetTitle("Momentum (GeV/c)");
  h3_SigmaP->GetYaxis()->SetTitle("Cross Section (mb)");
  h3_SigmaP->SetName("h3");
  TGraphErrors *h3_bg = new TGraphErrors( NumOfMomHist, X, dsigma4_bg, X_err, stat_err_dsigma4_bg );
  h3_bg->GetXaxis()->SetTitle("Momentum (GeV/c)");
  h3_bg->GetYaxis()->SetTitle("Cross Section (mb)");
  h3_bg->SetName("h3_bg");

  TFile *fout = new TFile("ForSigmaPScatCrossSection.root", "recreate");
  for(int i=0; i<NumOfMomHist; i++){
    for(int j=0; j<NbinCos; j++){
	  hdE_SigmaP[i][j]->Write();
	  hdsdE_SigmaP[i][j]->Write();
	  hdE_bg[i][j]->Write();
	  hdsdE_bg[i][j]->Write();
	}
	h0_SigmaP[i]->Write();
	h1_SigmaP[i]->Write();
	h2_SigmaP[i]->Write();
	h0_bg[i]->Write();
	h1_bg[i]->Write();
	h2_bg[i]->Write();
  }
  h3_SigmaP->Write();
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
