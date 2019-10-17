#define G4Simulation_cxx
#include "G4Simulation.h"
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
double getThetaCM( double m3, double m4, double cost, TLorentzVector PrimaryLv );
const double Deg2Rad = 3.14159265358/180;
const double Rad2Deg = 180/3.14159265358;

const int NbinCos = 20;
const int NumOfMomHist = 3;
const double MinMom[NumOfMomHist] = { 0.35, 0.45, 0.55 };
const double MaxMom[NumOfMomHist] = { 0.45, 0.55, 0.65 };
const int NbindE = 100;
#include "NPScatParam.txt"

// mass is MeV unit //
const double PionMass=0.1396;
const double KaonMass=0.4937;
const double ProtonMass=0.9383;
const double NeutronMass=0.9396;
const double SigmaMass=1.1974;

int main(){
  TChain *chain = new TChain("tree");
  // TString rootfile = Form("/home/had/kunpei/work/E40/analyzer2/dst_rootfile/run0%d_DstPiKAna_SM_Sigma.root",RunList[j]);
  TString rootfile = "G4Simulation.root";
  chain->Add(rootfile);

  G4Simulation obj( chain );
  obj.Loop();

  return 0;
}

void G4Simulation::Loop(){
  if (fChain == 0) return;
  fChain->SetBranchStatus("*",0);  
  fChain->SetBranchStatus("NNscatFlag",1);
  fChain->SetBranchStatus("momVectorProtonScat",1);//3
  fChain->SetBranchStatus("momVectorDecayNucleon",1);//3
  fChain->SetBranchStatus("momProtonScat",1);
  fChain->SetBranchStatus("momDecayNucleon",1);
  fChain->SetBranchStatus("nPi_CFT",1);
  fChain->SetBranchStatus("nP_CFT",1);
  fChain->SetBranchStatus("nSigma",1);
  fChain->SetBranchStatus("cdistDecayPi",1);
  fChain->SetBranchStatus("cdistNpScat",1);
  fChain->SetBranchStatus("Vertex_z",1);//nSigma
  fChain->SetBranchStatus("vertexDecayPi",1);//3
  fChain->SetBranchStatus("vertexNpScat",1);//3
  fChain->SetBranchStatus("EkinCorP",1);
  fChain->SetBranchStatus("scatNpEkinCal",1);

  TH1D *h0[NumOfMomHist][NbinCos], *h1[NumOfMomHist][NbinCos], *h2[NumOfMomHist][NbinCos], *h3[NumOfMomHist][NbinCos], *h4[NumOfMomHist][NbinCos];
  TH1D *hh0[NumOfMomHist], *hh1[NumOfMomHist], *hh2[NumOfMomHist], *hh3[NumOfMomHist], *hh4[NumOfMomHist], *hh5[NumOfMomHist], *hh6[NumOfMomHist];
  char buf[100], buf2[100];
  for(int i=0; i<NumOfMomHist; i++){
	for(int j=0; j<NbinCos; j++ ){
	  sprintf(buf, "h0_%d_%d", i,j);
	  sprintf(buf2, "dE (All) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	  h0[i][j] = new TH1D(buf, buf2, NbindE, -100, 100);
	  
	  sprintf(buf, "h1_%d_%d", i,j);
	  sprintf(buf2, "dE (Detected requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	  h1[i][j] = new TH1D(buf, buf2, NbindE, -100, 100);
	  
	  sprintf(buf, "h2_%d_%d", i,j);
	  sprintf(buf2, "dE (Detected && cdist requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	  h2[i][j] = new TH1D(buf, buf2, NbindE, -100, 100);
	  
	  sprintf(buf, "h3_%d_%d", i,j);
	  sprintf(buf2, "dE (Detected && dz requirement)(%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	  h3[i][j] = new TH1D(buf, buf2, NbindE, -100, 100);
	  
	  sprintf(buf, "h4_%d_%d", i,j);
	  sprintf(buf2, "dE (Detected && dz && cdist requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	  h4[i][j] = new TH1D(buf, buf2, NbindE, -100, 100);
	}
	sprintf(buf, "hh0_%d", i);
	sprintf(buf2, "costCM (All) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh0[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "hh1_%d", i);
	sprintf(buf2, "costCM (Detected requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh1[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "hh2_%d", i);
	sprintf(buf2, "costCM (Detected && cdist requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh2[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "hh3_%d", i);
	sprintf(buf2, "costCM (Detected && dz requirement)(%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh3[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "hh4_%d", i);
	sprintf(buf2, "costCM (Detected && dz && cdist requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh4[i] = new TH1D(buf, buf2, NbinCos, -1, 1);

	sprintf(buf, "hh5_%d", i);
	sprintf(buf2, "costCM (Efficiency1) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh5[i] = new TH1D(buf, buf2, NbinCos, -1, 1);

	sprintf(buf, "hh6_%d", i);
	sprintf(buf2, "costCM (Efficiency2) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh6[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
  }
  
  int nev = fChain->GetEntries();
  std::cout<<"[Total Entry:"<<nev<<"]"<<std::endl;
  
  for(int evenum =0; evenum<nev; ++evenum ){
	//for(int evenum =0; evenum<20000; ++evenum ){
	if( evenum%1000000==0) std::cout<<"["<<(long int)100*evenum/nev<<"/100]"<<std::endl;
	
	for( int i=0; i<3; i++ ){
	  momVectorProtonScat[i]=-9999;
	  momVectorDecayNucleon[i]=-9999;
	  vertexDecayPi[i]=-9999;
	  vertexNpScat[i]=-9999;
	}
	NNscatFlag=-1;
	momProtonScat=-9999;
	momDecayNucleon=-9999;
	cdistDecayPi=-9999;
	cdistNpScat=-9999;
	Vertex_z[0]=-9999;
	EkinCorP=-9999;
	scatNpEkinCal=-9999;
	nPi_CFT=-1;
	nP_CFT=-1;
	nSigma=-1;

	fChain->GetEntry(evenum);	
	if( NNscatFlag<0 ) continue;
	if( nSigma!=1 ) continue;
	
	TVector3 NBeamMom( momVectorDecayNucleon[0], momVectorDecayNucleon[1], momVectorDecayNucleon[2]);
	TVector3 ScatPMom( momVectorProtonScat[0], momVectorProtonScat[1], momVectorProtonScat[2]);
	TVector3 DecayPiVert( vertexDecayPi[0], vertexDecayPi[1], vertexDecayPi[2] );

	bool cdist_flag=false;
	if( cdistDecayPi<Limit_Sim_cdistDecayPi ) cdist_flag=true;
	bool dz_flag=false;
	if( vertexDecayPi[2]-Vertex_z[0]>Limit_Sim_dz1 ) dz_flag=true;
	bool target_flag=false;
	if( getTargetFlag( DecayPiVert ) ) target_flag=true;
      
	int index=-1;
	for ( int i=0; i<NumOfMomHist; i++ ) {
	  if ( MinMom[i] < NBeamMom.Mag()&&NBeamMom.Mag() < MaxMom[i] ){
		index = i;
		break;
	  }
	}
	if( !( 0 <= index&&index < NumOfMomHist ) ) continue;
	
	double thetaCM;
	double costLab = NBeamMom*ScatPMom/(NBeamMom.Mag()*ScatPMom.Mag());
	TLorentzVector LvN( NBeamMom, sqrt( NBeamMom.Mag()*NBeamMom.Mag()+NeutronMass*NeutronMass ));
	TLorentzVector LvP( 0, 0, 0, ProtonMass );
	TLorentzVector PrimaryLv = LvN + LvP;
	double thetaCM_P = getThetaCM( ProtonMass, NeutronMass, costLab, PrimaryLv );
	if( !( 0 <= thetaCM_P&&thetaCM_P <= 360) ) continue;
	else if( thetaCM_P+180<=360 ) thetaCM=thetaCM_P+180;
	else thetaCM = thetaCM_P-180;
	int cosbin=hh0[0]->FindBin(cos(thetaCM*Deg2Rad));

	double dE = EkinCorP - scatNpEkinCal;
	h0[index][cosbin-1]->Fill( dE );
	hh0[index]->Fill( cos(thetaCM*Deg2Rad) );
	if( nPi_CFT>0 && cdist_flag && dz_flag && target_flag && nP_CFT>0 && abs(dE)<Limit_Sim_DetectdE ){//Detected
	  h1[index][cosbin-1]->Fill( dE );
	  hh1[index]->Fill( cos(thetaCM*Deg2Rad) );
	  if( cdistNpScat<Limit_Sim_cdistNpScat ){
	  h2[index][cosbin-1]->Fill( dE );
	  hh2[index]->Fill( cos(thetaCM*Deg2Rad) );		
	  }
	  if( vertexNpScat[2]-vertexDecayPi[2]>Limit_Sim_dz2 ){
		h3[index][cosbin-1]->Fill( dE );
		hh3[index]->Fill( cos(thetaCM*Deg2Rad) );
	  }
	  if( cdistNpScat<Limit_Sim_cdistNpScat && vertexNpScat[2]-vertexDecayPi[2]>Limit_Sim_dz2 ){
		h4[index][cosbin-1]->Fill( dE );
		hh4[index]->Fill( cos(thetaCM*Deg2Rad) );	
   	  }
	}    
	
  }//event loop
	
  for(int i=0; i<NumOfMomHist; i++){
	for (int j=1; j<=NbinCos; j++){
	  double val1 = h4[i][j]->Integral(h4[i][j]->FindBin(-100),h4[i][j]->FindBin(100));
	  double val2 = h1[i][j]->Integral(h1[i][j]->FindBin(-100),h1[i][j]->FindBin(100));
	  if(val2 > 0){
		hh5[i]->SetBinContent(j, val1/val2);
		hh5[i]->SetBinError(j, 1./val2*sqrt(val1*(val2-val1)/val2));
	  }

	  val1 = h4[i][j]->Integral(h4[i][j]->FindBin(-Limit_Sim_CutdE),h4[i][j]->FindBin(Limit_Sim_CutdE));
	  if(val2 > 0){
		hh6[i]->SetBinContent(j, val1/val2);
		hh6[i]->SetBinError(j, 1./val2*sqrt(val1*(val2-val1)/val2));
	  }
	}
  }
    
  FILE *fp = fopen("NPScatCutEfficiency.txt","w");
  if(!fp) fprintf(stderr, "cannot open output file\n");
  
  fprintf(fp, "double Efficiency[2][NumOfMomHist][NbinCos] = { \n");
  for( int i=0; i<2; i++ ){
	fprintf(fp, "{ ");
	for( int j=0; j<NumOfMomHist; j++ ){
	  fprintf(fp, "{ ");
	  for( int k=1; k<=NbinCos; k++ ){
		if( i==0 ) fprintf( fp, " %.3lf", hh5[j]->GetBinContent(k) );
		if( i==1 ) fprintf( fp, " %.3lf", hh6[j]->GetBinContent(k) );
		if( k!=NbinCos ) fprintf(fp, ", ");
	  }
	  fprintf(fp, "}");
	  if( j!=NumOfMomHist-1 ) fprintf(fp, ", \n");
	}
	fprintf(fp, "}");
	if( i!=1 ) fprintf(fp, ", \n");
  }  
  fprintf(fp, "};\n");

  fprintf(fp, "double Efficiency_err[2][NumOfMomHist][NbinCos] = { \n");
  for( int i=0; i<2; i++ ){
	fprintf(fp, "{ ");
	for( int j=0; j<NumOfMomHist; j++ ){
	  fprintf(fp, "{ ");
	  for( int k=1; k<=NbinCos; k++ ){
		if( i==0 ) fprintf( fp, " %.3lf", hh5[j]->GetBinError(k) );
		if( i==1 ) fprintf( fp, " %.3lf", hh6[j]->GetBinError(k) );
		if( k!=NbinCos ) fprintf(fp, ", ");
	  }
	  fprintf(fp, "}");
	  if( j!=NumOfMomHist-1 ) fprintf(fp, ", \n");
	}
	fprintf(fp, "}");
	if( i!=1 ) fprintf(fp, ", \n");
  }  
  fprintf(fp, "};\n");
  
  fclose(fp);
 
  ///write rootfile////
  TFile *fout = new TFile("ForNPScatCutEfficiency.root", "recreate");
  for( int i=0; i<NumOfMomHist; i++ ){
	for( int j=0; j<NbinCos; j++ ){
	  h0[i][j]->Write();
	  h1[i][j]->Write();
	  h2[i][j]->Write();
	  h3[i][j]->Write();
	  h4[i][j]->Write();
	}			
	hh0[i]->Write();
	hh1[i]->Write();
	hh2[i]->Write();
	hh3[i]->Write();
	hh4[i]->Write();
	hh5[i]->Write();
	hh6[i]->Write();
  }
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
double getThetaCM( double m3, double m4, double cost, TLorentzVector PrimaryLv )
/*
  m3 = Scattered Particle mass (we know costLab and we want to know thetaCM);
  m4 = Other Particle  mass;
*/
{ 
  //TLorentzVector PrimaryLv = LvPi+LvP;
  double TotalEnergyCM = PrimaryLv.Mag();
  TVector3 beta( 1/PrimaryLv.E()*PrimaryLv.Vect() );

  //CM
  double TotalMomCM
	= 0.5*std::sqrt(( TotalEnergyCM*TotalEnergyCM - (m3+m4)*(m3+m4))
                    *( TotalEnergyCM*TotalEnergyCM - (m3-m4)*(m3-m4)))/TotalEnergyCM;

  double costLab = cost;
  double cottLab = costLab/std::sqrt(1.-costLab*costLab);
  double bt=beta.Mag(), gamma=1./std::sqrt(1.-bt*bt);
  double gbep=gamma*bt*std::sqrt(TotalMomCM*TotalMomCM+m3*m3)/TotalMomCM;
  double a  = gamma*gamma+cottLab*cottLab;
  double bp = gamma*gbep;
  double c  = gbep*gbep-cottLab*cottLab;
  double dd = bp*bp-a*c;
  
  if( dd<0. ){
    //std::cout<<"dd<0"<<std::endl;
    return -9999;
  } 
  
  double costCM = (std::sqrt(dd)-bp)/a;
  if( costCM>1. || costCM<-1. ){
    //std::cout<<"costCM>1. || costCM<-1."<<std::endl;
	return -9999;
  }

  double ThetaCM = acos(costCM)*Rad2Deg;

  /* double sintCM  = std::sqrt(1.-costCM*costCM);
  *pCalc = TotalMomCM*sintCM/std::sqrt(1.-costLab*costLab);
  */
  return ThetaCM;
}
