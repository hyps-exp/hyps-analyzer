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
double GetFlightLength( TVector3 v1, TVector3 v2, TVector3 v3 );
const double Deg2Rad = 3.14159265358/180;
const double Rad2Deg = 180/3.14159265358;

const int NbinCos = 20;
const int NumOfMomHist = 3;
const double MinMom[NumOfMomHist] = { 0.45, 0.55, 0.65 };
const double MaxMom[NumOfMomHist] = { 0.55, 0.65, 0.85 };
const int NbindE = 100;
#include "SigmaPScatParam.txt"

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
  fChain->SetBranchStatus("scatFlag",1);
  fChain->SetBranchStatus("thetaScatHypCM",1);
  fChain->SetBranchStatus("momHypBeam",1);
  fChain->SetBranchStatus("cdistScat",1);
  fChain->SetBranchStatus("cdistDecayPi",1);
  fChain->SetBranchStatus("nPi_CFT",1);
  fChain->SetBranchStatus("nP_CFT",1);
  fChain->SetBranchStatus("nSigma",1);
  fChain->SetBranchStatus("Vertex_x",1);//nSigma
  fChain->SetBranchStatus("Vertex_y",1);//nSigma
  fChain->SetBranchStatus("Vertex_z",1);//nSigma
  fChain->SetBranchStatus("vertexScat",1);//3
  fChain->SetBranchStatus("vertexDecayPi",1);//3
  fChain->SetBranchStatus("EkinCorP",1);
  fChain->SetBranchStatus("scatEkinCal",1);
  fChain->SetBranchStatus("scatNpEkinCal",1);
  fChain->SetBranchStatus("momLambda",1);
  fChain->SetBranchStatus("momCalLambda",1);

  TH1D *h0[NumOfMomHist][NbinCos], *h1[NumOfMomHist][NbinCos], *h2[NumOfMomHist][NbinCos], *h3[NumOfMomHist][NbinCos], *h4[NumOfMomHist][NbinCos],
	*h5[NumOfMomHist][NbinCos], *h6[NumOfMomHist][NbinCos], *h7[NumOfMomHist][NbinCos], *h8[NumOfMomHist][NbinCos], *h9[NumOfMomHist][NbinCos];
  TH1D *hh0[NumOfMomHist], *hh1[NumOfMomHist], *hh2[NumOfMomHist], *hh3[NumOfMomHist], *hh4[NumOfMomHist], *hh5[NumOfMomHist],
	*hh6[NumOfMomHist], *hh7[NumOfMomHist], *hh8[NumOfMomHist], *hh9[NumOfMomHist], *hh10[NumOfMomHist], *hh11[NumOfMomHist];

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
	  sprintf(buf2, "dE (Detected && cdist1 requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	  h2[i][j] = new TH1D(buf, buf2, NbindE, -100, 100);
	  
	  sprintf(buf, "h3_%d_%d", i,j);
	  sprintf(buf2, "dE (Detected && dz1 requirement)(%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	  h3[i][j] = new TH1D(buf, buf2, NbindE, -100, 100);
	  
	  sprintf(buf, "h4_%d_%d", i,j);
	  sprintf(buf2, "dE (Detected && cdist2 requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	  h4[i][j] = new TH1D(buf, buf2, NbindE, -100, 100);

	  sprintf(buf, "h5_%d_%d", i,j);
	  sprintf(buf2, "dE (Detected && dz2 requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	  h5[i][j] = new TH1D(buf, buf2, NbindE, -100, 100);

	  sprintf(buf, "h6_%d_%d", i,j);
	  sprintf(buf2, "dE (Detected && flength requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	  h6[i][j] = new TH1D(buf, buf2, NbindE, -100, 100);

	  sprintf(buf, "h7_%d_%d", i,j);
	  sprintf(buf2, "dE (Detected && NPCut requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	  h7[i][j] = new TH1D(buf, buf2, NbindE, -100, 100);

	  sprintf(buf, "h8_%d_%d", i,j);
	  sprintf(buf2, "dE (Detected && LambdaNCut requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	  h8[i][j] = new TH1D(buf, buf2, NbindE, -100, 100);

	  sprintf(buf, "h9_%d_%d", i,j);
	  sprintf(buf2, "dE (Detected && AllCut requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	  h9[i][j] = new TH1D(buf, buf2, NbindE, -100, 100);
	}
	sprintf(buf, "hh0_%d", i);
	sprintf(buf2, "costCM (All) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh0[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "hh1_%d", i);
	sprintf(buf2, "costCM (Detected requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh1[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "hh2_%d", i);
	sprintf(buf2, "costCM (Detected && cdist1 requirement)(%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh2[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "hh3_%d", i);
	sprintf(buf2, "costCM (Detected && dz1 requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh3[i] = new TH1D(buf, buf2, NbinCos, -1, 1);

	sprintf(buf, "hh4_%d", i);
	sprintf(buf2, "costCM (Detected && cdist2 requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh4[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "hh5_%d", i);
	sprintf(buf2, "costCM (Detected && dz2 requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh5[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "hh6_%d", i);
	sprintf(buf2, "costCM (Detected && flength requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh6[i] = new TH1D(buf, buf2, NbinCos, -1, 1);

	sprintf(buf, "hh7_%d", i);
	sprintf(buf2, "costCM (Detected && NPCut requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh7[i] = new TH1D(buf, buf2, NbinCos, -1, 1);

	sprintf(buf, "hh8_%d", i);
	sprintf(buf2, "costCM (Detected && LambdaNCut requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh8[i] = new TH1D(buf, buf2, NbinCos, -1, 1);

	sprintf(buf, "hh9_%d", i);
	sprintf(buf2, "costCM (Detected && AllCut requirement) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh9[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "hh10_%d", i);
	sprintf(buf2, "costCM (Efficiency1) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh10[i] = new TH1D(buf, buf2, NbinCos, -1, 1);

	sprintf(buf, "hh11_%d", i);
	sprintf(buf2, "costCM (Efficiency2) (%.2lf < p < %.2lf)", MinMom[i], MaxMom[i]);
	hh11[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
  }
  
  int nev = fChain->GetEntries();
  std::cout<<"[Total Entry:"<<nev<<"]"<<std::endl;
  
  for(int evenum =0; evenum<nev; ++evenum ){
	//for(int evenum =0; evenum<20000; ++evenum ){
	if( evenum%1000000==0) std::cout<<"["<<(long int)100*evenum/nev<<"/100]"<<std::endl;
	
	for( int i=0; i<3; i++ ){
	  vertexScat[i]=-9999;
	  vertexDecayPi[i]=-9999;
	}
	scatFlag=-1;
	thetaScatHypCM=-9999;
	cdistScat=-9999;
	cdistDecayPi=-9999;
	Vertex_x[0]=-9999;
	Vertex_y[0]=-9999;
	Vertex_z[0]=-9999;
	EkinCorP=-9999;
	scatEkinCal=-9999;
	scatNpEkinCal=-9999;
	momHypBeam=-9999;
	momLambda=-9999;
	momCalLambda=-9999;
	nPi_CFT=-1;
	nP_CFT=-1;
	nSigma=-1;

	fChain->GetEntry(evenum);	
	if( scatFlag!=1 ) continue;
	if( nSigma!=1 ) continue;

	TVector3 Vertex( Vertex_x[0], Vertex_y[0], Vertex_z[0] );
	TVector3 VertexScat( vertexScat[0], vertexScat[1], vertexScat[2] );
	TVector3 VertexDecay( vertexDecayPi[0], vertexDecayPi[1], vertexDecayPi[2] );

	bool target_flag=false;
	if( getTargetFlag( Vertex ) ) target_flag=true;
	bool cdist1_flag=false;
	if( cdistScat<Limit_Sim_cdistScat ) cdist1_flag=true;
	bool dz1_flag=false;
	if( vertexScat[2]-Vertex_z[0]>Limit_Sim_dz1 ) dz1_flag=true;
	bool cdist2_flag=false;
	if( cdistDecayPi<Limit_Sim_cdistDecayPi ) cdist2_flag=true;
	bool dz2_flag=false;
	if( vertexDecayPi[2]-vertexScat[2]>Limit_Sim_dz2 ) dz2_flag=true;
	bool flength_flag=false;
	if( GetFlightLength( Vertex, VertexScat, VertexDecay )<Limit_Sim_flength ) flength_flag=true;
	bool NPCut_flag=false;
	if( abs(EkinCorP-scatNpEkinCal)>Limit_Sim_NPCut ) NPCut_flag=true;
	bool LambdaNCut_flag=false;
	if( momLambda-momCalLambda<Limit_Sim_LambdaNCut ) LambdaNCut_flag=true;

	int index=-1;
	for ( int i=0; i<NumOfMomHist; i++ ) {
	  if ( MinMom[i] < momHypBeam&&momHypBeam < MaxMom[i] ){
		index = i;
		break;
	  }
	}
	if( !( 0 <= index&&index < NumOfMomHist ) ) continue;

	double thetaCM = thetaScatHypCM;
	int cosbin=hh0[0]->FindBin(cos(thetaCM*Deg2Rad));
	
	double dE = EkinCorP - scatEkinCal;
	h0[index][cosbin-1]->Fill( dE );
	hh0[index]->Fill( cos(thetaCM*Deg2Rad) );
	if( target_flag && nPi_CFT>0 && nP_CFT>0 && abs(dE)<Limit_Sim_DetectdE ){//Detected
	  h1[index][cosbin-1]->Fill( dE );
	  hh1[index]->Fill( cos(thetaCM*Deg2Rad) );
	  if( cdist1_flag ){
	  h2[index][cosbin-1]->Fill( dE );
	  hh2[index]->Fill( cos(thetaCM*Deg2Rad) );		
	  }
	  if( dz1_flag ){
		h3[index][cosbin-1]->Fill( dE );
		hh3[index]->Fill( cos(thetaCM*Deg2Rad) );
	  }
	  if( cdist2_flag ){
		h4[index][cosbin-1]->Fill( dE );
		hh4[index]->Fill( cos(thetaCM*Deg2Rad) );	
   	  }
	  if( dz2_flag ){
		h5[index][cosbin-1]->Fill( dE );
		hh5[index]->Fill( cos(thetaCM*Deg2Rad) );
	  }
	  if( flength_flag ){
		h6[index][cosbin-1]->Fill( dE );
		hh6[index]->Fill( cos(thetaCM*Deg2Rad) );
	  }    
	  if( NPCut_flag ){
		h7[index][cosbin-1]->Fill( dE );
		hh7[index]->Fill( cos(thetaCM*Deg2Rad) );
	  }
	  if( LambdaNCut_flag ){
		h8[index][cosbin-1]->Fill( dE );
		hh8[index]->Fill( cos(thetaCM*Deg2Rad) );
	  }
	  if( cdist1_flag && dz1_flag && cdist2_flag && dz2_flag && flength_flag && NPCut_flag /*&& LambdaNCut*/ ){
		h9[index][cosbin-1]->Fill( dE );
		hh9[index]->Fill( cos(thetaCM*Deg2Rad) );
	  }
	}
	
  }//event loop
	
  for(int i=0; i<NumOfMomHist; i++){
	for (int j=1; j<=NbinCos; j++){
	  double val1 = h9[i][j]->Integral(h9[i][j]->FindBin(-100),h9[i][j]->FindBin(100));
	  double val2 = h1[i][j]->Integral(h1[i][j]->FindBin(-100),h1[i][j]->FindBin(100));
	  if(val2 > 0){
		hh10[i]->SetBinContent(j, val1/val2);
		hh10[i]->SetBinError(j, 1./val2*sqrt(val1*(val2-val1)/val2));
	  }

	  val1 = h9[i][j]->Integral(h9[i][j]->FindBin(-Limit_Sim_CutdE),h9[i][j]->FindBin(Limit_Sim_CutdE));
	  if(val2 > 0){
		hh11[i]->SetBinContent(j, val1/val2);
		hh11[i]->SetBinError(j, 1./val2*sqrt(val1*(val2-val1)/val2));
	  }
	}
  }
    
  FILE *fp = fopen("SigmaPScatCutEfficiency.txt","w");
  if(!fp) fprintf(stderr, "cannot open output file\n");
  
  fprintf(fp, "double Efficiency[2][NumOfMomHist][NbinCos] = { \n");
  for( int i=0; i<2; i++ ){
	fprintf(fp, "{ ");
	for( int j=0; j<NumOfMomHist; j++ ){
	  fprintf(fp, "{ ");
	  for( int k=1; k<=NbinCos; k++ ){
		if( i==0 ) fprintf( fp, " %.3lf", hh10[j]->GetBinContent(k) );
		if( i==1 ) fprintf( fp, " %.3lf", hh11[j]->GetBinContent(k) );
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
		if( i==0 ) fprintf( fp, " %.3lf", hh10[j]->GetBinError(k) );
		if( i==1 ) fprintf( fp, " %.3lf", hh11[j]->GetBinError(k) );
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
  TFile *fout = new TFile("ForSigmaPScatCutEfficiency.root", "recreate");
  for( int i=0; i<NumOfMomHist; i++ ){
	for( int j=0; j<NbinCos; j++ ){
	  h0[i][j]->Write();
	  h1[i][j]->Write();
	  h2[i][j]->Write();
	  h3[i][j]->Write();
	  h4[i][j]->Write();
	  h5[i][j]->Write();
	  h6[i][j]->Write();
	  h7[i][j]->Write();
	  h8[i][j]->Write();
	  h9[i][j]->Write();
	}			
	hh0[i]->Write();
	hh1[i]->Write();
	hh2[i]->Write();
	hh3[i]->Write();
	hh4[i]->Write();
	hh5[i]->Write();
	hh6[i]->Write();
	hh7[i]->Write();
	hh8[i]->Write();
	hh9[i]->Write();
	hh10[i]->Write();
	hh11[i]->Write();
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
//////////////////////////////////////////////////////////////////////////
double GetFlightLength( TVector3 v1, TVector3 v2, TVector3 v3 ){
  TVector3 VertexLength1 = v2-v1;
  TVector3 VertexLength2 = v3-v2;
  double flightlength = VertexLength1.Mag()+VertexLength2.Mag();
  return flightlength;
}
