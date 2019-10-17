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
const int NumOfMom = 20;
const double MinMom = 0.45;
const double MaxMom = 0.85;
#include "LambdaNConvParam.txt"

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
  fChain->SetBranchStatus("nPi_CFT",1);
  fChain->SetBranchStatus("nP_CFT",1);
  fChain->SetBranchStatus("Vertex_x",1);//nSigma
  fChain->SetBranchStatus("Vertex_y",1);//nSigma
  fChain->SetBranchStatus("Vertex_z",1);//nSigma
  fChain->SetBranchStatus("momLambda",1);
  fChain->SetBranchStatus("momCalLambda",1);

  TH1D *h0[NumOfMom], *h1[NumOfMom], *h2[NumOfMom], *h3[NumOfMom], *h4[NumOfMom], *h5[NumOfMom];
  char buf[100], buf2[100];
  for (int i=0; i<NumOfMom; i++) {
	double low_range = MinMom+0.02*i;
	double up_range  = low_range+0.02;
	
	sprintf(buf, "h0_%d", i);
	sprintf(buf2, "Generated Num (All) (%.2lf < p < %.2lf)", low_range, up_range);
	h0[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	  
	sprintf(buf, "h1_%d", i);
	sprintf(buf2, "Generated Num (SigmaBeam requirement) (%.2lf < p < %.2lf)", low_range, up_range);
	h1[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "h2_%d", i);
	sprintf(buf2, "Detected Num (%.2lf < p < %.2lf)", low_range, up_range);
	h2[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "h3_%d", i);
	sprintf(buf2, "Detected Num (abs(dp)<10) (%.2lf < p < %.2lf)", low_range, up_range);
	h3[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	  
	sprintf(buf, "h4_%d", i);
	sprintf(buf2, "Efficiency (%.2lf < p < %.2lf)", low_range, up_range);
	h4[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "h5_%d", i);
	sprintf(buf2, "Efficiency (abs(dp)<10) (%.2lf < p < %.2lf)", low_range, up_range);
	h5[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
  }
  TH2D *H0 = new TH2D("H0","<Generate Num (All)> momentum vs costCM", NbinCos, -1, 1, NumOfMom, MinMom, MaxMom);
  TH2D *H1 = new TH2D("H1","<Generate Num (SigmaBeam requirement)> momentum vs costCM", NbinCos, -1, 1, NumOfMom, MinMom, MaxMom);
  TH2D *H2 = new TH2D("H2","<Detected Num> momentum vs costCM", NbinCos, -1, 1, NumOfMom, MinMom, MaxMom);
  TH2D *H3 = new TH2D("H3","<Detected Num (abs(dp)<10)> momentum vs costCM", NbinCos, -1, 1, NumOfMom, MinMom, MaxMom);
  TH2D *H4 = new TH2D("H4","<Efficiency> momentum vs costCM", NbinCos, -1, 1, NumOfMom, MinMom, MaxMom);
  TH2D *H5 = new TH2D("H5","<Efficiency (abs(dp)<10)> momentum vs costCM", NbinCos, -1, 1, NumOfMom, MinMom, MaxMom);
  
  int nev = fChain->GetEntries();
  std::cout<<"[Total Entry:"<<nev<<"]"<<std::endl;
  
  for(int evenum =0; evenum<nev; ++evenum ){
	//for(int evenum =0; evenum<20000; ++evenum ){
	if( evenum%1000000==0) std::cout<<"["<<(long int)100*evenum/nev<<"/100]"<<std::endl;
	
	scatFlag=-1;
	thetaScatHypCM=-9999;
	Vertex_x[0]=-9999;
	Vertex_y[0]=-9999;
	Vertex_z[0]=-9999;
	EkinCorP=-9999;
	scatEkinCal=-9999;
	momHypBeam=-9999;
	nPi_CFT=-1;
	nP_CFT=-1;

	fChain->GetEntry(evenum);	
	if( scatFlag!=1 ) continue;
	
	bool target_flag=false;
	TVector3 Vertex( Vertex_x[0], Vertex_y[0], Vertex_z[0] );
	if( getTargetFlag( Vertex ) ) target_flag=true;
      
	int index=-1;
	for ( int i=0; i<NumOfMom; i++ ) {
	  double low_range = MinMom+0.02*i;
	  double up_range  = low_range+0.02;
	  if ( low_range < momHypBeam&&momHypBeam < up_range ){
		index = i;
		break;
	  }
	}
	if( !( 0 <= index&&index < NumOfMom ) ) continue;
	
	double thetaCM = thetaScatHypCM;
	double dp = momLambda - momCalLambda;
	h0[index]->Fill( cos(thetaCM*Deg2Rad) );
	H0->Fill( cos(thetaCM*Deg2Rad) ,momHypBeam );
	if( target_flag ){
	  h1[index]->Fill( cos(thetaCM*Deg2Rad) );
	  H1->Fill( cos(thetaCM*Deg2Rad), momHypBeam );
	  if( nPi_CFT>0 && nP_CFT>0 && dp>Limit_Sim_Detectdp){
		h2[index]->Fill( cos(thetaCM*Deg2Rad) );
		H2->Fill( cos(thetaCM*Deg2Rad), momHypBeam );
		if( dp>-0.1 ){
		  h3[index]->Fill( cos(thetaCM*Deg2Rad) );
		  H3->Fill( cos(thetaCM*Deg2Rad), momHypBeam );
		}
	  }
	}    
	
  }//event loop
	
  for (int i=0; i<NumOfMom; i++) {
	for (int j=1; j<=NbinCos; j++) {
	  double val1 = h2[i]->GetBinContent(j);
	  double val2 = h1[i]->GetBinContent(j);
	  
	  if (val2 > 0) {
		h4[i]->SetBinContent(j, val1/val2);
		h4[i]->SetBinError(j, 1./val2*sqrt(val1*(val2-val1)/val2));
		H4->SetBinContent(j,i,val1/val2);
	  }
	  
	  val1 = h3[i]->GetBinContent(j);
	  val2 = h1[i]->GetBinContent(j);
      
	  if (val2 > 0) {
		h5[i]->SetBinContent(j, val1/val2);
		h5[i]->SetBinError(j, 1./val2*sqrt(val1*(val2-val1)/val2));
		H5->SetBinContent(j,i,val1/val2);
	  }
	}
  }
  
  
  FILE *fp = fopen("LambdaNConvAcceptance.txt","w");
  if(!fp) fprintf(stderr, "cannot open output file\n");
  
  fprintf(fp, "double Acceptance[NumOfMom][NbinCos] = { \n");
  for( int i=0; i<NumOfMom; i++ ){
	fprintf(fp, "{ ");
	for( int j=1; j<=NbinCos; j++ ){
	  fprintf( fp, " %.3lf", h4[i]->GetBinContent(j) );
	  if( j!=NbinCos ) fprintf(fp, ", ");
	}
	fprintf(fp, "}");
	if(i != NumOfMom-1) fprintf(fp, ", \n");
  }
  fprintf(fp, "};\n");
  
  fprintf(fp, "double Acceptance_err[NumOfMom][NbinCos] = { \n");
  for(int i=0; i<NumOfMom; i++) {
	fprintf(fp, "{ ");
	for(int j=1; j<=NbinCos; j++) {
	  fprintf(fp, " %.3lf", h4[i]->GetBinError(j));
	  if(j != NbinCos) fprintf(fp, ", ");
	  }
	fprintf(fp, "}");
	if(i != NumOfMom-1) fprintf(fp, ", \n");
  }
  fprintf(fp, "};\n");
  
  fclose(fp);
 
  ///write rootfile////
  TFile *fout = new TFile("ForLambdaNConvAcceptance.root", "recreate");
  for( int i=0; i<NumOfMom; i++ ){
	h0[i]->Write();
	h1[i]->Write();
	h2[i]->Write();
	h3[i]->Write();
	h4[i]->Write();
	h5[i]->Write();
  }			
  H0->Write();
  H1->Write();
  H2->Write();
  H3->Write();
  H4->Write();
  H5->Write();
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
