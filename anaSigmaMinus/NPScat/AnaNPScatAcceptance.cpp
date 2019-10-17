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
const int NumOfMom = 30;
const double MinMom = 0.35;
const double MaxMom = 0.65;
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
  fChain->SetBranchStatus("cdistDecayPi",1);
  fChain->SetBranchStatus("Vertex_z",1);//nSigma
  fChain->SetBranchStatus("vertexDecayPi",1);//3
  fChain->SetBranchStatus("EkinCorP",1);
  fChain->SetBranchStatus("scatNpEkinCal",1);

  TH1D *h0[NumOfMom], *h1[NumOfMom], *h2[NumOfMom], *h3[NumOfMom], *h4[NumOfMom], *h5[NumOfMom];
  char buf[100], buf2[100];
  for (int i=0; i<NumOfMom; i++) {
	double low_range = MinMom+0.01*i;
	double up_range  = low_range+0.01;
	
	sprintf(buf, "h0_%d", i);
	sprintf(buf2, "Generated Num (All) (%.2lf < p < %.2lf)", low_range, up_range);
	h0[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	  
	sprintf(buf, "h1_%d", i);
	sprintf(buf2, "Generated Num (NBeam requirement) (%.2lf < p < %.2lf)", low_range, up_range);
	h1[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "h2_%d", i);
	sprintf(buf2, "Detected Num (%.2lf < p < %.2lf)", low_range, up_range);
	h2[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "h3_%d", i);
	sprintf(buf2, "Detected Num (abs(dE)<10) (%.2lf < p < %.2lf)", low_range, up_range);
	h3[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	  
	sprintf(buf, "h4_%d", i);
	sprintf(buf2, "Efficiency (%.2lf < p < %.2lf)", low_range, up_range);
	h4[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
	
	sprintf(buf, "h5_%d", i);
	sprintf(buf2, "Efficiency (abs(dE)<10) (%.2lf < p < %.2lf)", low_range, up_range);
	h5[i] = new TH1D(buf, buf2, NbinCos, -1, 1);
  }
  TH2D *H0 = new TH2D("H0","<Generate Num (All)> momentum vs costCM", NbinCos, -1, 1, NumOfMom, MinMom, MaxMom);
  TH2D *H1 = new TH2D("H1","<Generate Num (NBeam requirement)> momentum vs costCM", NbinCos, -1, 1, NumOfMom, MinMom, MaxMom);
  TH2D *H2 = new TH2D("H2","<Detected Num> momentum vs costCM", NbinCos, -1, 1, NumOfMom, MinMom, MaxMom);
  TH2D *H3 = new TH2D("H3","<Detected Num (abs(dE)<10)> momentum vs costCM", NbinCos, -1, 1, NumOfMom, MinMom, MaxMom);
  TH2D *H4 = new TH2D("H4","<Efficiency> momentum vs costCM", NbinCos, -1, 1, NumOfMom, MinMom, MaxMom);
  TH2D *H5 = new TH2D("H5","<Efficiency (abs(dE)<10)> momentum vs costCM", NbinCos, -1, 1, NumOfMom, MinMom, MaxMom);
  
  int nev = fChain->GetEntries();
  std::cout<<"[Total Entry:"<<nev<<"]"<<std::endl;
  
  for(int evenum =0; evenum<nev; ++evenum ){
	//for(int evenum =0; evenum<20000; ++evenum ){
	if( evenum%1000000==0) std::cout<<"["<<(long int)100*evenum/nev<<"/100]"<<std::endl;
	
	for( int i=0; i<3; i++ ){
	  momVectorProtonScat[i]=-9999;
	  momVectorDecayNucleon[i]=-9999;
	  vertexDecayPi[i]=-9999;
	}
	NNscatFlag=-1;
	momProtonScat=-9999;
	momDecayNucleon=-9999;
	cdistDecayPi=-9999;
	Vertex_z[0]=-9999;
	EkinCorP=-9999;
	scatNpEkinCal=-9999;
	nPi_CFT=-1;
	nP_CFT=-1;

	fChain->GetEntry(evenum);	
	if( NNscatFlag<0 ) continue;
	
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
	for ( int i=0; i<NumOfMom; i++ ) {
	  double low_range = MinMom+0.01*i;
	  double up_range  = low_range+0.01;
	  if ( low_range < NBeamMom.Mag()&&NBeamMom.Mag() < up_range ){
		index = i;
		break;
	  }
	}
	if( !( 0 <= index&&index < NumOfMom ) ) continue;
	
	double thetaCM;
	double costLab = NBeamMom*ScatPMom/(NBeamMom.Mag()*ScatPMom.Mag());
	TLorentzVector LvN( NBeamMom, sqrt( NBeamMom.Mag()*NBeamMom.Mag()+NeutronMass*NeutronMass ));
	TLorentzVector LvP( 0, 0, 0, ProtonMass );
	TLorentzVector PrimaryLv = LvN + LvP;
	double thetaCM_P = getThetaCM( ProtonMass, NeutronMass, costLab, PrimaryLv );
	if( !( 0 <= thetaCM_P&&thetaCM_P <= 360) ) continue;
	else if( thetaCM_P+180<=360 ) thetaCM=thetaCM_P+180;
	else thetaCM = thetaCM_P-180;

	double dE = EkinCorP - scatNpEkinCal;
	h0[index]->Fill( cos(thetaCM*Deg2Rad) );
	H0->Fill( cos(thetaCM*Deg2Rad) ,NBeamMom.Mag() );
	if( nPi_CFT>0 && cdist_flag && dz_flag && target_flag ){
	  h1[index]->Fill( cos(thetaCM*Deg2Rad) );
	  H1->Fill( cos(thetaCM*Deg2Rad), NBeamMom.Mag() );
	  if( nPi_CFT>0 && nP_CFT>0 && abs(dE)<Limit_Sim_DetectdE ){
		h2[index]->Fill( cos(thetaCM*Deg2Rad) );
		H2->Fill( cos(thetaCM*Deg2Rad), NBeamMom.Mag() );
		if( fabs(dE)<10 ){
		  h3[index]->Fill( cos(thetaCM*Deg2Rad) );
		  H3->Fill( cos(thetaCM*Deg2Rad), NBeamMom.Mag() );
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
  
  
  FILE *fp = fopen("NPScatAcceptance.txt","w");
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
  TFile *fout = new TFile("ForNPScatAcceptance.root", "recreate");
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
