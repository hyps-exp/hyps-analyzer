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

double SigmaFlightLength( TVector3 vertex, TVector3 pSigma );
bool decayCheck( double ctau, double momentum, double mass, double dx);
double calcEnergyDeposit( double momentum, double mass, double dx );
double calc_dE_dx( double beta ); 
bool scatteringCheck( double rate, double dx );
double MakeRandom( double start, double end );
bool getTargetFlag( TVector3 pos );
const double pi = 3.14159265358;

const int NumOfMom = 20;
const int NumOfMomHist = 3;
const double MinMom[NumOfMomHist] = { 0.45, 0.55, 0.65 };
const double MaxMom[NumOfMomHist] = { 0.55, 0.65, 0.85 };

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
  fChain->SetBranchStatus("*",0);  // disable all branches
  fChain->SetBranchStatus("nPiKCatch",1);  // activate branchname
  fChain->SetBranchStatus("MissMassCorrDE",1);  // activate branchname
  fChain->SetBranchStatus("chisqrKuramapik",1);  // activate branchname
  fChain->SetBranchStatus("pKuramapik",1);  // activate branchname
  fChain->SetBranchStatus("qKuramapik",1);  // activate branchname
  fChain->SetBranchStatus("ccm2",1);  // activate branchname
  fChain->SetBranchStatus("KURAMAPID",1);  // activate branchname
  fChain->SetBranchStatus("vtx_KURAMA",1);  // activate branchname
  fChain->SetBranchStatus("vty_KURAMA",1);  // activate branchname
  fChain->SetBranchStatus("vtz_KURAMA",1);  // activate branchname
  fChain->SetBranchStatus("MissMomcal",1);  // activate branchname
  fChain->SetBranchStatus("MissMomxcal",1);  // activate branchname
  fChain->SetBranchStatus("MissMomycal",1);  // activate branchname
  fChain->SetBranchStatus("MissMomzcal",1);  // activate branchname

  const int MAXnPiKCatch = 24;

  int nev = fChain->GetEntries();
  std::cout<<"[Total Entry:"<<nev<<"]"<<std::endl;

  TVector3 v;//mm
  TVector3 p;//MeV

  double FlightLen[NumOfMom];
  for( int i=0; i<NumOfMom; i++ ) FlightLen[i]=0;

  for(int evenum =0; evenum<nev; ++evenum ){
	//for(int evenum =0; evenum<50000; ++evenum ){
	if( evenum%20000==0) std::cout<<"["<<(long int)100*evenum/nev<<"/100]"<<std::endl;
	
    nPiKCatch=-9999.;
    for( int i=0; i<MAXnPiKCatch; i++ ){
	  chisqrKuramapik[i]=-9999.;
	  pKuramapik[i]=-9999.;
	  qKuramapik[i]=-9999.;
	  MissMassCorrDE[i]=-9999.;
	  ccm2[i]=-9999.;
	  MissMomcal[i]=-9999.;
	  MissMomxcal[i]=-9999.;
	  MissMomycal[i]=-9999.;
	  MissMomzcal[i]=-9999.;
	  vtx_KURAMA[i]=-9999.;
	  vty_KURAMA[i]=-9999.;
	  vtz_KURAMA[i]=-9999.;
    }
	bool flag=false;
    	
    fChain->GetEntry(evenum);
    for( int iPiK=0; iPiK<nPiKCatch; iPiK++ ){
	  if( chisqrKuramapik[iPiK]<50 && pKuramapik[iPiK]<0.89 && qKuramapik[iPiK]>0 && 0.15<ccm2[iPiK]&&ccm2[iPiK]<0.35 && 1.16<MissMassCorrDE[iPiK]&&MissMassCorrDE[iPiK]<1.25 ){
		if( flag ) continue;
		else flag=true;
		v.SetXYZ(vtx_KURAMA[iPiK],vty_KURAMA[iPiK],vtz_KURAMA[iPiK]);
		p.SetXYZ(1000*MissMomxcal[iPiK],1000*MissMomycal[iPiK],1000*MissMomzcal[iPiK]);
		if( getTargetFlag(v) ){
		  double l = SigmaFlightLength( v,p );
		  for ( int i=0; i<NumOfMom; i++ ) {
			double low_range = MinMom[0]+0.02*i;
			double up_range  = low_range+0.02;
			if ( low_range < MissMomcal[iPiK]&&MissMomcal[iPiK] < up_range ) {
			  FlightLen[i]+=l/10.;
			  break;
			}
		  }
		  
		}
	  }//sigma event
	}
  }//event loop

  FILE *fp = fopen("SigmaBeam.txt","w");
  if(!fp) fprintf(stderr, "cannot open output file\n");
  
  fprintf(fp, "double FlightLen[NumOfMom] = { \n");
  for( int i=0; i<NumOfMom; i++ ){
	  fprintf( fp, " %.3lf", FlightLen[i] );
	  if(i != NumOfMom-1) fprintf(fp, ", \n");
  }
  fprintf(fp, "};\n");
    
  fclose(fp);

}
/////////////////////////////////////////////////////////////////////////
double SigmaFlightLength( TVector3 vertex, TVector3 pSigma ){

  double ctau=44.34; /*mm*/

  double p_sigma = pSigma.Mag();//MeV
  double m_sigma = 1197.4; //MeV
 
  double dx = 0.1; //mm
  double totalx=0.0;  //mm

  TVector3 localSigmaPos = vertex;
  bool flagDecay=false;
  bool flagScattering=false;
  bool flagTgt=true;

  int nIteration=0;
  while (1) {
    totalx += dx;
    localSigmaPos.SetX( localSigmaPos.X() + pSigma.X()*dx/pSigma.Mag() );
    localSigmaPos.SetY( localSigmaPos.Y() + pSigma.Y()*dx/pSigma.Mag() );
    localSigmaPos.SetZ( localSigmaPos.Z() + pSigma.Z()*dx/pSigma.Mag() );

    p_sigma = calcEnergyDeposit( p_sigma, m_sigma, dx );
    flagDecay = decayCheck(ctau, p_sigma, m_sigma, dx );
    if (fabs(p_sigma)<0.000001){
      std::cout<<"stop in mateerial"<<std::endl;
      flagDecay = true;
    }
    if (flagDecay){
      //std::cout<<"decay"<<std::endl;//cout
      break;
    }
    flagTgt = getTargetFlag(localSigmaPos);
    if (!flagTgt){
      //std::cout<<"out of target"<<std::endl;//cout
      break;
    }
    nIteration++;
    if (nIteration>5000) break;
  }
  return totalx;
}
/////////////////////////////////////////////////////////////////////////
bool decayCheck( double ctau, double momentum, double mass, double dx){

  double lambda=ctau*momentum/mass;
  double value = exp(-dx/lambda);

  double x = MakeRandom(0,1);
  if (x > value) return true;
  else return false;  
}
/////////////////////////////////////////////////////////////////////////
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
double calcEnergyDeposit( double momentum, double mass, double dx ){

  double E0 = sqrt(momentum*momentum+mass*mass); // MeV
  double beta = momentum/sqrt(momentum*momentum+mass*mass);
  
  double dEdx = calc_dE_dx( beta );// MeV/cm
  if (dEdx < 0.) return -1;
  
  double dE = dEdx * dx / 10.; // dx (mm), MeV
  double E = E0 - dE;          // MeV

  if (E-mass<0.005)  return 0.0; //particle stop in terget
  double p = sqrt(E*E-mass*mass); //MeV/c

  return p;
}
///////////////////////////////////////////////////////////////////////////
double calc_dE_dx( double beta ){
  double value;
  const double C=0.1535; /*MeVcm^2/g*/
  const double m_e=0.511;
  double logterm;
  double gamma_2;
  double W_max;
  int z = 1;

  double rho = 0.0708; /* g/cm^3 */
  double I = 21.8;   /* eV */
  double Z_A = 0.99212;

  gamma_2=1/(1-pow(beta,2.0));

  W_max=2.0*m_e*pow(beta,2.0)*gamma_2;
  logterm=log(2.0*m_e*gamma_2*pow(beta,2.0)*W_max*pow(10.0,12.0)/pow(I,2.0))-2.0*pow(beta,2.0);

  value=C*rho*Z_A*pow((double)z,2.0)*logterm/pow(beta,2.0);

  return value;
}
/////////////////////////////////////////////////////////////////////////
bool scatteringCheck( double rate, double dx ){
  double value = rate*dx;
  double x = MakeRandom(0,1);
  if (x < value) return true;
  else return false;  
}
/////////////////////////////////////////////////////////////////////////
double MakeRandom( double start, double end ){
  std::random_device rnd;  //this is seed value
  std::mt19937 mt(rnd());  // this is random nnumber generator
  std::uniform_real_distribution<> rand( start, end ); 
  
  return rand(mt);
}
