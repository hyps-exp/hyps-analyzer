//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Oct 12 02:06:28 2019 by ROOT version 6.16/00
// from TTree tree/tree of Sks
// found on file: run13955_ana.root
//////////////////////////////////////////////////////////

#ifndef G4Simulation_h
#define G4Simulation_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class G4Simulation {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        thetaMeson;
   Double_t        phiMeson;
   Double_t        thetaMesonCM;
   Double_t        phiMesonCM;
   Double_t        thetaScatHypCM;
   Double_t        phiScatHypCM;
   Double_t        thetaScatPLab;
   Double_t        thetaScatHypLab;
   Double_t        momVectorScatMeson[3];
   Double_t        momVectorHypBeam[3];
   Double_t        momVectorHypScat[3];
   Double_t        momVectorProtonScat[3];
   Double_t        momVectorDecayPi[3];
   Double_t        momVectorDecayNucleon[3];
   Double_t        momScatMeson;
   Double_t        momHypBeam;
   Double_t        momHypScat;
   Double_t        momProtonScat;
   Double_t        momDecayPi;
   Double_t        momDecayNucleon;
   Double_t        primaryVertex[3];
   Double_t        scatPos0[3];
   Double_t        NNscatPos[3];
   Double_t        PiNscatPos[3];
   Double_t        decayPos[3];
   Int_t           decayFlag;
   Int_t           scatFlag;
   Int_t           scatTarget;
   Int_t           NNscatFlag;
   Int_t           NNscatTarget;
   Int_t           PiNscatFlag;
   Int_t           PiNscatTarget;
   Int_t           ntSks;
   Double_t        p[4];   //[ntSks]
   Int_t           ntSksPart;
   Double_t        pSksPart[2];   //[ntSksPart]
   Double_t        m2[2];   //[ntSksPart]
   Int_t           nK;
   Double_t        pK[1];   //[nK]
   Double_t        MissMass[1];   //[nK]
   Double_t        u0SdcIn;
   Double_t        v0SdcIn;
   Double_t        resSftU;
   Double_t        resSftV;
   Int_t           nSigma;
   Double_t        pKCal[1];   //[nSigma]
   Double_t        theta[1];   //[nSigma]
   Double_t        pSigmaCal[1];   //[nSigma]
   Double_t        Vertex_x[1];   //[nSigma]
   Double_t        Vertex_y[1];   //[nSigma]
   Double_t        Vertex_z[1];   //[nSigma]
   Int_t           nPi_CFT;
   Int_t           nK_CFT;
   Int_t           nP_CFT;
   Double_t        EkinP;
   Double_t        EkinCorP;
   Double_t        vertexDecayPi[3];
   Double_t        cdistDecayPi;
   Double_t        thetaDecayPi;
   Double_t        decayPiMomCal;
   Double_t        decayPiMomVecCal;
   Double_t        decayNMomCal;
   Double_t        decayNMomVecCal[3];
   Double_t        vertexScat[3];
   Double_t        hypBeamVec[3];
   Double_t        cdistScat;
   Double_t        thetaScat;
   Double_t        scatMomCal;
   Double_t        scatEkinCal;
   Double_t        thetaScatCM;
   Double_t        thetaDecayPi2;
   Double_t        decayPiMomCal2;
   Double_t        vertexNpScat[3];
   Double_t        cdistNpScat;
   Double_t        thetaNpScat;
   Double_t        scatNpMomCal;
   Double_t        scatNpEkinCal;
   Double_t        thetaNpScatCM;
   Double_t        cdistLambdaDecay;
   Double_t        thetaLambdaDecay;
   Double_t        vertexLambdaDecay[3];
   Double_t        momPiFromLambda;
   Double_t        momProtonFromLambda;
   Double_t        momLambda;
   Double_t        momVecPiFromLambda;
   Double_t        momVecProtonFromLambda;
   Double_t        momVecLambda;
   Double_t        cdistLambdaNConv;
   Double_t        thetaLambdaNConv;
   Double_t        vertexLambdaNConv[3];
   Double_t        momCalLambda;
   Double_t        momCalLambda2;
   Double_t        thetaCMLambdaNConv;
   Double_t        vertexPiPScat[3];
   Double_t        cdistPiPScat;
   Double_t        piBeamMom;
   Double_t        thetaPiPScat;
   Double_t        scatPiMomCal;
   Double_t        scatPiMomCalFromE;
   Int_t           ntCFT;
   Double_t        xDirCFT[4];   //[ntCFT]
   Double_t        yDirCFT[4];   //[ntCFT]
   Double_t        zDirCFT[4];   //[ntCFT]
   Double_t        xPos0CFT[4];   //[ntCFT]
   Double_t        yPos0CFT[4];   //[ntCFT]
   Double_t        zPos0CFT[4];   //[ntCFT]
   Double_t        zBGOCFT[4];   //[ntCFT]
   Double_t        BGO_Edep[4];   //[ntCFT]
   Double_t        TotalEdep[4];   //[ntCFT]
   Double_t        CFT_TotalEdep[4];   //[ntCFT]
   Double_t        CFT_NormTotalEdep[4];   //[ntCFT]
   Double_t        PiV_Edep[4];   //[ntCFT]
   Double_t        xDirCFT_P[2];   //[nP_CFT]
   Double_t        yDirCFT_P[2];   //[nP_CFT]
   Double_t        zDirCFT_P[2];   //[nP_CFT]
   Double_t        xPos0CFT_P[2];   //[nP_CFT]
   Double_t        yPos0CFT_P[2];   //[nP_CFT]
   Double_t        zPos0CFT_P[2];   //[nP_CFT]
   Double_t        zBGOCFT_P[2];   //[nP_CFT]
   Double_t        CFTVtx_x_P[2];   //[nP_CFT]
   Double_t        CFTVtx_y_P[2];   //[nP_CFT]
   Double_t        CFTVtx_z_P[2];   //[nP_CFT]
   Double_t        BGO_Edep_P[2];   //[nP_CFT]
   Double_t        TotalEdep_P[2];   //[nP_CFT]
   Double_t        CFT_TotalEdep_P[2];   //[nP_CFT]
   Double_t        CFT_NormTotalEdep_P[2];   //[nP_CFT]
   Double_t        PiV_Edep_P[2];   //[nP_CFT]
   Double_t        xDirCFT_Pi[3];   //[nPi_CFT]
   Double_t        yDirCFT_Pi[3];   //[nPi_CFT]
   Double_t        zDirCFT_Pi[3];   //[nPi_CFT]
   Double_t        xPos0CFT_Pi[3];   //[nPi_CFT]
   Double_t        yPos0CFT_Pi[3];   //[nPi_CFT]
   Double_t        zPos0CFT_Pi[3];   //[nPi_CFT]
   Double_t        zBGOCFT_Pi[3];   //[nPi_CFT]
   Double_t        CFTVtx_x_Pi[3];   //[nPi_CFT]
   Double_t        CFTVtx_y_Pi[3];   //[nPi_CFT]
   Double_t        CFTVtx_z_Pi[3];   //[nPi_CFT]
   Double_t        u0CFT_Pi[3];   //[nPi_CFT]
   Double_t        v0CFT_Pi[3];   //[nPi_CFT]
   Double_t        BGO_Edep_Pi[3];   //[nPi_CFT]
   Double_t        TotalEdep_Pi[3];   //[nPi_CFT]
   Double_t        CFT_TotalEdep_Pi[3];   //[nPi_CFT]
   Double_t        CFT_NormTotalEdep_Pi[3];   //[nPi_CFT]
   Double_t        PiV_Edep_Pi[3];   //[nPi_CFT]

   // List of branches
   TBranch        *b_thetaMeson;   //!
   TBranch        *b_phiMeson;   //!
   TBranch        *b_thetaMesonCM;   //!
   TBranch        *b_phiMesonCM;   //!
   TBranch        *b_thetaScatHypCM;   //!
   TBranch        *b_phiScatHypCM;   //!
   TBranch        *b_thetaScatPLab;   //!
   TBranch        *b_thetaScatHypLab;   //!
   TBranch        *b_momVectorScatMeson;   //!
   TBranch        *b_momVectorHypBeam;   //!
   TBranch        *b_momVectorHypScat;   //!
   TBranch        *b_momVectorProtonScat;   //!
   TBranch        *b_momVectorDecayPi;   //!
   TBranch        *b_momVectorDecayNucleon;   //!
   TBranch        *b_momScatMeson;   //!
   TBranch        *b_momHypBeam;   //!
   TBranch        *b_momHypScat;   //!
   TBranch        *b_momProtonScat;   //!
   TBranch        *b_momDecayPi;   //!
   TBranch        *b_momDecayNucleon;   //!
   TBranch        *b_primaryVertex;   //!
   TBranch        *b_scatPos0;   //!
   TBranch        *b_NNscatPos;   //!
   TBranch        *b_PiNscatPos;   //!
   TBranch        *b_decayPos;   //!
   TBranch        *b_decayFlag;   //!
   TBranch        *b_scatFlag;   //!
   TBranch        *b_scatTarget;   //!
   TBranch        *b_NNscatFlag;   //!
   TBranch        *b_NNscatTarget;   //!
   TBranch        *b_PiNscatFlag;   //!
   TBranch        *b_PiNscatTarget;   //!
   TBranch        *b_ntSks;   //!
   TBranch        *b_p;   //!
   TBranch        *b_ntSksPart;   //!
   TBranch        *b_pSksPart;   //!
   TBranch        *b_m2;   //!
   TBranch        *b_nK;   //!
   TBranch        *b_pK;   //!
   TBranch        *b_MissMass;   //!
   TBranch        *b_u0SdcIn;   //!
   TBranch        *b_v0SdcIn;   //!
   TBranch        *b_resSftU;   //!
   TBranch        *b_resSftV;   //!
   TBranch        *b_nSigma;   //!
   TBranch        *b_pKCal;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_pSigmaCal;   //!
   TBranch        *b_Vertex_x;   //!
   TBranch        *b_Vertex_y;   //!
   TBranch        *b_Vertex_z;   //!
   TBranch        *b_nPi_CFT;   //!
   TBranch        *b_nK_CFT;   //!
   TBranch        *b_nP_CFT;   //!
   TBranch        *b_EkinP;   //!
   TBranch        *b_EkinCorP;   //!
   TBranch        *b_vertexDecayPi;   //!
   TBranch        *b_cdistDecayPi;   //!
   TBranch        *b_thetaDecayPi;   //!
   TBranch        *b_decayPiMomCal;   //!
   TBranch        *b_decayPiMomVecCal;   //!
   TBranch        *b_decayNMomCal;   //!
   TBranch        *b_decayNMomVecCal;   //!
   TBranch        *b_vertexScat;   //!
   TBranch        *b_hypBeamVec;   //!
   TBranch        *b_cdistScat;   //!
   TBranch        *b_thetaScat;   //!
   TBranch        *b_scatMomCal;   //!
   TBranch        *b_scatEkinCal;   //!
   TBranch        *b_thetaScatCM;   //!
   TBranch        *b_thetaDecayPi2;   //!
   TBranch        *b_decayPiMomCal2;   //!
   TBranch        *b_vertexNpScat;   //!
   TBranch        *b_cdistNpScat;   //!
   TBranch        *b_thetaNpScat;   //!
   TBranch        *b_scatNpMomCal;   //!
   TBranch        *b_scatNpEkinCal;   //!
   TBranch        *b_thetaNpScatCM;   //!
   TBranch        *b_cdistLambdaDecay;   //!
   TBranch        *b_thetaLambdaDecay;   //!
   TBranch        *b_vertexLambdaDecay;   //!
   TBranch        *b_momPiFromLambda;   //!
   TBranch        *b_momProtonFromLambda;   //!
   TBranch        *b_momLambda;   //!
   TBranch        *b_momVecPiFromLambda;   //!
   TBranch        *b_momVecProtonFromLambda;   //!
   TBranch        *b_momVecLambda;   //!
   TBranch        *b_cdistLambdaNConv;   //!
   TBranch        *b_thetaLambdaNConv;   //!
   TBranch        *b_vertexLambdaNConv;   //!
   TBranch        *b_momCalLambda;   //!
   TBranch        *b_momCalLambda2;   //!
   TBranch        *b_thetaCMLambdaNConv;   //!
   TBranch        *b_vertexPiPScat;   //!
   TBranch        *b_cdistPiPScat;   //!
   TBranch        *b_piBeamMom;   //!
   TBranch        *b_thetaPiPScat;   //!
   TBranch        *b_scatPiMomCal;   //!
   TBranch        *b_scatPiMomCalFromE;   //!
   TBranch        *b_ntCFT;   //!
   TBranch        *b_xDirCFT;   //!
   TBranch        *b_yDirCFT;   //!
   TBranch        *b_zDirCFT;   //!
   TBranch        *b_xPos0CFT;   //!
   TBranch        *b_yPos0CFT;   //!
   TBranch        *b_zPos0CFT;   //!
   TBranch        *b_zBGOCFT;   //!
   TBranch        *b_BGO_Edep;   //!
   TBranch        *b_TotalEdep;   //!
   TBranch        *b_CFT_TotalEdep;   //!
   TBranch        *b_CFT_NormTotalEdep;   //!
   TBranch        *b_PiV_Edep;   //!
   TBranch        *b_xDirCFT_P;   //!
   TBranch        *b_yDirCFT_P;   //!
   TBranch        *b_zDirCFT_P;   //!
   TBranch        *b_xPos0CFT_P;   //!
   TBranch        *b_yPos0CFT_P;   //!
   TBranch        *b_zPos0CFT_P;   //!
   TBranch        *b_zBGOCFT_P;   //!
   TBranch        *b_CFTVtx_x_P;   //!
   TBranch        *b_CFTVtx_y_P;   //!
   TBranch        *b_CFTVtx_z_P;   //!
   TBranch        *b_BGO_Edep_P;   //!
   TBranch        *b_TotalEdep_P;   //!
   TBranch        *b_CFT_TotalEdep_P;   //!
   TBranch        *b_CFT_NormTotalEdep_P;   //!
   TBranch        *b_PiV_Edep_P;   //!
   TBranch        *b_xDirCFT_Pi;   //!
   TBranch        *b_yDirCFT_Pi;   //!
   TBranch        *b_zDirCFT_Pi;   //!
   TBranch        *b_xPos0CFT_Pi;   //!
   TBranch        *b_yPos0CFT_Pi;   //!
   TBranch        *b_zPos0CFT_Pi;   //!
   TBranch        *b_zBGOCFT_Pi;   //!
   TBranch        *b_CFTVtx_x_Pi;   //!
   TBranch        *b_CFTVtx_y_Pi;   //!
   TBranch        *b_CFTVtx_z_Pi;   //!
   TBranch        *b_u0CFT_Pi;   //!
   TBranch        *b_v0CFT_Pi;   //!
   TBranch        *b_BGO_Edep_Pi;   //!
   TBranch        *b_TotalEdep_Pi;   //!
   TBranch        *b_CFT_TotalEdep_Pi;   //!
   TBranch        *b_CFT_NormTotalEdep_Pi;   //!
   TBranch        *b_PiV_Edep_Pi;   //!

   G4Simulation(TTree *tree=0);
   virtual ~G4Simulation();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef G4Simulation_cxx
G4Simulation::G4Simulation(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("run13955_ana.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("run13955_ana.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

G4Simulation::~G4Simulation()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t G4Simulation::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t G4Simulation::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void G4Simulation::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("thetaMeson", &thetaMeson, &b_thetaMeson);
   fChain->SetBranchAddress("phiMeson", &phiMeson, &b_phiMeson);
   fChain->SetBranchAddress("thetaMesonCM", &thetaMesonCM, &b_thetaMesonCM);
   fChain->SetBranchAddress("phiMesonCM", &phiMesonCM, &b_phiMesonCM);
   fChain->SetBranchAddress("thetaScatHypCM", &thetaScatHypCM, &b_thetaScatHypCM);
   fChain->SetBranchAddress("phiScatHypCM", &phiScatHypCM, &b_phiScatHypCM);
   fChain->SetBranchAddress("thetaScatPLab", &thetaScatPLab, &b_thetaScatPLab);
   fChain->SetBranchAddress("thetaScatHypLab", &thetaScatHypLab, &b_thetaScatHypLab);
   fChain->SetBranchAddress("momVectorScatMeson", momVectorScatMeson, &b_momVectorScatMeson);
   fChain->SetBranchAddress("momVectorHypBeam", momVectorHypBeam, &b_momVectorHypBeam);
   fChain->SetBranchAddress("momVectorHypScat", momVectorHypScat, &b_momVectorHypScat);
   fChain->SetBranchAddress("momVectorProtonScat", momVectorProtonScat, &b_momVectorProtonScat);
   fChain->SetBranchAddress("momVectorDecayPi", momVectorDecayPi, &b_momVectorDecayPi);
   fChain->SetBranchAddress("momVectorDecayNucleon", momVectorDecayNucleon, &b_momVectorDecayNucleon);
   fChain->SetBranchAddress("momScatMeson", &momScatMeson, &b_momScatMeson);
   fChain->SetBranchAddress("momHypBeam", &momHypBeam, &b_momHypBeam);
   fChain->SetBranchAddress("momHypScat", &momHypScat, &b_momHypScat);
   fChain->SetBranchAddress("momProtonScat", &momProtonScat, &b_momProtonScat);
   fChain->SetBranchAddress("momDecayPi", &momDecayPi, &b_momDecayPi);
   fChain->SetBranchAddress("momDecayNucleon", &momDecayNucleon, &b_momDecayNucleon);
   fChain->SetBranchAddress("primaryVertex", primaryVertex, &b_primaryVertex);
   fChain->SetBranchAddress("scatPos0", scatPos0, &b_scatPos0);
   fChain->SetBranchAddress("NNscatPos", NNscatPos, &b_NNscatPos);
   fChain->SetBranchAddress("PiNscatPos", PiNscatPos, &b_PiNscatPos);
   fChain->SetBranchAddress("decayPos", decayPos, &b_decayPos);
   fChain->SetBranchAddress("decayFlag", &decayFlag, &b_decayFlag);
   fChain->SetBranchAddress("scatFlag", &scatFlag, &b_scatFlag);
   fChain->SetBranchAddress("scatTarget", &scatTarget, &b_scatTarget);
   fChain->SetBranchAddress("NNscatFlag", &NNscatFlag, &b_NNscatFlag);
   fChain->SetBranchAddress("NNscatTarget", &NNscatTarget, &b_NNscatTarget);
   fChain->SetBranchAddress("PiNscatFlag", &PiNscatFlag, &b_PiNscatFlag);
   fChain->SetBranchAddress("PiNscatTarget", &PiNscatTarget, &b_PiNscatTarget);
   fChain->SetBranchAddress("ntSks", &ntSks, &b_ntSks);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("ntSksPart", &ntSksPart, &b_ntSksPart);
   fChain->SetBranchAddress("pSksPart", pSksPart, &b_pSksPart);
   fChain->SetBranchAddress("m2", m2, &b_m2);
   fChain->SetBranchAddress("nK", &nK, &b_nK);
   fChain->SetBranchAddress("pK", pK, &b_pK);
   fChain->SetBranchAddress("MissMass", MissMass, &b_MissMass);
   fChain->SetBranchAddress("u0SdcIn", &u0SdcIn, &b_u0SdcIn);
   fChain->SetBranchAddress("v0SdcIn", &v0SdcIn, &b_v0SdcIn);
   fChain->SetBranchAddress("resSftU", &resSftU, &b_resSftU);
   fChain->SetBranchAddress("resSftV", &resSftV, &b_resSftV);
   fChain->SetBranchAddress("nSigma", &nSigma, &b_nSigma);
   fChain->SetBranchAddress("pKCal", pKCal, &b_pKCal);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("pSigmaCal", pSigmaCal, &b_pSigmaCal);
   fChain->SetBranchAddress("Vertex_x", Vertex_x, &b_Vertex_x);
   fChain->SetBranchAddress("Vertex_y", Vertex_y, &b_Vertex_y);
   fChain->SetBranchAddress("Vertex_z", Vertex_z, &b_Vertex_z);
   fChain->SetBranchAddress("nPi_CFT", &nPi_CFT, &b_nPi_CFT);
   fChain->SetBranchAddress("nK_CFT", &nK_CFT, &b_nK_CFT);
   fChain->SetBranchAddress("nP_CFT", &nP_CFT, &b_nP_CFT);
   fChain->SetBranchAddress("EkinP", &EkinP, &b_EkinP);
   fChain->SetBranchAddress("EkinCorP", &EkinCorP, &b_EkinCorP);
   fChain->SetBranchAddress("vertexDecayPi", vertexDecayPi, &b_vertexDecayPi);
   fChain->SetBranchAddress("cdistDecayPi", &cdistDecayPi, &b_cdistDecayPi);
   fChain->SetBranchAddress("thetaDecayPi", &thetaDecayPi, &b_thetaDecayPi);
   fChain->SetBranchAddress("decayPiMomCal", &decayPiMomCal, &b_decayPiMomCal);
   fChain->SetBranchAddress("decayPiMomVecCal", &decayPiMomVecCal, &b_decayPiMomVecCal);
   fChain->SetBranchAddress("decayNMomCal", &decayNMomCal, &b_decayNMomCal);
   fChain->SetBranchAddress("decayNMomVecCal", decayNMomVecCal, &b_decayNMomVecCal);
   fChain->SetBranchAddress("vertexScat", vertexScat, &b_vertexScat);
   fChain->SetBranchAddress("hypBeamVec", hypBeamVec, &b_hypBeamVec);
   fChain->SetBranchAddress("cdistScat", &cdistScat, &b_cdistScat);
   fChain->SetBranchAddress("thetaScat", &thetaScat, &b_thetaScat);
   fChain->SetBranchAddress("scatMomCal", &scatMomCal, &b_scatMomCal);
   fChain->SetBranchAddress("scatEkinCal", &scatEkinCal, &b_scatEkinCal);
   fChain->SetBranchAddress("thetaScatCM", &thetaScatCM, &b_thetaScatCM);
   fChain->SetBranchAddress("thetaDecayPi2", &thetaDecayPi2, &b_thetaDecayPi2);
   fChain->SetBranchAddress("decayPiMomCal2", &decayPiMomCal2, &b_decayPiMomCal2);
   fChain->SetBranchAddress("vertexNpScat", vertexNpScat, &b_vertexNpScat);
   fChain->SetBranchAddress("cdistNpScat", &cdistNpScat, &b_cdistNpScat);
   fChain->SetBranchAddress("thetaNpScat", &thetaNpScat, &b_thetaNpScat);
   fChain->SetBranchAddress("scatNpMomCal", &scatNpMomCal, &b_scatNpMomCal);
   fChain->SetBranchAddress("scatNpEkinCal", &scatNpEkinCal, &b_scatNpEkinCal);
   fChain->SetBranchAddress("thetaNpScatCM", &thetaNpScatCM, &b_thetaNpScatCM);
   fChain->SetBranchAddress("cdistLambdaDecay", &cdistLambdaDecay, &b_cdistLambdaDecay);
   fChain->SetBranchAddress("thetaLambdaDecay", &thetaLambdaDecay, &b_thetaLambdaDecay);
   fChain->SetBranchAddress("vertexLambdaDecay", vertexLambdaDecay, &b_vertexLambdaDecay);
   fChain->SetBranchAddress("momPiFromLambda", &momPiFromLambda, &b_momPiFromLambda);
   fChain->SetBranchAddress("momProtonFromLambda", &momProtonFromLambda, &b_momProtonFromLambda);
   fChain->SetBranchAddress("momLambda", &momLambda, &b_momLambda);
   fChain->SetBranchAddress("momVecPiFromLambda", &momVecPiFromLambda, &b_momVecPiFromLambda);
   fChain->SetBranchAddress("momVecProtonFromLambda", &momVecProtonFromLambda, &b_momVecProtonFromLambda);
   fChain->SetBranchAddress("momVecLambda", &momVecLambda, &b_momVecLambda);
   fChain->SetBranchAddress("cdistLambdaNConv", &cdistLambdaNConv, &b_cdistLambdaNConv);
   fChain->SetBranchAddress("thetaLambdaNConv", &thetaLambdaNConv, &b_thetaLambdaNConv);
   fChain->SetBranchAddress("vertexLambdaNConv", vertexLambdaNConv, &b_vertexLambdaNConv);
   fChain->SetBranchAddress("momCalLambda", &momCalLambda, &b_momCalLambda);
   fChain->SetBranchAddress("momCalLambda2", &momCalLambda2, &b_momCalLambda2);
   fChain->SetBranchAddress("thetaCMLambdaNConv", &thetaCMLambdaNConv, &b_thetaCMLambdaNConv);
   fChain->SetBranchAddress("vertexPiPScat", vertexPiPScat, &b_vertexPiPScat);
   fChain->SetBranchAddress("cdistPiPScat", &cdistPiPScat, &b_cdistPiPScat);
   fChain->SetBranchAddress("piBeamMom", &piBeamMom, &b_piBeamMom);
   fChain->SetBranchAddress("thetaPiPScat", &thetaPiPScat, &b_thetaPiPScat);
   fChain->SetBranchAddress("scatPiMomCal", &scatPiMomCal, &b_scatPiMomCal);
   fChain->SetBranchAddress("scatPiMomCalFromE", &scatPiMomCalFromE, &b_scatPiMomCalFromE);
   fChain->SetBranchAddress("ntCFT", &ntCFT, &b_ntCFT);
   fChain->SetBranchAddress("xDirCFT", xDirCFT, &b_xDirCFT);
   fChain->SetBranchAddress("yDirCFT", yDirCFT, &b_yDirCFT);
   fChain->SetBranchAddress("zDirCFT", zDirCFT, &b_zDirCFT);
   fChain->SetBranchAddress("xPos0CFT", xPos0CFT, &b_xPos0CFT);
   fChain->SetBranchAddress("yPos0CFT", yPos0CFT, &b_yPos0CFT);
   fChain->SetBranchAddress("zPos0CFT", zPos0CFT, &b_zPos0CFT);
   fChain->SetBranchAddress("zBGOCFT", zBGOCFT, &b_zBGOCFT);
   fChain->SetBranchAddress("BGO_Edep", BGO_Edep, &b_BGO_Edep);
   fChain->SetBranchAddress("TotalEdep", TotalEdep, &b_TotalEdep);
   fChain->SetBranchAddress("CFT_TotalEdep", CFT_TotalEdep, &b_CFT_TotalEdep);
   fChain->SetBranchAddress("CFT_NormTotalEdep", CFT_NormTotalEdep, &b_CFT_NormTotalEdep);
   fChain->SetBranchAddress("PiV_Edep", PiV_Edep, &b_PiV_Edep);
   fChain->SetBranchAddress("xDirCFT_P", xDirCFT_P, &b_xDirCFT_P);
   fChain->SetBranchAddress("yDirCFT_P", yDirCFT_P, &b_yDirCFT_P);
   fChain->SetBranchAddress("zDirCFT_P", zDirCFT_P, &b_zDirCFT_P);
   fChain->SetBranchAddress("xPos0CFT_P", xPos0CFT_P, &b_xPos0CFT_P);
   fChain->SetBranchAddress("yPos0CFT_P", yPos0CFT_P, &b_yPos0CFT_P);
   fChain->SetBranchAddress("zPos0CFT_P", zPos0CFT_P, &b_zPos0CFT_P);
   fChain->SetBranchAddress("zBGOCFT_P", zBGOCFT_P, &b_zBGOCFT_P);
   fChain->SetBranchAddress("CFTVtx_x_P", CFTVtx_x_P, &b_CFTVtx_x_P);
   fChain->SetBranchAddress("CFTVtx_y_P", CFTVtx_y_P, &b_CFTVtx_y_P);
   fChain->SetBranchAddress("CFTVtx_z_P", CFTVtx_z_P, &b_CFTVtx_z_P);
   fChain->SetBranchAddress("BGO_Edep_P", BGO_Edep_P, &b_BGO_Edep_P);
   fChain->SetBranchAddress("TotalEdep_P", TotalEdep_P, &b_TotalEdep_P);
   fChain->SetBranchAddress("CFT_TotalEdep_P", CFT_TotalEdep_P, &b_CFT_TotalEdep_P);
   fChain->SetBranchAddress("CFT_NormTotalEdep_P", CFT_NormTotalEdep_P, &b_CFT_NormTotalEdep_P);
   fChain->SetBranchAddress("PiV_Edep_P", PiV_Edep_P, &b_PiV_Edep_P);
   fChain->SetBranchAddress("xDirCFT_Pi", xDirCFT_Pi, &b_xDirCFT_Pi);
   fChain->SetBranchAddress("yDirCFT_Pi", yDirCFT_Pi, &b_yDirCFT_Pi);
   fChain->SetBranchAddress("zDirCFT_Pi", zDirCFT_Pi, &b_zDirCFT_Pi);
   fChain->SetBranchAddress("xPos0CFT_Pi", xPos0CFT_Pi, &b_xPos0CFT_Pi);
   fChain->SetBranchAddress("yPos0CFT_Pi", yPos0CFT_Pi, &b_yPos0CFT_Pi);
   fChain->SetBranchAddress("zPos0CFT_Pi", zPos0CFT_Pi, &b_zPos0CFT_Pi);
   fChain->SetBranchAddress("zBGOCFT_Pi", zBGOCFT_Pi, &b_zBGOCFT_Pi);
   fChain->SetBranchAddress("CFTVtx_x_Pi", CFTVtx_x_Pi, &b_CFTVtx_x_Pi);
   fChain->SetBranchAddress("CFTVtx_y_Pi", CFTVtx_y_Pi, &b_CFTVtx_y_Pi);
   fChain->SetBranchAddress("CFTVtx_z_Pi", CFTVtx_z_Pi, &b_CFTVtx_z_Pi);
   fChain->SetBranchAddress("u0CFT_Pi", u0CFT_Pi, &b_u0CFT_Pi);
   fChain->SetBranchAddress("v0CFT_Pi", v0CFT_Pi, &b_v0CFT_Pi);
   fChain->SetBranchAddress("BGO_Edep_Pi", BGO_Edep_Pi, &b_BGO_Edep_Pi);
   fChain->SetBranchAddress("TotalEdep_Pi", TotalEdep_Pi, &b_TotalEdep_Pi);
   fChain->SetBranchAddress("CFT_TotalEdep_Pi", CFT_TotalEdep_Pi, &b_CFT_TotalEdep_Pi);
   fChain->SetBranchAddress("CFT_NormTotalEdep_Pi", CFT_NormTotalEdep_Pi, &b_CFT_NormTotalEdep_Pi);
   fChain->SetBranchAddress("PiV_Edep_Pi", PiV_Edep_Pi, &b_PiV_Edep_Pi);
   Notify();
}

Bool_t G4Simulation::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void G4Simulation::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t G4Simulation::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef G4Simulation_cxx
