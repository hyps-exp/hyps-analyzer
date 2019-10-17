//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct  8 23:33:26 2019 by ROOT version 6.16/00
// from TTree pik/tree of PiKAna
// found on file: ../dst_rootfile/run06444_DstPiKAna_SM_Sigma.root
//////////////////////////////////////////////////////////

#ifndef DstPiKAna_SM_h
#define DstPiKAna_SM_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DstPiKAna_SM {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runnum;
   Int_t           evnum;
   Int_t           spill;
   Int_t           trignhits;
   Int_t           trigpat[17];   //[trignhits]
   Int_t           trigflag[32];
   Int_t           nhBh1;
   Int_t           csBh1[5];   //[nhBh1]
   Double_t        Bh1Seg[5];   //[nhBh1]
   Double_t        tBh1[5];   //[nhBh1]
   Double_t        dtBh1[5];   //[nhBh1]
   Double_t        deBh1[5];   //[nhBh1]
   Double_t        btof[5];   //[nhBh1]
   Int_t           nhBh2;
   Int_t           csBh2[3];   //[nhBh2]
   Double_t        Bh2Seg[3];   //[nhBh2]
   Double_t        tBh2[3];   //[nhBh2]
   Double_t        t0Bh2[3];   //[nhBh2]
   Double_t        dtBh2[3];   //[nhBh2]
   Double_t        deBh2[3];   //[nhBh2]
   Int_t           nhTof;
   Int_t           csTof[10];   //[nhTof]
   Double_t        TofSeg[10];   //[nhTof]
   Double_t        tTof[10];   //[nhTof]
   Double_t        dtTof[10];   //[nhTof]
   Double_t        deTof[10];   //[nhTof]
   Int_t           bft_ncl;
   Int_t           nlBcOut;
   Int_t           ntBcOut;
   Int_t           nhBcOut[5];   //[ntBcOut]
   Double_t        chisqrBcOut[5];   //[ntBcOut]
   Double_t        x0BcOut[5];   //[ntBcOut]
   Double_t        y0BcOut[5];   //[ntBcOut]
   Double_t        u0BcOut[5];   //[ntBcOut]
   Double_t        v0BcOut[5];   //[ntBcOut]
   Double_t        xtgtBcOut[5];   //[ntBcOut]
   Double_t        ytgtBcOut[5];   //[ntBcOut]
   Double_t        xbh2BcOut[5];   //[ntBcOut]
   Double_t        ybh2BcOut[5];   //[ntBcOut]
   Int_t           ntK18;
   Int_t           nhK18[5];   //[ntK18]
   Double_t        chisqrK18[5];   //[ntK18]
   Double_t        pK18[5];   //[ntK18]
   Double_t        xtgtK18[5];   //[ntK18]
   Double_t        ytgtK18[5];   //[ntK18]
   Double_t        utgtK18[5];   //[ntK18]
   Double_t        vtgtK18[5];   //[ntK18]
   Double_t        thetaK18[5];   //[ntK18]
   Int_t           nlSdcIn;
   Int_t           ntSdcIn;
   Int_t           nhSdcIn[4];   //[ntSdcIn]
   Double_t        chisqrSdcIn[4];   //[ntSdcIn]
   Double_t        x0SdcIn[4];   //[ntSdcIn]
   Double_t        y0SdcIn[4];   //[ntSdcIn]
   Double_t        u0SdcIn[4];   //[ntSdcIn]
   Double_t        v0SdcIn[4];   //[ntSdcIn]
   Int_t           nlSdcOut;
   Int_t           ntSdcOut;
   Int_t           nhSdcOut[4];   //[ntSdcOut]
   Double_t        chisqrSdcOut[4];   //[ntSdcOut]
   Double_t        x0SdcOut[4];   //[ntSdcOut]
   Double_t        y0SdcOut[4];   //[ntSdcOut]
   Double_t        u0SdcOut[4];   //[ntSdcOut]
   Double_t        v0SdcOut[4];   //[ntSdcOut]
   Int_t           ntKurama;
   Int_t           nhKurama[4];   //[ntKurama]
   Double_t        stof[4];   //[ntKurama]
   Double_t        cstof[4];   //[ntKurama]
   Double_t        path[4];   //[ntKurama]
   Double_t        chisqrKurama[4];   //[ntKurama]
   Double_t        pKurama[4];   //[ntKurama]
   Double_t        qKurama[4];   //[ntKurama]
   Double_t        m2[4];   //[ntKurama]
   Double_t        cm2[4];   //[ntKurama]
   Double_t        xtgtKurama[4];   //[ntKurama]
   Double_t        ytgtKurama[4];   //[ntKurama]
   Double_t        utgtKurama[4];   //[ntKurama]
   Double_t        vtgtKurama[4];   //[ntKurama]
   Double_t        thetaKurama[4];   //[ntKurama]
   Double_t        xtofKurama[4];   //[ntKurama]
   Double_t        ytofKurama[4];   //[ntKurama]
   Double_t        utofKurama[4];   //[ntKurama]
   Double_t        vtofKurama[4];   //[ntKurama]
   Double_t        tofsegKurama[4];   //[ntKurama]
   Double_t        best_deTof[4];   //[ntKurama]
   Double_t        best_TofSeg[4];   //[ntKurama]
   Int_t           ntCFT;
   Int_t           ntProton;
   Int_t           ntOther;
   Double_t        theta_cft[5];   //[ntCFT]
   Double_t        phi_cft[5];   //[ntCFT]
   Double_t        vtx_cft[5];   //[ntCFT]
   Double_t        vty_cft[5];   //[ntCFT]
   Double_t        vtz_cft[5];   //[ntCFT]
   Double_t        Total_dE[5];   //[ntCFT]
   Double_t        Total_dEphi[5];   //[ntCFT]
   Double_t        Total_dEuv[5];   //[ntCFT]
   Double_t        Total_dE_max[5];   //[ntCFT]
   Double_t        Total_dEphi_max[5];   //[ntCFT]
   Double_t        Total_dEuv_max[5];   //[ntCFT]
   Int_t           nhit_phi[5];   //[ntCFT]
   Int_t           nhit_uv[5];   //[ntCFT]
   Int_t           segBGOt[5];   //[ntCFT]
   Double_t        energyBGOt[5];   //[ntCFT]
   Int_t           segPiIDt[5];   //[ntCFT]
   Int_t           protonflagt[5];   //[ntCFT]
   Int_t           BGOnohitt[5];   //[ntCFT]
   Int_t           simKuramat[5];   //[ntCFT]
   Int_t           tracknoproton[3];   //[ntProton]
   Int_t           tracknoother[5];   //[ntOther]
   Int_t           nPi;
   Int_t           nK;
   Int_t           nPiK;
   Int_t           nCatch;
   Int_t           nPiKCatch;
   Double_t        vtx_KURAMA[36];   //[nPiKCatch]
   Double_t        vty_KURAMA[36];   //[nPiKCatch]
   Double_t        vtz_KURAMA[36];   //[nPiKCatch]
   Double_t        ccm2[36];   //[nPiKCatch]
   Double_t        cpath[36];   //[nPiKCatch]
   Double_t        chisqrKuramapik[36];   //[nPiKCatch]
   Double_t        pKuramapik[36];   //[nPiKCatch]
   Double_t        qKuramapik[36];   //[nPiKCatch]
   Int_t           simKuramatpik[36];   //[nPiKCatch]
   Int_t           Kflag[36];   //[nPiKCatch]
   Double_t        closedist_KURAMA[36];   //[nPiKCatch]
   Double_t        theta[36];   //[nPiKCatch]
   Double_t        MissMass[36];   //[nPiKCatch]
   Double_t        MissMassCorr[36];   //[nPiKCatch]
   Double_t        MissMassCorrDE[36];   //[nPiKCatch]
   Double_t        MissMass2CorrDE[36];   //[nPiKCatch]
   Double_t        MissMassPiPi[36];   //[nPiKCatch]
   Double_t        MissMassPiPiCorrDE[36];   //[nPiKCatch]
   Double_t        MissMass2PiPiCorrDE[36];   //[nPiKCatch]
   Double_t        pCorrDE_PiPi[36];   //[nPiKCatch]
   Double_t        pCalcCorrDE_PiPi[36];   //[nPiKCatch]
   Double_t        MissMassPiP[36];   //[nPiKCatch]
   Double_t        MissMassPiPCorrDE[36];   //[nPiKCatch]
   Double_t        MissMass2PiPCorrDE[36];   //[nPiKCatch]
   Double_t        pCorrDE_PiP[36];   //[nPiKCatch]
   Double_t        pCalcCorrDE_PiP[36];   //[nPiKCatch]
   Double_t        thetaCM[36];   //[nPiKCatch]
   Double_t        thetaCMCorrDE[36];   //[nPiKCatch]
   Double_t        thetaCMCorrDE_PiPi[36];   //[nPiKCatch]
   Double_t        thetaCMCorrDE_PiP[36];   //[nPiKCatch]
   Double_t        costCM[36];   //[nPiKCatch]
   Int_t           KURAMAPID[36];   //[nPiKCatch]
   Double_t        MissMom[36];   //[nPiKCatch]
   Double_t        MissMomx[36];   //[nPiKCatch]
   Double_t        MissMomy[36];   //[nPiKCatch]
   Double_t        MissMomz[36];   //[nPiKCatch]
   Double_t        MissMomcal[36];   //[nPiKCatch]
   Double_t        MissMomxcal[36];   //[nPiKCatch]
   Double_t        MissMomycal[36];   //[nPiKCatch]
   Double_t        MissMomzcal[36];   //[nPiKCatch]
   Double_t        SigmaBeta[36];   //[nPiKCatch]
   Double_t        xpi[36];   //[nPiKCatch]
   Double_t        ypi[36];   //[nPiKCatch]
   Double_t        upi[36];   //[nPiKCatch]
   Double_t        vpi[36];   //[nPiKCatch]
   Double_t        xk[36];   //[nPiKCatch]
   Double_t        yk[36];   //[nPiKCatch]
   Double_t        uk[36];   //[nPiKCatch]
   Double_t        vk[36];   //[nPiKCatch]
   Double_t        uc[36];   //[nPiKCatch]
   Double_t        vc[36];   //[nPiKCatch]
   Double_t        pOrg[36];   //[nPiKCatch]
   Double_t        pCalc[36];   //[nPiKCatch]
   Double_t        pCalcCorrDE[36];   //[nPiKCatch]
   Double_t        pCorr[36];   //[nPiKCatch]
   Double_t        pCorrDE[36];   //[nPiKCatch]
   Double_t        pCorrDE_PiBeam[36];   //[nPiKCatch]
   Double_t        vtx_K18Catch[36];   //[nPiKCatch]
   Double_t        vty_K18Catch[36];   //[nPiKCatch]
   Double_t        vtz_K18Catch[36];   //[nPiKCatch]
   Double_t        closedist_K18Catch[36];   //[nPiKCatch]
   Double_t        theta_K18Catch[36];   //[nPiKCatch]
   Double_t        vertex_distance[36];   //[nPiKCatch]
   Double_t        vtx_XCatch[36];   //[nPiKCatch]
   Double_t        vty_XCatch[36];   //[nPiKCatch]
   Double_t        vtz_XCatch[36];   //[nPiKCatch]
   Double_t        closedist_XCatch[36];   //[nPiKCatch]
   Double_t        theta_XCatch[36];   //[nPiKCatch]
   Double_t        vtx_XCatchCorrDE[36];   //[nPiKCatch]
   Double_t        vty_XCatchCorrDE[36];   //[nPiKCatch]
   Double_t        vtz_XCatchCorrDE[36];   //[nPiKCatch]
   Double_t        closedist_XCatchCorrDE[36];   //[nPiKCatch]
   Double_t        theta_XCatchCorrDE[36];   //[nPiKCatch]
   Double_t        vtx_WCatch[36];   //[nPiKCatch]
   Double_t        vty_WCatch[36];   //[nPiKCatch]
   Double_t        vtz_WCatch[36];   //[nPiKCatch]
   Double_t        closedist_WCatch[36];   //[nPiKCatch]
   Double_t        PiBeamPScatPMom;
   Double_t        PiBeamPScatPiMom;
   Double_t        theta_PiBeamPScatPi;
   Int_t           PiBeamPScatNo;
   Double_t        NBeamMom;
   Double_t        NBeamMomx;
   Double_t        NBeamMomy;
   Double_t        NBeamMomz;
   Double_t        NBeamLength;
   Double_t        vtx_NBeam;
   Double_t        vty_NBeam;
   Double_t        vtz_NBeam;
   Double_t        vdist_NBeam;
   Double_t        cdistNBeam;
   Int_t           pitrno_NBeam;
   Int_t           pikno_NBeam;
   Double_t        DeltaE_NPScat;
   Double_t        ProtonMom_NPScat;
   Double_t        DecayNeutronMom;
   Double_t        DecayPionMom;
   Double_t        MissMassSigmaP_NPScat;
   Double_t        vtx_Decay2np;
   Double_t        vty_Decay2np;
   Double_t        vtz_Decay2np;
   Double_t        vtx_NPScat;
   Double_t        vty_NPScat;
   Double_t        vtz_NPScat;
   Double_t        vdist1_NPScat;
   Double_t        vdist2_NPScat;
   Double_t        cdistDecay2np;
   Double_t        cdistNPScat;
   Double_t        theta_NPScat;
   Double_t        thetaCM_NPScat;
   Int_t           ptrno_NPScat;
   Int_t           pitrno_NPScat;
   Int_t           pikno_NPScat;
   Int_t           priority_NPScat;
   Double_t        DeltaE_SigmaPScat;
   Double_t        MissMassSigmaP_SigmaPScat;
   Double_t        vtx_SigmaPScat;
   Double_t        vty_SigmaPScat;
   Double_t        vtz_SigmaPScat;
   Double_t        vdist_SigmaPScat;
   Double_t        cdistSigmaPScat;
   Int_t           ptrno_SigmaPScat;
   Int_t           pikno_SigmaPScat;
   Double_t        DeltaE_SigmaPScat2npi;
   Double_t        ProtonMom_SigmaPScat2npi;
   Double_t        MissMassSigmaP_SigmaPScat2npi;
   Double_t        vtx_ScatSigmaDecay2npi;
   Double_t        vty_ScatSigmaDecay2npi;
   Double_t        vtz_ScatSigmaDecay2npi;
   Double_t        cdistScatSigmaDecay2npi;
   Double_t        vdist_SigmaPScat2npi;
   Double_t        theta_SigmaPScat2npi;
   Double_t        thetaCM_SigmaPScat2npi;
   Int_t           ptrno_SigmaPScat2npi;
   Int_t           pitrno_SigmaPScat2npi;
   Int_t           pikno_SigmaPScat2npi;
   Int_t           priority_SigmaPScat2npi;
   Double_t        DeltaP_LambdaNConv;
   Double_t        NeutronMom_LambdaNConv;
   Double_t        LambdaMom;
   Double_t        thetaCM_LamdaNConv;
   Double_t        vtx_LambdaNConv;
   Double_t        vty_LambdaNConv;
   Double_t        vtz_LambdaNConv;
   Double_t        vtx_LambdaDecay;
   Double_t        vty_LambdaDecay;
   Double_t        vtz_LambdaDecay;
   Double_t        vdist1_LambdaNConv;
   Double_t        vdist2_LambdaNConv;
   Double_t        cdistLambdaNConv;
   Double_t        cdistLambdaDecay;
   Int_t           ptrno_LambdaNConv;
   Int_t           pitrno_LambdaNConv;
   Int_t           pikno_LambdaNConv;
   Int_t           priority_LambdaNConv;
   Double_t        DeltaP_PiPScat;
   Double_t        DecayNeutronMom2;
   Double_t        DecayPionMom2;
   Double_t        vtx_Decay2pip;
   Double_t        vty_Decay2pip;
   Double_t        vtz_Decay2pip;
   Double_t        vtx_PiPScat;
   Double_t        vty_PiPScat;
   Double_t        vtz_PiPScat;
   Double_t        vdist1_PiPScat;
   Double_t        vdist2_PiPScat;
   Double_t        cdistDecay2pip;
   Double_t        cdistPiPScat;
   Int_t           ptrno_PiPScat;
   Int_t           pitrno_PiPScat;
   Int_t           pikno_PiPScat;
   Int_t           priority_PiPScat;

   // List of branches
   TBranch        *b_runnum;   //!
   TBranch        *b_evnum;   //!
   TBranch        *b_spill;   //!
   TBranch        *b_trignhits;   //!
   TBranch        *b_trigpat;   //!
   TBranch        *b_trigflag;   //!
   TBranch        *b_nhBh1;   //!
   TBranch        *b_csBh1;   //!
   TBranch        *b_Bh1Seg;   //!
   TBranch        *b_tBh1;   //!
   TBranch        *b_dtBh1;   //!
   TBranch        *b_deBh1;   //!
   TBranch        *b_btof;   //!
   TBranch        *b_nhBh2;   //!
   TBranch        *b_csBh2;   //!
   TBranch        *b_Bh2Seg;   //!
   TBranch        *b_tBh2;   //!
   TBranch        *b_t0Bh2;   //!
   TBranch        *b_dtBh2;   //!
   TBranch        *b_deBh2;   //!
   TBranch        *b_nhTof;   //!
   TBranch        *b_csTof;   //!
   TBranch        *b_TofSeg;   //!
   TBranch        *b_tTof;   //!
   TBranch        *b_dtTof;   //!
   TBranch        *b_deTof;   //!
   TBranch        *b_bft_ncl;   //!
   TBranch        *b_nlBcOut;   //!
   TBranch        *b_ntBcOut;   //!
   TBranch        *b_nhBcOut;   //!
   TBranch        *b_chisqrBcOut;   //!
   TBranch        *b_x0BcOut;   //!
   TBranch        *b_y0BcOut;   //!
   TBranch        *b_u0BcOut;   //!
   TBranch        *b_v0BcOut;   //!
   TBranch        *b_xtgtBcOut;   //!
   TBranch        *b_ytgtBcOut;   //!
   TBranch        *b_xbh2BcOut;   //!
   TBranch        *b_ybh2BcOut;   //!
   TBranch        *b_ntK18;   //!
   TBranch        *b_nhK18;   //!
   TBranch        *b_chisqrK18;   //!
   TBranch        *b_pK18;   //!
   TBranch        *b_xtgtK18;   //!
   TBranch        *b_ytgtK18;   //!
   TBranch        *b_utgtK18;   //!
   TBranch        *b_vtgtK18;   //!
   TBranch        *b_thetaK18;   //!
   TBranch        *b_nlSdcIn;   //!
   TBranch        *b_ntSdcIn;   //!
   TBranch        *b_nhSdcIn;   //!
   TBranch        *b_chisqrSdcIn;   //!
   TBranch        *b_x0SdcIn;   //!
   TBranch        *b_y0SdcIn;   //!
   TBranch        *b_u0SdcIn;   //!
   TBranch        *b_v0SdcIn;   //!
   TBranch        *b_nlSdcOut;   //!
   TBranch        *b_ntSdcOut;   //!
   TBranch        *b_nhSdcOut;   //!
   TBranch        *b_chisqrSdcOut;   //!
   TBranch        *b_x0SdcOut;   //!
   TBranch        *b_y0SdcOut;   //!
   TBranch        *b_u0SdcOut;   //!
   TBranch        *b_v0SdcOut;   //!
   TBranch        *b_ntKurama;   //!
   TBranch        *b_nhKurama;   //!
   TBranch        *b_stof;   //!
   TBranch        *b_cstof;   //!
   TBranch        *b_path;   //!
   TBranch        *b_chisqrKurama;   //!
   TBranch        *b_pKurama;   //!
   TBranch        *b_qKurama;   //!
   TBranch        *b_m2;   //!
   TBranch        *b_cm2;   //!
   TBranch        *b_xtgtKurama;   //!
   TBranch        *b_ytgtKurama;   //!
   TBranch        *b_utgtKurama;   //!
   TBranch        *b_vtgtKurama;   //!
   TBranch        *b_thetaKurama;   //!
   TBranch        *b_xtofKurama;   //!
   TBranch        *b_ytofKurama;   //!
   TBranch        *b_utofKurama;   //!
   TBranch        *b_vtofKurama;   //!
   TBranch        *b_tofsegKurama;   //!
   TBranch        *b_best_deTof;   //!
   TBranch        *b_best_TofSeg;   //!
   TBranch        *b_ntCFT;   //!
   TBranch        *b_ntProton;   //!
   TBranch        *b_ntOther;   //!
   TBranch        *b_theta_cft;   //!
   TBranch        *b_phi_cft;   //!
   TBranch        *b_vtx_cft;   //!
   TBranch        *b_vty_cft;   //!
   TBranch        *b_vtz_cft;   //!
   TBranch        *b_Total_dE;   //!
   TBranch        *b_Total_dEphi;   //!
   TBranch        *b_Total_dEuv;   //!
   TBranch        *b_Total_dE_max;   //!
   TBranch        *b_Total_dEphi_max;   //!
   TBranch        *b_Total_dEuv_max;   //!
   TBranch        *b_nhit_phi;   //!
   TBranch        *b_nhit_uv;   //!
   TBranch        *b_segBGOt;   //!
   TBranch        *b_energyBGOt;   //!
   TBranch        *b_segPiIDt;   //!
   TBranch        *b_protonflagt;   //!
   TBranch        *b_BGOnohitt;   //!
   TBranch        *b_simKuramat;   //!
   TBranch        *b_tracknoproton;   //!
   TBranch        *b_tracknoother;   //!
   TBranch        *b_nPi;   //!
   TBranch        *b_nK;   //!
   TBranch        *b_nPiK;   //!
   TBranch        *b_nCatch;   //!
   TBranch        *b_nPiKCatch;   //!
   TBranch        *b_vtx_KURAMA;   //!
   TBranch        *b_vty_KURAMA;   //!
   TBranch        *b_vtz_KURAMA;   //!
   TBranch        *b_ccm2;   //!
   TBranch        *b_cpath;   //!
   TBranch        *b_chisqrKuramapik;   //!
   TBranch        *b_pKuramapik;   //!
   TBranch        *b_qKuramapik;   //!
   TBranch        *b_simKuramatpik;   //!
   TBranch        *b_Kflag;   //!
   TBranch        *b_closedist_KURAMA;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_MissMass;   //!
   TBranch        *b_MissMassCorr;   //!
   TBranch        *b_MissMassCorrDE;   //!
   TBranch        *b_MissMass2CorrDE;   //!
   TBranch        *b_MissMassPiPi;   //!
   TBranch        *b_MissMassPiPiCorrDE;   //!
   TBranch        *b_MissMass2PiPiCorrDE;   //!
   TBranch        *b_pCorrDE_PiPi;   //!
   TBranch        *b_pCalcCorrDE_PiPi;   //!
   TBranch        *b_MissMassPiP;   //!
   TBranch        *b_MissMassPiPCorrDE;   //!
   TBranch        *b_MissMass2PiPCorrDE;   //!
   TBranch        *b_pCorrDE_PiP;   //!
   TBranch        *b_pCalcCorrDE_PiP;   //!
   TBranch        *b_thetaCM;   //!
   TBranch        *b_thetaCMCorrDE;   //!
   TBranch        *b_thetaCMCorrDE_PiPi;   //!
   TBranch        *b_thetaCMCorrDE_PiP;   //!
   TBranch        *b_costCM;   //!
   TBranch        *b_KURAMAPID;   //!
   TBranch        *b_MissMom;   //!
   TBranch        *b_MissMomx;   //!
   TBranch        *b_MissMomy;   //!
   TBranch        *b_MissMomz;   //!
   TBranch        *b_MissMomcal;   //!
   TBranch        *b_MissMomxcal;   //!
   TBranch        *b_MissMomycal;   //!
   TBranch        *b_MissMomzcal;   //!
   TBranch        *b_SigmaBeta;   //!
   TBranch        *b_xpi;   //!
   TBranch        *b_ypi;   //!
   TBranch        *b_upi;   //!
   TBranch        *b_vpi;   //!
   TBranch        *b_xk;   //!
   TBranch        *b_yk;   //!
   TBranch        *b_uk;   //!
   TBranch        *b_vk;   //!
   TBranch        *b_uc;   //!
   TBranch        *b_vc;   //!
   TBranch        *b_pOrg;   //!
   TBranch        *b_pCalc;   //!
   TBranch        *b_pCalcCorrDE;   //!
   TBranch        *b_pCorr;   //!
   TBranch        *b_pCorrDE;   //!
   TBranch        *b_pCorrDE_PiBeam;   //!
   TBranch        *b_vtx_K18Catch;   //!
   TBranch        *b_vty_K18Catch;   //!
   TBranch        *b_vtz_K18Catch;   //!
   TBranch        *b_closedist_K18Catch;   //!
   TBranch        *b_theta_K18Catch;   //!
   TBranch        *b_vertex_distance;   //!
   TBranch        *b_vtx_XCatch;   //!
   TBranch        *b_vty_XCatch;   //!
   TBranch        *b_vtz_XCatch;   //!
   TBranch        *b_closedist_XCatch;   //!
   TBranch        *b_theta_XCatch;   //!
   TBranch        *b_vtx_XCatchCorrDE;   //!
   TBranch        *b_vty_XCatchCorrDE;   //!
   TBranch        *b_vtz_XCatchCorrDE;   //!
   TBranch        *b_closedist_XCatchCorrDE;   //!
   TBranch        *b_theta_XCatchCorrDE;   //!
   TBranch        *b_vtx_WCatch;   //!
   TBranch        *b_vty_WCatch;   //!
   TBranch        *b_vtz_WCatch;   //!
   TBranch        *b_closedist_WCatch;   //!
   TBranch        *b_PiBeamPScatPMom;   //!
   TBranch        *b_PiBeamPScatPiMom;   //!
   TBranch        *b_theta_PiBeamPScatPi;   //!
   TBranch        *b_PiBeamPScatNo;   //!
   TBranch        *b_NBeamMom;   //!
   TBranch        *b_NBeamMomx;   //!
   TBranch        *b_NBeamMomy;   //!
   TBranch        *b_NBeamMomz;   //!
   TBranch        *b_NBeamLength;   //!
   TBranch        *b_vtx_NBeam;   //!
   TBranch        *b_vty_NBeam;   //!
   TBranch        *b_vtz_NBeam;   //!
   TBranch        *b_vdist_NBeam;   //!
   TBranch        *b_cdistNBeam;   //!
   TBranch        *b_pitrno_NBeam;   //!
   TBranch        *b_pikno_NBeam;   //!
   TBranch        *b_DeltaE_NPScat;   //!
   TBranch        *b_ProtonMom_NPScat;   //!
   TBranch        *b_DecayNeutronMom;   //!
   TBranch        *b_DecayPionMom;   //!
   TBranch        *b_MissMassSigmaP_NPScat;   //!
   TBranch        *b_vtx_Decay2np;   //!
   TBranch        *b_vty_Decay2np;   //!
   TBranch        *b_vtz_Decay2np;   //!
   TBranch        *b_vtx_NPScat;   //!
   TBranch        *b_vty_NPScat;   //!
   TBranch        *b_vtz_NPScat;   //!
   TBranch        *b_vdist1_NPScat;   //!
   TBranch        *b_vdist2_NPScat;   //!
   TBranch        *b_cdistDecay2np;   //!
   TBranch        *b_cdistNPScat;   //!
   TBranch        *b_theta_NPScat;   //!
   TBranch        *b_thetaCM_NPScat;   //!
   TBranch        *b_ptrno_NPScat;   //!
   TBranch        *b_pitrno_NPScat;   //!
   TBranch        *b_pikno_NPScat;   //!
   TBranch        *b_priority_NPScat;   //!
   TBranch        *b_DeltaE_SigmaPScat;   //!
   TBranch        *b_MissMassSigmaP_SigmaPScat;   //!
   TBranch        *b_vtx_SigmaPScat;   //!
   TBranch        *b_vty_SigmaPScat;   //!
   TBranch        *b_vtz_SigmaPScat;   //!
   TBranch        *b_vdist_SigmaPScat;   //!
   TBranch        *b_cdistSigmaPScat;   //!
   TBranch        *b_ptrno_SigmaPScat;   //!
   TBranch        *b_pikno_SigmaPScat;   //!
   TBranch        *b_DeltaE_SigmaPScat2npi;   //!
   TBranch        *b_ProtonMom_SigmaPScat2npi;   //!
   TBranch        *b_MissMassSigmaP_SigmaPScat2npi;   //!
   TBranch        *b_vtx_ScatSigmaDecay2npi;   //!
   TBranch        *b_vty_ScatSigmaDecay2npi;   //!
   TBranch        *b_vtz_ScatSigmaDecay2npi;   //!
   TBranch        *b_cdistScatSigmaDecay2npi;   //!
   TBranch        *b_vdist_SigmaPScat2npi;   //!
   TBranch        *b_theta_SigmaPScat2npi;   //!
   TBranch        *b_thetaCM_SigmaPScat2npi;   //!
   TBranch        *b_ptrno_SigmaPScat2npi;   //!
   TBranch        *b_pitrno_SigmaPScat2npi;   //!
   TBranch        *b_pikno_SigmaPScat2npi;   //!
   TBranch        *b_priority_SigmaPScat2npi;   //!
   TBranch        *b_DeltaP_LambdaNConv;   //!
   TBranch        *b_NeutronMom_LambdaNConv;   //!
   TBranch        *b_LambdaMom;   //!
   TBranch        *b_thetaCM_LamdaNConv;   //!
   TBranch        *b_vtx_LambdaNConv;   //!
   TBranch        *b_vty_LambdaNConv;   //!
   TBranch        *b_vtz_LambdaNConv;   //!
   TBranch        *b_vtx_LambdaDecay;   //!
   TBranch        *b_vty_LambdaDecay;   //!
   TBranch        *b_vtz_LambdaDecay;   //!
   TBranch        *b_vdist1_LambdaNConv;   //!
   TBranch        *b_vdist2_LambdaNConv;   //!
   TBranch        *b_cdistLambdaNConv;   //!
   TBranch        *b_cdistLambdaDecay;   //!
   TBranch        *b_ptrno_LambdaNConv;   //!
   TBranch        *b_pitrno_LambdaNConv;   //!
   TBranch        *b_pikno_LambdaNConv;   //!
   TBranch        *b_priority_LambdaNConv;   //!
   TBranch        *b_DeltaP_PiPScat;   //!
   TBranch        *b_DecayNeutronMom2;   //!
   TBranch        *b_DecayPionMom2;   //!
   TBranch        *b_vtx_Decay2pip;   //!
   TBranch        *b_vty_Decay2pip;   //!
   TBranch        *b_vtz_Decay2pip;   //!
   TBranch        *b_vtx_PiPScat;   //!
   TBranch        *b_vty_PiPScat;   //!
   TBranch        *b_vtz_PiPScat;   //!
   TBranch        *b_vdist1_PiPScat;   //!
   TBranch        *b_vdist2_PiPScat;   //!
   TBranch        *b_cdistDecay2pip;   //!
   TBranch        *b_cdistPiPScat;   //!
   TBranch        *b_ptrno_PiPScat;   //!
   TBranch        *b_pitrno_PiPScat;   //!
   TBranch        *b_pikno_PiPScat;   //!
   TBranch        *b_priority_PiPScat;   //!

   DstPiKAna_SM(TTree *tree=0);
   virtual ~DstPiKAna_SM();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DstPiKAna_SM_cxx
DstPiKAna_SM::DstPiKAna_SM(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../dst_rootfile/run06444_DstPiKAna_SM_Sigma.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../dst_rootfile/run06444_DstPiKAna_SM_Sigma.root");
      }
      f->GetObject("pik",tree);

   }
   Init(tree);
}

DstPiKAna_SM::~DstPiKAna_SM()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DstPiKAna_SM::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DstPiKAna_SM::LoadTree(Long64_t entry)
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

void DstPiKAna_SM::Init(TTree *tree)
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

   fChain->SetBranchAddress("runnum", &runnum, &b_runnum);
   fChain->SetBranchAddress("evnum", &evnum, &b_evnum);
   fChain->SetBranchAddress("spill", &spill, &b_spill);
   fChain->SetBranchAddress("trignhits", &trignhits, &b_trignhits);
   fChain->SetBranchAddress("trigpat", trigpat, &b_trigpat);
   fChain->SetBranchAddress("trigflag", trigflag, &b_trigflag);
   fChain->SetBranchAddress("nhBh1", &nhBh1, &b_nhBh1);
   fChain->SetBranchAddress("csBh1", csBh1, &b_csBh1);
   fChain->SetBranchAddress("Bh1Seg", Bh1Seg, &b_Bh1Seg);
   fChain->SetBranchAddress("tBh1", tBh1, &b_tBh1);
   fChain->SetBranchAddress("dtBh1", dtBh1, &b_dtBh1);
   fChain->SetBranchAddress("deBh1", deBh1, &b_deBh1);
   fChain->SetBranchAddress("btof", btof, &b_btof);
   fChain->SetBranchAddress("nhBh2", &nhBh2, &b_nhBh2);
   fChain->SetBranchAddress("csBh2", csBh2, &b_csBh2);
   fChain->SetBranchAddress("Bh2Seg", Bh2Seg, &b_Bh2Seg);
   fChain->SetBranchAddress("tBh2", tBh2, &b_tBh2);
   fChain->SetBranchAddress("t0Bh2", t0Bh2, &b_t0Bh2);
   fChain->SetBranchAddress("dtBh2", dtBh2, &b_dtBh2);
   fChain->SetBranchAddress("deBh2", deBh2, &b_deBh2);
   fChain->SetBranchAddress("nhTof", &nhTof, &b_nhTof);
   fChain->SetBranchAddress("csTof", csTof, &b_csTof);
   fChain->SetBranchAddress("TofSeg", TofSeg, &b_TofSeg);
   fChain->SetBranchAddress("tTof", tTof, &b_tTof);
   fChain->SetBranchAddress("dtTof", dtTof, &b_dtTof);
   fChain->SetBranchAddress("deTof", deTof, &b_deTof);
   fChain->SetBranchAddress("bft_ncl", &bft_ncl, &b_bft_ncl);
   fChain->SetBranchAddress("nlBcOut", &nlBcOut, &b_nlBcOut);
   fChain->SetBranchAddress("ntBcOut", &ntBcOut, &b_ntBcOut);
   fChain->SetBranchAddress("nhBcOut", nhBcOut, &b_nhBcOut);
   fChain->SetBranchAddress("chisqrBcOut", chisqrBcOut, &b_chisqrBcOut);
   fChain->SetBranchAddress("x0BcOut", x0BcOut, &b_x0BcOut);
   fChain->SetBranchAddress("y0BcOut", y0BcOut, &b_y0BcOut);
   fChain->SetBranchAddress("u0BcOut", u0BcOut, &b_u0BcOut);
   fChain->SetBranchAddress("v0BcOut", v0BcOut, &b_v0BcOut);
   fChain->SetBranchAddress("xtgtBcOut", xtgtBcOut, &b_xtgtBcOut);
   fChain->SetBranchAddress("ytgtBcOut", ytgtBcOut, &b_ytgtBcOut);
   fChain->SetBranchAddress("xbh2BcOut", xbh2BcOut, &b_xbh2BcOut);
   fChain->SetBranchAddress("ybh2BcOut", ybh2BcOut, &b_ybh2BcOut);
   fChain->SetBranchAddress("ntK18", &ntK18, &b_ntK18);
   fChain->SetBranchAddress("nhK18", nhK18, &b_nhK18);
   fChain->SetBranchAddress("chisqrK18", chisqrK18, &b_chisqrK18);
   fChain->SetBranchAddress("pK18", pK18, &b_pK18);
   fChain->SetBranchAddress("xtgtK18", xtgtK18, &b_xtgtK18);
   fChain->SetBranchAddress("ytgtK18", ytgtK18, &b_ytgtK18);
   fChain->SetBranchAddress("utgtK18", utgtK18, &b_utgtK18);
   fChain->SetBranchAddress("vtgtK18", vtgtK18, &b_vtgtK18);
   fChain->SetBranchAddress("thetaK18", thetaK18, &b_thetaK18);
   fChain->SetBranchAddress("nlSdcIn", &nlSdcIn, &b_nlSdcIn);
   fChain->SetBranchAddress("ntSdcIn", &ntSdcIn, &b_ntSdcIn);
   fChain->SetBranchAddress("nhSdcIn", nhSdcIn, &b_nhSdcIn);
   fChain->SetBranchAddress("chisqrSdcIn", chisqrSdcIn, &b_chisqrSdcIn);
   fChain->SetBranchAddress("x0SdcIn", x0SdcIn, &b_x0SdcIn);
   fChain->SetBranchAddress("y0SdcIn", y0SdcIn, &b_y0SdcIn);
   fChain->SetBranchAddress("u0SdcIn", u0SdcIn, &b_u0SdcIn);
   fChain->SetBranchAddress("v0SdcIn", v0SdcIn, &b_v0SdcIn);
   fChain->SetBranchAddress("nlSdcOut", &nlSdcOut, &b_nlSdcOut);
   fChain->SetBranchAddress("ntSdcOut", &ntSdcOut, &b_ntSdcOut);
   fChain->SetBranchAddress("nhSdcOut", nhSdcOut, &b_nhSdcOut);
   fChain->SetBranchAddress("chisqrSdcOut", chisqrSdcOut, &b_chisqrSdcOut);
   fChain->SetBranchAddress("x0SdcOut", x0SdcOut, &b_x0SdcOut);
   fChain->SetBranchAddress("y0SdcOut", y0SdcOut, &b_y0SdcOut);
   fChain->SetBranchAddress("u0SdcOut", u0SdcOut, &b_u0SdcOut);
   fChain->SetBranchAddress("v0SdcOut", v0SdcOut, &b_v0SdcOut);
   fChain->SetBranchAddress("ntKurama", &ntKurama, &b_ntKurama);
   fChain->SetBranchAddress("nhKurama", nhKurama, &b_nhKurama);
   fChain->SetBranchAddress("stof", stof, &b_stof);
   fChain->SetBranchAddress("cstof", cstof, &b_cstof);
   fChain->SetBranchAddress("path", path, &b_path);
   fChain->SetBranchAddress("chisqrKurama", chisqrKurama, &b_chisqrKurama);
   fChain->SetBranchAddress("pKurama", pKurama, &b_pKurama);
   fChain->SetBranchAddress("qKurama", qKurama, &b_qKurama);
   fChain->SetBranchAddress("m2", m2, &b_m2);
   fChain->SetBranchAddress("cm2", cm2, &b_cm2);
   fChain->SetBranchAddress("xtgtKurama", xtgtKurama, &b_xtgtKurama);
   fChain->SetBranchAddress("ytgtKurama", ytgtKurama, &b_ytgtKurama);
   fChain->SetBranchAddress("utgtKurama", utgtKurama, &b_utgtKurama);
   fChain->SetBranchAddress("vtgtKurama", vtgtKurama, &b_vtgtKurama);
   fChain->SetBranchAddress("thetaKurama", thetaKurama, &b_thetaKurama);
   fChain->SetBranchAddress("xtofKurama", xtofKurama, &b_xtofKurama);
   fChain->SetBranchAddress("ytofKurama", ytofKurama, &b_ytofKurama);
   fChain->SetBranchAddress("utofKurama", utofKurama, &b_utofKurama);
   fChain->SetBranchAddress("vtofKurama", vtofKurama, &b_vtofKurama);
   fChain->SetBranchAddress("tofsegKurama", tofsegKurama, &b_tofsegKurama);
   fChain->SetBranchAddress("best_deTof", best_deTof, &b_best_deTof);
   fChain->SetBranchAddress("best_TofSeg", best_TofSeg, &b_best_TofSeg);
   fChain->SetBranchAddress("ntCFT", &ntCFT, &b_ntCFT);
   fChain->SetBranchAddress("ntProton", &ntProton, &b_ntProton);
   fChain->SetBranchAddress("ntOther", &ntOther, &b_ntOther);
   fChain->SetBranchAddress("theta_cft", theta_cft, &b_theta_cft);
   fChain->SetBranchAddress("phi_cft", phi_cft, &b_phi_cft);
   fChain->SetBranchAddress("vtx_cft", vtx_cft, &b_vtx_cft);
   fChain->SetBranchAddress("vty_cft", vty_cft, &b_vty_cft);
   fChain->SetBranchAddress("vtz_cft", vtz_cft, &b_vtz_cft);
   fChain->SetBranchAddress("Total_dE", Total_dE, &b_Total_dE);
   fChain->SetBranchAddress("Total_dEphi", Total_dEphi, &b_Total_dEphi);
   fChain->SetBranchAddress("Total_dEuv", Total_dEuv, &b_Total_dEuv);
   fChain->SetBranchAddress("Total_dE_max", Total_dE_max, &b_Total_dE_max);
   fChain->SetBranchAddress("Total_dEphi_max", Total_dEphi_max, &b_Total_dEphi_max);
   fChain->SetBranchAddress("Total_dEuv_max", Total_dEuv_max, &b_Total_dEuv_max);
   fChain->SetBranchAddress("nhit_phi", nhit_phi, &b_nhit_phi);
   fChain->SetBranchAddress("nhit_uv", nhit_uv, &b_nhit_uv);
   fChain->SetBranchAddress("segBGOt", segBGOt, &b_segBGOt);
   fChain->SetBranchAddress("energyBGOt", energyBGOt, &b_energyBGOt);
   fChain->SetBranchAddress("segPiIDt", segPiIDt, &b_segPiIDt);
   fChain->SetBranchAddress("protonflagt", protonflagt, &b_protonflagt);
   fChain->SetBranchAddress("BGOnohitt", BGOnohitt, &b_BGOnohitt);
   fChain->SetBranchAddress("simKuramat", simKuramat, &b_simKuramat);
   fChain->SetBranchAddress("tracknoproton", tracknoproton, &b_tracknoproton);
   fChain->SetBranchAddress("tracknoother", tracknoother, &b_tracknoother);
   fChain->SetBranchAddress("nPi", &nPi, &b_nPi);
   fChain->SetBranchAddress("nK", &nK, &b_nK);
   fChain->SetBranchAddress("nPiK", &nPiK, &b_nPiK);
   fChain->SetBranchAddress("nCatch", &nCatch, &b_nCatch);
   fChain->SetBranchAddress("nPiKCatch", &nPiKCatch, &b_nPiKCatch);
   fChain->SetBranchAddress("vtx_KURAMA", vtx_KURAMA, &b_vtx_KURAMA);
   fChain->SetBranchAddress("vty_KURAMA", vty_KURAMA, &b_vty_KURAMA);
   fChain->SetBranchAddress("vtz_KURAMA", vtz_KURAMA, &b_vtz_KURAMA);
   fChain->SetBranchAddress("ccm2", ccm2, &b_ccm2);
   fChain->SetBranchAddress("cpath", cpath, &b_cpath);
   fChain->SetBranchAddress("chisqrKuramapik", chisqrKuramapik, &b_chisqrKuramapik);
   fChain->SetBranchAddress("pKuramapik", pKuramapik, &b_pKuramapik);
   fChain->SetBranchAddress("qKuramapik", qKuramapik, &b_qKuramapik);
   fChain->SetBranchAddress("simKuramatpik", simKuramatpik, &b_simKuramatpik);
   fChain->SetBranchAddress("Kflag", Kflag, &b_Kflag);
   fChain->SetBranchAddress("closedist_KURAMA", closedist_KURAMA, &b_closedist_KURAMA);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("MissMass", MissMass, &b_MissMass);
   fChain->SetBranchAddress("MissMassCorr", MissMassCorr, &b_MissMassCorr);
   fChain->SetBranchAddress("MissMassCorrDE", MissMassCorrDE, &b_MissMassCorrDE);
   fChain->SetBranchAddress("MissMass2CorrDE", MissMass2CorrDE, &b_MissMass2CorrDE);
   fChain->SetBranchAddress("MissMassPiPi", MissMassPiPi, &b_MissMassPiPi);
   fChain->SetBranchAddress("MissMassPiPiCorrDE", MissMassPiPiCorrDE, &b_MissMassPiPiCorrDE);
   fChain->SetBranchAddress("MissMass2PiPiCorrDE", MissMass2PiPiCorrDE, &b_MissMass2PiPiCorrDE);
   fChain->SetBranchAddress("pCorrDE_PiPi", pCorrDE_PiPi, &b_pCorrDE_PiPi);
   fChain->SetBranchAddress("pCalcCorrDE_PiPi", pCalcCorrDE_PiPi, &b_pCalcCorrDE_PiPi);
   fChain->SetBranchAddress("MissMassPiP", MissMassPiP, &b_MissMassPiP);
   fChain->SetBranchAddress("MissMassPiPCorrDE", MissMassPiPCorrDE, &b_MissMassPiPCorrDE);
   fChain->SetBranchAddress("MissMass2PiPCorrDE", MissMass2PiPCorrDE, &b_MissMass2PiPCorrDE);
   fChain->SetBranchAddress("pCorrDE_PiP", pCorrDE_PiP, &b_pCorrDE_PiP);
   fChain->SetBranchAddress("pCalcCorrDE_PiP", pCalcCorrDE_PiP, &b_pCalcCorrDE_PiP);
   fChain->SetBranchAddress("thetaCM", thetaCM, &b_thetaCM);
   fChain->SetBranchAddress("thetaCMCorrDE", thetaCMCorrDE, &b_thetaCMCorrDE);
   fChain->SetBranchAddress("thetaCMCorrDE_PiPi", thetaCMCorrDE_PiPi, &b_thetaCMCorrDE_PiPi);
   fChain->SetBranchAddress("thetaCMCorrDE_PiP", thetaCMCorrDE_PiP, &b_thetaCMCorrDE_PiP);
   fChain->SetBranchAddress("costCM", costCM, &b_costCM);
   fChain->SetBranchAddress("KURAMAPID", KURAMAPID, &b_KURAMAPID);
   fChain->SetBranchAddress("MissMom", MissMom, &b_MissMom);
   fChain->SetBranchAddress("MissMomx", MissMomx, &b_MissMomx);
   fChain->SetBranchAddress("MissMomy", MissMomy, &b_MissMomy);
   fChain->SetBranchAddress("MissMomz", MissMomz, &b_MissMomz);
   fChain->SetBranchAddress("MissMomcal", MissMomcal, &b_MissMomcal);
   fChain->SetBranchAddress("MissMomxcal", MissMomxcal, &b_MissMomxcal);
   fChain->SetBranchAddress("MissMomycal", MissMomycal, &b_MissMomycal);
   fChain->SetBranchAddress("MissMomzcal", MissMomzcal, &b_MissMomzcal);
   fChain->SetBranchAddress("SigmaBeta", SigmaBeta, &b_SigmaBeta);
   fChain->SetBranchAddress("xpi", xpi, &b_xpi);
   fChain->SetBranchAddress("ypi", ypi, &b_ypi);
   fChain->SetBranchAddress("upi", upi, &b_upi);
   fChain->SetBranchAddress("vpi", vpi, &b_vpi);
   fChain->SetBranchAddress("xk", xk, &b_xk);
   fChain->SetBranchAddress("yk", yk, &b_yk);
   fChain->SetBranchAddress("uk", uk, &b_uk);
   fChain->SetBranchAddress("vk", vk, &b_vk);
   fChain->SetBranchAddress("uc", uc, &b_uc);
   fChain->SetBranchAddress("vc", vc, &b_vc);
   fChain->SetBranchAddress("pOrg", pOrg, &b_pOrg);
   fChain->SetBranchAddress("pCalc", pCalc, &b_pCalc);
   fChain->SetBranchAddress("pCalcCorrDE", pCalcCorrDE, &b_pCalcCorrDE);
   fChain->SetBranchAddress("pCorr", pCorr, &b_pCorr);
   fChain->SetBranchAddress("pCorrDE", pCorrDE, &b_pCorrDE);
   fChain->SetBranchAddress("pCorrDE_PiBeam", pCorrDE_PiBeam, &b_pCorrDE_PiBeam);
   fChain->SetBranchAddress("vtx_K18Catch", vtx_K18Catch, &b_vtx_K18Catch);
   fChain->SetBranchAddress("vty_K18Catch", vty_K18Catch, &b_vty_K18Catch);
   fChain->SetBranchAddress("vtz_K18Catch", vtz_K18Catch, &b_vtz_K18Catch);
   fChain->SetBranchAddress("closedist_K18Catch", closedist_K18Catch, &b_closedist_K18Catch);
   fChain->SetBranchAddress("theta_K18Catch", theta_K18Catch, &b_theta_K18Catch);
   fChain->SetBranchAddress("vertex_distance", vertex_distance, &b_vertex_distance);
   fChain->SetBranchAddress("vtx_XCatch", vtx_XCatch, &b_vtx_XCatch);
   fChain->SetBranchAddress("vty_XCatch", vty_XCatch, &b_vty_XCatch);
   fChain->SetBranchAddress("vtz_XCatch", vtz_XCatch, &b_vtz_XCatch);
   fChain->SetBranchAddress("closedist_XCatch", closedist_XCatch, &b_closedist_XCatch);
   fChain->SetBranchAddress("theta_XCatch", theta_XCatch, &b_theta_XCatch);
   fChain->SetBranchAddress("vtx_XCatchCorrDE", vtx_XCatchCorrDE, &b_vtx_XCatchCorrDE);
   fChain->SetBranchAddress("vty_XCatchCorrDE", vty_XCatchCorrDE, &b_vty_XCatchCorrDE);
   fChain->SetBranchAddress("vtz_XCatchCorrDE", vtz_XCatchCorrDE, &b_vtz_XCatchCorrDE);
   fChain->SetBranchAddress("closedist_XCatchCorrDE", closedist_XCatchCorrDE, &b_closedist_XCatchCorrDE);
   fChain->SetBranchAddress("theta_XCatchCorrDE", theta_XCatchCorrDE, &b_theta_XCatchCorrDE);
   fChain->SetBranchAddress("vtx_WCatch", vtx_WCatch, &b_vtx_WCatch);
   fChain->SetBranchAddress("vty_WCatch", vty_WCatch, &b_vty_WCatch);
   fChain->SetBranchAddress("vtz_WCatch", vtz_WCatch, &b_vtz_WCatch);
   fChain->SetBranchAddress("closedist_WCatch", closedist_WCatch, &b_closedist_WCatch);
   fChain->SetBranchAddress("PiBeamPScatPMom", &PiBeamPScatPMom, &b_PiBeamPScatPMom);
   fChain->SetBranchAddress("PiBeamPScatPiMom", &PiBeamPScatPiMom, &b_PiBeamPScatPiMom);
   fChain->SetBranchAddress("theta_PiBeamPScatPi", &theta_PiBeamPScatPi, &b_theta_PiBeamPScatPi);
   fChain->SetBranchAddress("PiBeamPScatNo", &PiBeamPScatNo, &b_PiBeamPScatNo);
   fChain->SetBranchAddress("NBeamMom", &NBeamMom, &b_NBeamMom);
   fChain->SetBranchAddress("NBeamMomx", &NBeamMomx, &b_NBeamMomx);
   fChain->SetBranchAddress("NBeamMomy", &NBeamMomy, &b_NBeamMomy);
   fChain->SetBranchAddress("NBeamMomz", &NBeamMomz, &b_NBeamMomz);
   fChain->SetBranchAddress("NBeamLength", &NBeamLength, &b_NBeamLength);
   fChain->SetBranchAddress("vtx_NBeam", &vtx_NBeam, &b_vtx_NBeam);
   fChain->SetBranchAddress("vty_NBeam", &vty_NBeam, &b_vty_NBeam);
   fChain->SetBranchAddress("vtz_NBeam", &vtz_NBeam, &b_vtz_NBeam);
   fChain->SetBranchAddress("vdist_NBeam", &vdist_NBeam, &b_vdist_NBeam);
   fChain->SetBranchAddress("cdistNBeam", &cdistNBeam, &b_cdistNBeam);
   fChain->SetBranchAddress("pitrno_NBeam", &pitrno_NBeam, &b_pitrno_NBeam);
   fChain->SetBranchAddress("pikno_NBeam", &pikno_NBeam, &b_pikno_NBeam);
   fChain->SetBranchAddress("DeltaE_NPScat", &DeltaE_NPScat, &b_DeltaE_NPScat);
   fChain->SetBranchAddress("ProtonMom_NPScat", &ProtonMom_NPScat, &b_ProtonMom_NPScat);
   fChain->SetBranchAddress("DecayNeutronMom", &DecayNeutronMom, &b_DecayNeutronMom);
   fChain->SetBranchAddress("DecayPionMom", &DecayPionMom, &b_DecayPionMom);
   fChain->SetBranchAddress("MissMassSigmaP_NPScat", &MissMassSigmaP_NPScat, &b_MissMassSigmaP_NPScat);
   fChain->SetBranchAddress("vtx_Decay2np", &vtx_Decay2np, &b_vtx_Decay2np);
   fChain->SetBranchAddress("vty_Decay2np", &vty_Decay2np, &b_vty_Decay2np);
   fChain->SetBranchAddress("vtz_Decay2np", &vtz_Decay2np, &b_vtz_Decay2np);
   fChain->SetBranchAddress("vtx_NPScat", &vtx_NPScat, &b_vtx_NPScat);
   fChain->SetBranchAddress("vty_NPScat", &vty_NPScat, &b_vty_NPScat);
   fChain->SetBranchAddress("vtz_NPScat", &vtz_NPScat, &b_vtz_NPScat);
   fChain->SetBranchAddress("vdist1_NPScat", &vdist1_NPScat, &b_vdist1_NPScat);
   fChain->SetBranchAddress("vdist2_NPScat", &vdist2_NPScat, &b_vdist2_NPScat);
   fChain->SetBranchAddress("cdistDecay2np", &cdistDecay2np, &b_cdistDecay2np);
   fChain->SetBranchAddress("cdistNPScat", &cdistNPScat, &b_cdistNPScat);
   fChain->SetBranchAddress("theta_NPScat", &theta_NPScat, &b_theta_NPScat);
   fChain->SetBranchAddress("thetaCM_NPScat", &thetaCM_NPScat, &b_thetaCM_NPScat);
   fChain->SetBranchAddress("ptrno_NPScat", &ptrno_NPScat, &b_ptrno_NPScat);
   fChain->SetBranchAddress("pitrno_NPScat", &pitrno_NPScat, &b_pitrno_NPScat);
   fChain->SetBranchAddress("pikno_NPScat", &pikno_NPScat, &b_pikno_NPScat);
   fChain->SetBranchAddress("priority_NPScat", &priority_NPScat, &b_priority_NPScat);
   fChain->SetBranchAddress("DeltaE_SigmaPScat", &DeltaE_SigmaPScat, &b_DeltaE_SigmaPScat);
   fChain->SetBranchAddress("MissMassSigmaP_SigmaPScat", &MissMassSigmaP_SigmaPScat, &b_MissMassSigmaP_SigmaPScat);
   fChain->SetBranchAddress("vtx_SigmaPScat", &vtx_SigmaPScat, &b_vtx_SigmaPScat);
   fChain->SetBranchAddress("vty_SigmaPScat", &vty_SigmaPScat, &b_vty_SigmaPScat);
   fChain->SetBranchAddress("vtz_SigmaPScat", &vtz_SigmaPScat, &b_vtz_SigmaPScat);
   fChain->SetBranchAddress("vdist_SigmaPScat", &vdist_SigmaPScat, &b_vdist_SigmaPScat);
   fChain->SetBranchAddress("cdistSigmaPScat", &cdistSigmaPScat, &b_cdistSigmaPScat);
   fChain->SetBranchAddress("ptrno_SigmaPScat", &ptrno_SigmaPScat, &b_ptrno_SigmaPScat);
   fChain->SetBranchAddress("pikno_SigmaPScat", &pikno_SigmaPScat, &b_pikno_SigmaPScat);
   fChain->SetBranchAddress("DeltaE_SigmaPScat2npi", &DeltaE_SigmaPScat2npi, &b_DeltaE_SigmaPScat2npi);
   fChain->SetBranchAddress("ProtonMom_SigmaPScat2npi", &ProtonMom_SigmaPScat2npi, &b_ProtonMom_SigmaPScat2npi);
   fChain->SetBranchAddress("MissMassSigmaP_SigmaPScat2npi", &MissMassSigmaP_SigmaPScat2npi, &b_MissMassSigmaP_SigmaPScat2npi);
   fChain->SetBranchAddress("vtx_ScatSigmaDecay2npi", &vtx_ScatSigmaDecay2npi, &b_vtx_ScatSigmaDecay2npi);
   fChain->SetBranchAddress("vty_ScatSigmaDecay2npi", &vty_ScatSigmaDecay2npi, &b_vty_ScatSigmaDecay2npi);
   fChain->SetBranchAddress("vtz_ScatSigmaDecay2npi", &vtz_ScatSigmaDecay2npi, &b_vtz_ScatSigmaDecay2npi);
   fChain->SetBranchAddress("cdistScatSigmaDecay2npi", &cdistScatSigmaDecay2npi, &b_cdistScatSigmaDecay2npi);
   fChain->SetBranchAddress("vdist_SigmaPScat2npi", &vdist_SigmaPScat2npi, &b_vdist_SigmaPScat2npi);
   fChain->SetBranchAddress("theta_SigmaPScat2npi", &theta_SigmaPScat2npi, &b_theta_SigmaPScat2npi);
   fChain->SetBranchAddress("thetaCM_SigmaPScat2npi", &thetaCM_SigmaPScat2npi, &b_thetaCM_SigmaPScat2npi);
   fChain->SetBranchAddress("ptrno_SigmaPScat2npi", &ptrno_SigmaPScat2npi, &b_ptrno_SigmaPScat2npi);
   fChain->SetBranchAddress("pitrno_SigmaPScat2npi", &pitrno_SigmaPScat2npi, &b_pitrno_SigmaPScat2npi);
   fChain->SetBranchAddress("pikno_SigmaPScat2npi", &pikno_SigmaPScat2npi, &b_pikno_SigmaPScat2npi);
   fChain->SetBranchAddress("priority_SigmaPScat2npi", &priority_SigmaPScat2npi, &b_priority_SigmaPScat2npi);
   fChain->SetBranchAddress("DeltaP_LambdaNConv", &DeltaP_LambdaNConv, &b_DeltaP_LambdaNConv);
   fChain->SetBranchAddress("NeutronMom_LambdaNConv", &NeutronMom_LambdaNConv, &b_NeutronMom_LambdaNConv);
   fChain->SetBranchAddress("LambdaMom", &LambdaMom, &b_LambdaMom);
   fChain->SetBranchAddress("thetaCM_LamdaNConv", &thetaCM_LamdaNConv, &b_thetaCM_LamdaNConv);
   fChain->SetBranchAddress("vtx_LambdaNConv", &vtx_LambdaNConv, &b_vtx_LambdaNConv);
   fChain->SetBranchAddress("vty_LambdaNConv", &vty_LambdaNConv, &b_vty_LambdaNConv);
   fChain->SetBranchAddress("vtz_LambdaNConv", &vtz_LambdaNConv, &b_vtz_LambdaNConv);
   fChain->SetBranchAddress("vtx_LambdaDecay", &vtx_LambdaDecay, &b_vtx_LambdaDecay);
   fChain->SetBranchAddress("vty_LambdaDecay", &vty_LambdaDecay, &b_vty_LambdaDecay);
   fChain->SetBranchAddress("vtz_LambdaDecay", &vtz_LambdaDecay, &b_vtz_LambdaDecay);
   fChain->SetBranchAddress("vdist1_LambdaNConv", &vdist1_LambdaNConv, &b_vdist1_LambdaNConv);
   fChain->SetBranchAddress("vdist2_LambdaNConv", &vdist2_LambdaNConv, &b_vdist2_LambdaNConv);
   fChain->SetBranchAddress("cdistLambdaNConv", &cdistLambdaNConv, &b_cdistLambdaNConv);
   fChain->SetBranchAddress("cdistLambdaDecay", &cdistLambdaDecay, &b_cdistLambdaDecay);
   fChain->SetBranchAddress("ptrno_LambdaNConv", &ptrno_LambdaNConv, &b_ptrno_LambdaNConv);
   fChain->SetBranchAddress("pitrno_LambdaNConv", &pitrno_LambdaNConv, &b_pitrno_LambdaNConv);
   fChain->SetBranchAddress("pikno_LambdaNConv", &pikno_LambdaNConv, &b_pikno_LambdaNConv);
   fChain->SetBranchAddress("priority_LambdaNConv", &priority_LambdaNConv, &b_priority_LambdaNConv);
   fChain->SetBranchAddress("DeltaP_PiPScat", &DeltaP_PiPScat, &b_DeltaP_PiPScat);
   fChain->SetBranchAddress("DecayNeutronMom2", &DecayNeutronMom2, &b_DecayNeutronMom2);
   fChain->SetBranchAddress("DecayPionMom2", &DecayPionMom2, &b_DecayPionMom2);
   fChain->SetBranchAddress("vtx_Decay2pip", &vtx_Decay2pip, &b_vtx_Decay2pip);
   fChain->SetBranchAddress("vty_Decay2pip", &vty_Decay2pip, &b_vty_Decay2pip);
   fChain->SetBranchAddress("vtz_Decay2pip", &vtz_Decay2pip, &b_vtz_Decay2pip);
   fChain->SetBranchAddress("vtx_PiPScat", &vtx_PiPScat, &b_vtx_PiPScat);
   fChain->SetBranchAddress("vty_PiPScat", &vty_PiPScat, &b_vty_PiPScat);
   fChain->SetBranchAddress("vtz_PiPScat", &vtz_PiPScat, &b_vtz_PiPScat);
   fChain->SetBranchAddress("vdist1_PiPScat", &vdist1_PiPScat, &b_vdist1_PiPScat);
   fChain->SetBranchAddress("vdist2_PiPScat", &vdist2_PiPScat, &b_vdist2_PiPScat);
   fChain->SetBranchAddress("cdistDecay2pip", &cdistDecay2pip, &b_cdistDecay2pip);
   fChain->SetBranchAddress("cdistPiPScat", &cdistPiPScat, &b_cdistPiPScat);
   fChain->SetBranchAddress("ptrno_PiPScat", &ptrno_PiPScat, &b_ptrno_PiPScat);
   fChain->SetBranchAddress("pitrno_PiPScat", &pitrno_PiPScat, &b_pitrno_PiPScat);
   fChain->SetBranchAddress("pikno_PiPScat", &pikno_PiPScat, &b_pikno_PiPScat);
   fChain->SetBranchAddress("priority_PiPScat", &priority_PiPScat, &b_priority_PiPScat);
   Notify();
}

Bool_t DstPiKAna_SM::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DstPiKAna_SM::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DstPiKAna_SM::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DstPiKAna_SM_cxx
