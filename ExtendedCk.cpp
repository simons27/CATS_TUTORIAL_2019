
#include "Basics.h"
#include "ExtendedCk.h"
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"

#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"

#include "TGraph.h"
#include "TString.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"

void Ck_pL_Ledni_Usmani(){
    const unsigned NumMomBins = 60;
    const double kMin = 0;
    const double kMax = 240;

    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource,2,true);
    SOURCE_PARS.SetParameter(0,1.2);

    CATS Kitty_pL;
    Kitty_pL.SetMomBins(NumMomBins,kMin,kMax);
    Kitty_pL.SetAnaSource(GaussSource,SOURCE_PARS);
    Kitty_pL.SetAnaSource(0,SOURCE_PARS.GetParameter(0));
    Kitty_pL.SetUseAnalyticSource(true);
    Kitty_pL.SetMomentumDependentSource(false);
    Kitty_pL.SetThetaDependentSource(false);
    Kitty_pL.SetExcludeFailedBins(false);
    Kitty_pL.SetQ1Q2(0);
    Kitty_pL.SetPdgId(2212, 3122);
    Kitty_pL.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );
    Kitty_pL.SetNumChannels(2);
    Kitty_pL.SetNumPW(0,1);
    Kitty_pL.SetNumPW(1,1);
    Kitty_pL.SetSpin(0,0);
    Kitty_pL.SetSpin(1,1);
    Kitty_pL.SetChannelWeight(0, 1./4.);
    Kitty_pL.SetChannelWeight(1, 3./4.);

    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double PotPars1S0[8]={pL_UsmaniOli,0,0,0,0,0,0,0};
    double PotPars3S1[8]={pL_UsmaniOli,0,0,0,0,1,0,1};
    CATSparameters cPotPars1S0(CATSparameters::tPotential,8,true);
    cPotPars1S0.SetParameters(PotPars1S0);
    CATSparameters cPotPars3S1(CATSparameters::tPotential,8,true);
    cPotPars3S1.SetParameters(PotPars3S1);

    Kitty_pL.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0);
    Kitty_pL.SetShortRangePotential(1,0,fDlmPot,cPotPars3S1);
    //Kitty_pL.SetMaxNumThreads(4);
    Kitty_pL.KillTheCat();

    DLM_Ck Ck_Usmani(1,0,Kitty_pL);
    Ck_Usmani.Update();

    DLM_Ck Ck_Lednicky(1,4,NumMomBins,kMin,kMax,Lednicky_SingletTriplet);
    Ck_Lednicky.SetSourcePar(0,SOURCE_PARS.GetParameter(0));
    Ck_Lednicky.SetPotPar(0,2.88);
    Ck_Lednicky.SetPotPar(1,2.92);
    Ck_Lednicky.SetPotPar(2,1.66);
    Ck_Lednicky.SetPotPar(3,3.78);
    Ck_Lednicky.Update();

    RootFile_DlmCk("Ck_pL_Ledni_Usmani.root","gCk_pL_Usmani",&Ck_Usmani);
    RootFile_DlmCk("Ck_pL_Ledni_Usmani.root","gCk_pL_Lednicky",&Ck_Lednicky);
}

void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, DLM_Ck* CkToPlot){
    TFile* RootFile = new TFile(RootFileName,"update");
    if(!RootFile) RootFile = new TFile(RootFileName,"recreate");
    const unsigned NumBins = CkToPlot->GetNbins();
    TGraph graph;
    graph.SetName(GraphName);
    graph.Set(NumBins);
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        graph.SetPoint(uBin,CkToPlot->GetBinCenter(0,uBin),CkToPlot->GetBinContent(uBin));
    }
    graph.Write("",TObject::kOverwrite);
    delete RootFile;
}

void Ck_pp_Decomposition(){

}
