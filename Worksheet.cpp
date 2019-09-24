
#include "Worksheet.h"
#include "Basics.h"
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"
#include "TGraph.h"
#include "TString.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"

CATS* KITTY_CATS_WORKFIT_PL;

double CATS_WORKFIT_PL(double* x, double* pars){
    double& Norm = pars[0];
    double& LambdaPar = pars[1];
    double& SourceSize = pars[2];
    //set the radius to the fit value
    //the last parameter says its a small change of the radius, which does not require a new computing grid (saves time)
    //however, the last step requires making sure a good initial value of the radius
    KITTY_CATS_WORKFIT_PL->SetAnaSource(0,SourceSize,true);
    //useful tip: this makes CATS to shut up and not flood your screen. Only errors will be displayed. Use nSilent to suppress even those
    KITTY_CATS_WORKFIT_PL->SetNotifications(CATS::nError);
    KITTY_CATS_WORKFIT_PL->KillTheCat();
    return Norm*(LambdaPar*KITTY_CATS_WORKFIT_PL->EvalCorrFun(*x)+1.-LambdaPar);
}


void Worksheet_ProtonLambda(){

    const unsigned NumMomBins = 60;
    const double kMin = 0;
    const double kMax = 240;

    CATSparameters POT_PARS_1S0(CATSparameters::tPotential,4,true);
    POT_PARS_1S0.SetParameter(0,0);
    POT_PARS_1S0.SetParameter(1,2137);
    POT_PARS_1S0.SetParameter(2,0.5);
    POT_PARS_1S0.SetParameter(3,0.2);
    CATSparameters POT_PARS_3S1(CATSparameters::tPotential,4,true);
    POT_PARS_3S1.SetParameter(0,1);
    POT_PARS_3S1.SetParameter(1,2137);
    POT_PARS_3S1.SetParameter(2,0.5);
    POT_PARS_3S1.SetParameter(3,0.2);

    CATS Kitty_pL;

    //DEFINE YOUR CATS OBJECT HERE
    //USE the potential function:
    Kitty_pL.SetShortRangePotential(0,0,Basics_Potential_Usmani,POT_PARS_1S0);
    Kitty_pL.SetShortRangePotential(1,0,Basics_Potential_Usmani,POT_PARS_3S1);

    Kitty_pL.KillTheCat();

    KITTY_CATS_WORKFIT_PL = &Kitty_pL;

    TFile* fInput = new TFile("../Files/DummyProtonLambda.root","read");
    TH1F* hInput = (TH1F*)fInput->Get("hDummyProtonLambda");
    TF1* FITTER = new TF1("FITTER",CATS_WORKFIT_PL,kMin,kMax,3);
    FITTER->SetParameter(0,1);FITTER->SetParLimits(0,0.5,2);
    FITTER->SetParameter(1,0.5);FITTER->SetParLimits(1,0,1);
    //! SET THE SOURCE PAR
    //FITTER->SetParameter(2,SOURCE_PARS.GetParameter(0));FITTER->SetParLimits(2,0.5,3);


    //FITTER->SetParameter(3,POT_PARS_1S0.GetParameter(1));FITTER->SetParLimits(3,0.99*POT_PARS_1S0.GetParameter(1),1.01*POT_PARS_1S0.GetParameter(1));
    //FITTER->SetParameter(4,POT_PARS_1S0.GetParameter(2));FITTER->SetParLimits(4,0.99*POT_PARS_1S0.GetParameter(2),1.01*POT_PARS_1S0.GetParameter(2));
    //FITTER->SetParameter(5,POT_PARS_1S0.GetParameter(3));FITTER->SetParLimits(5,0.99*POT_PARS_1S0.GetParameter(3),1.01*POT_PARS_1S0.GetParameter(3));

    //FITTER->FixParameter(1,0.6);
    //FITTER->FixParameter(3,POT_PARS_1S0.GetParameter(1));
    //FITTER->FixParameter(4,POT_PARS_1S0.GetParameter(2));
    //FITTER->FixParameter(5,POT_PARS_1S0.GetParameter(3));


    hInput->Fit(FITTER,"S, N, R, M");

    delete FITTER;
    delete fInput;

}
