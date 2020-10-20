
#include "Basics.h"
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"
#include "TGraph.h"
#include "TString.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
//#include "TCanvas.h"

CATS* KITTY_CATS_FIT_PL;

double Basics_Potential_Usmani(double* pars){
    double r = pars[0];
    //pars[1] is the momentum
    //Values for the potential
    double vbar = 6.2;//6.2
    double vsigma = 0.25;//0.25

    double x=r*0.7;
    double vc = pars[3]/(1+exp((r-pars[4])/pars[5]));
    double tpi = (1.0+3.0/x+3.0/(x*x)) * (exp(-x)/x) * pow(1.-exp(-2.*r*r),2.);

    double v = 0.;

    double& Spin = pars[2];
    if (Spin == 0) v = vc - (vbar + 0.75*vsigma)*tpi*tpi;//Usmani singlet
    else if (Spin == 1)  v = vc - (vbar - 0.25*vsigma)*tpi*tpi;//Usmani triplet
    else printf ("wrong polarization\n");

    return v;
}
double Basics_Source_Gauss(double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];
    double& Size = Pars[3];
    return 4.*Pi*Radius*Radius*pow(4.*Pi*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}
double Basics_PionCoulombCorrection(const bool& Identical){

}

//computes the expectation based on quantum statistics only
//the formula used is from chapter 4.1 in Phys. Rev. C 96 (2017) 064908 (ATLAS paper on pipi correlations in p-Pb)
TGraph* Basics_PiPiTheory(const bool& Identical, const bool& WithCoulomb){
    //load from a Mathematica output file
    FILE *InFileCBE;
    const TString CBEname = "../Mathematica/tab_txCBE.dat"; // changed to match actual exec. location
    InFileCBE = fopen(CBEname.Data(), "r");
    if(!InFileCBE){
        printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", CBEname.Data());
        return NULL;
    }
    char*  cdummy = new char [512];
    float kstar;
    float ckval;
    TGraph gCBE;
    gCBE.SetName("gCBE");
    unsigned NumPointsCBE=0;
    while(!feof(InFileCBE)){
        if(!fgets(cdummy, 511, InFileCBE)) continue;
        sscanf(cdummy, "%f %f",&kstar, &ckval);
        gCBE.SetPoint(NumPointsCBE,kstar,ckval);
        NumPointsCBE++;
    }
    fclose(InFileCBE);

    FILE *InFileK;
    const TString Kname = "../Mathematica/tab_txCoulombSame.dat";
    InFileK = fopen(Kname.Data(), "r");
    TGraph gK;
    gK.SetName("gK");
    unsigned NumPointsK=0;
    if(!InFileK){
        printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", Kname.Data());
        return NULL;
    }
    while(!feof(InFileK)){
        if(!fgets(cdummy, 511, InFileK)) continue;
        sscanf(cdummy, "%f %f",&kstar, &ckval);
        gK.SetPoint(NumPointsK,kstar,ckval);
        NumPointsK++;
    }
    fclose(InFileK);

    TGraph* gCk = new TGraph();
    gCk->SetName("gPiPiTheory");
    double dkstar,dckval;
    for(unsigned uBin=0; uBin<NumPointsCBE; uBin++){
        gCBE.GetPoint(uBin,dkstar,dckval);
        gCk->SetPoint(uBin,dkstar,dckval*gK.Eval(dkstar));
    }

    delete [] cdummy;
    return gCk;
}

TGraph* Basics_PiPiCATS(const bool& Identical, const bool& WithCoulomb){

    const unsigned NumMomBins = 200;
    const double kMin = 0;
    const double kMax = 100;

    CATS PionKitty;
    //(#bins,min,max) or (#bins,BinRangeArray[] as in ROOT)
    PionKitty.SetMomBins(NumMomBins,kMin,kMax);

    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource,1,true);
    //set the first and only par (source size)
    SOURCE_PARS.SetParameter(0,1.2);
    //say to CATS which Source function to use, and with which parameter set
    PionKitty.SetAnaSource(Basics_Source_Gauss,SOURCE_PARS);
    PionKitty.SetUseAnalyticSource(true);

    //standard settings for a CATS object which has no strong interaction potential included
    PionKitty.SetNumChannels(1);
    //#which channel, how many PWs
    PionKitty.SetNumPW(0,0);
    //which channel, spin value
    PionKitty.SetSpin(0,0);
    //which channel, weight
    PionKitty.SetChannelWeight(0, 1);

    //include the coulomb interaction. Q1Q2 is the multiplied charge numbers of the two particles
    PionKitty.SetQ1Q2(1);
    //the PdgId is needed when using a source from a transport model
    //actually the only important thing here is if the particles are identical
    if(Identical) PionKitty.SetPdgId(211, 211);
    else PionKitty.SetPdgId(-211, 211);

    if(WithCoulomb) PionKitty.SetRedMass( 0.5*Mass_Pic );
    else if(Identical) PionKitty.SetRedMass( 0.5*Mass_Pi0 );
    else PionKitty.SetRedMass( (Mass_Pi0*Mass_Pic)/(Mass_Pi0+Mass_Pic) );

    PionKitty.KillTheCat();

    TGraph* grCk = new TGraph();
    grCk->SetName("grCk");
    grCk->Set(NumMomBins);
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        grCk->SetPoint(uMom,PionKitty.GetMomentum(uMom),PionKitty.GetCorrFun(uMom));
    }

    return grCk;
}


double CATS_FIT_PL(double* x, double* pars){
   	// printf("CATS_FIT_PL has been called\n");
	double& Norm = pars[0];
    double& LambdaPar = pars[1];
    double& SourceSize = pars[2];
    //set the radius to the fit value
    //the last parameter says its a small change of the radius, which does not require a new computing grid (saves time)
    //however, the last step requires making sure a good initial value of the radius
    KITTY_CATS_FIT_PL->SetAnaSource(0,SourceSize,true);
    //useful tip: this makes CATS to shut up and not flood your screen. Only errors will be displayed. Use nSilent to suppress even those
    KITTY_CATS_FIT_PL->SetNotifications(CATS::nError);
    KITTY_CATS_FIT_PL->KillTheCat();
    return Norm*(LambdaPar*KITTY_CATS_FIT_PL->EvalCorrFun(*x)+1.-LambdaPar);
}

void Basics_ProtonLambda(){

    const unsigned NumMomBins = 60;
    const double kMin = 0;
    const double kMax = 240;
	class TCanvas;
	Bool_t UseDummy = false;
	Bool_t UseCanvas = true;
	Double_t VarFactorLower = 0.5;
	Double_t VarFactorUpper = 2.0;
	TString InputFileName = "OutputFiles/CSV3_OutputPL01.root";
	TString InputHistoName = "OutputHisto";
	TString OutputFileName = "OutputFilePLFit.root"; // TODO: use TString::Replace von InputFile
	if (UseDummy) {
		TString InputFileName = "OutputFiles/DummyProtonLambda.root";
		TString InputHistoName = "hDummyProtonLambda";
	}
    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource,1,true);
    SOURCE_PARS.SetParameter(0,1.3);

    CATS Kitty_pL;
    Kitty_pL.SetMomBins(NumMomBins,kMin,kMax);
    Kitty_pL.SetAnaSource(Basics_Source_Gauss,SOURCE_PARS);
    Kitty_pL.SetAnaSource(0,SOURCE_PARS.GetParameter(0));
    Kitty_pL.SetUseAnalyticSource(true);
    Kitty_pL.SetMomentumDependentSource(false);
    Kitty_pL.SetThetaDependentSource(false);
    //should you include in the result any bins, where the Schroedinger solver failed
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
    //WhichChannel,WhichPW,PotentialFunction,Parameters
    Kitty_pL.SetShortRangePotential(0,0,Basics_Potential_Usmani,POT_PARS_1S0);
    Kitty_pL.SetShortRangePotential(1,0,Basics_Potential_Usmani,POT_PARS_3S1);
    Kitty_pL.SetMaxNumThreads(4);
    Kitty_pL.KillTheCat();

    KITTY_CATS_FIT_PL = &Kitty_pL;

    TFile* fInput = new TFile(InputFileName,"read");	// replaced "Files" with "OutputFiles"
    TH1F* hInput = (TH1F*)fInput->Get(InputHistoName);
	//
	//adding section for the GeV->MeV rescaling:
	printf("hInput range: %f - %f\n", hInput->GetXaxis()->GetXmin(), hInput->GetXaxis()->GetXmax());
	printf("hInput Nbins: %d\n", hInput->GetNbinsX());
	/*
	TH1F* hInputRescaled = new TH1F("hInputRescaled", "hInputRescaled", (hInput->GetNbinsX())/3, 0, 500);
    printf("Created Rescaled Histo\n");
    for(unsigned uBin=1; uBin<=(hInput->GetNbinsX())/3; uBin++) hInputRescaled->SetBinContent(uBin, hInput->GetBinContent(uBin));
    printf("Filled the Rescaled Histo\n");
    for(unsigned uBin=1; uBin<=(hInput->GetNbinsX())/3; uBin++) hInputRescaled->SetBinError(uBin, hInput->GetBinError(uBin));
	*/
	//

	TH1F* hInputRescaled = new TH1F("hInputRescaled", "hInputRescaled", (hInput->GetNbinsX())/12, 0, 250);
    printf("Created Rescaled Histo\n");
    for(unsigned uBin=1; uBin<=(hInput->GetNbinsX())/12; uBin++) hInputRescaled->SetBinContent(uBin, hInput->GetBinContent(uBin));
    printf("Filled the Rescaled Histo\n");
    for(unsigned uBin=1; uBin<=(hInput->GetNbinsX())/12; uBin++) hInputRescaled->SetBinError(uBin, hInput->GetBinError(uBin));
	//

    TF1* FITTER = new TF1("FITTER",CATS_FIT_PL,kMin,kMax,6);
    FITTER->SetParameter(0,1);FITTER->SetParLimits(0,0.5,2);
    FITTER->SetParameter(1,0.5);FITTER->SetParLimits(1,0,1);
    FITTER->SetParameter(2,SOURCE_PARS.GetParameter(0));FITTER->SetParLimits(2,0.5,3);
    
	FITTER->SetParameter(3,POT_PARS_1S0.GetParameter(1));FITTER->SetParLimits(3,VarFactorLower*POT_PARS_1S0.GetParameter(1),VarFactorUpper*POT_PARS_1S0.GetParameter(1)); // originally at 0.99 and 1.01 respectively -> made wider
    //FITTER->SetParameter(4,POT_PARS_1S0.GetParameter(2));FITTER->SetParLimits(4,VarFactorLower*POT_PARS_1S0.GetParameter(2),VarFactorUpper*POT_PARS_1S0.GetParameter(2));
    //FITTER->SetParameter(5,POT_PARS_1S0.GetParameter(3));FITTER->SetParLimits(5,VarFactorLower*POT_PARS_1S0.GetParameter(3),VarFactorUpper*POT_PARS_1S0.GetParameter(3));

    //FITTER->FixParameter(1,0.6);

    //FITTER->FixParameter(3,POT_PARS_1S0.GetParameter(1));
    FITTER->FixParameter(4,POT_PARS_1S0.GetParameter(2));
    FITTER->FixParameter(5,POT_PARS_1S0.GetParameter(3));


    //hInput->Fit(FITTER,"S, N, R, M");
    hInputRescaled->Fit(FITTER,"S, N, R, M");
    //hInputRescaled->Fit(FITTER,"R, M");
	//Descriptions of Options from: https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html
	//"S": The result of the fit is returned in the TFitResultPtr
	//"N": Do not store the graphics function, donot draw
	//"R": Use the range specified in the function range
	//"M": Improve fit results, by using the IMPROVE algorithm of TMinuit
	
	TFile* fOut = TFile::Open(OutputFileName, "RECREATE");
	fOut->cd();
	hInput->Write();
	hInputRescaled->Write();
	FITTER->Write();
	fOut->Close();


	delete fOut;
    delete FITTER;
    delete fInput;
}
