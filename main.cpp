#include<iostream>
#include "Basics.h"
#include "ExtendedCk.h"

#include "TROOT.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"

void ComparePionPion(){
    TGraph* grCATS = Basics_PiPiCATS(1,1);
    TGraph* grTheory = Basics_PiPiTheory(1,1);
    TFile* OutputFile = new TFile("ComparePionPion.root","recreate");
    grCATS->Write();
    grTheory->Write();

    delete grCATS;
    delete grTheory;
    delete OutputFile;
}

int main(int argc, char *argv[]){

    printf("\nHello to this lovely CATS tutorial!\n");
    printf(" To find a bug: continue with this tutorial\n");
    printf("-------------------------------------------\n");

    //ComparePionPion();
    //printf("Done with ComparePionPion\n");

    //Basics_ProtonLambda();
	//printf("Done with Basics_ProtonLambda\n");

    Basics_LambdaLambda();
	printf("Done with Basics_LambdaLambda\n");
 
	//Check_Potential();
	//printf("Done with Check_Potential");
	
	//Ck_pL_Ledni_Usmani();
	//printf("Done with Ck_pL_Ledni_Usmani\n");
 
    //Ck_pp_Decomposition("Gauss");
    //Ck_pp_Decomposition("CoreReso");

    return 0;
}
