#include "TOV_library.h"
#include "read_par_eos.h"
#include "read_par_tov.h"

int main(){

    int usePP=0, useTab=0, useTASPP=0;

    struct PP eosPP; struct TASPP eosTASPP; struct Tab eosTab;
    read_par_eos(&eosPP, &eosTASPP, &eosTab, "eos.par", &usePP, &useTab, &useTASPP);
    int saveRadial, checkAccuracy, sequence; double *seqdetails;
    read_par_tov("tov.par", &saveRadial, &checkAccuracy, &sequence, &seqdetails);
    double Mmax, Lambda, error, radius;
    if(sequence){
        if(usePP)
            starSequence_PP(seqdetails[0], seqdetails[1], seqdetails[2], checkAccuracy, eosPP);
        if(useTASPP)
            starSequence_TASPP(seqdetails[0], seqdetails[1], seqdetails[2], checkAccuracy, eosTASPP);
        if(useTab)
            starSequence_TAB(seqdetails[0], seqdetails[1], seqdetails[2], checkAccuracy, eosTab);
    }

    else{
        if(usePP)
            radius=solveStarMRPL_PP(&Mmax, &Lambda, seqdetails[0], saveRadial, checkAccuracy, &error, eosPP);
        if(useTASPP)
            radius=solveStarMRPL_TASPP(&Mmax, &Lambda, seqdetails[0], saveRadial, checkAccuracy, &error, eosTASPP);
        if(useTab)
            radius=solveStarMRPL_TAB(&Mmax, &Lambda, seqdetails[0], saveRadial, checkAccuracy, &error, eosTab);
        printf("The NS calculated has:\nradius(km):%.15e \nmass(Ms):%.15e \nLambda():%.15e \n",radius, Mmax, Lambda);
        if(checkAccuracy) printf("The maximum relative difference from constant is: %.15e\n", error);
    }




    return 0;
}