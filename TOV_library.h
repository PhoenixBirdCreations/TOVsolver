#include "EOS_cold.h"

#define SIGN(a) (a)>0 ? 1 : -1

#define Ggrav 6.67408e-8
#define c2 898755178736817640000.0   //cgs units

double solveStarMRP_PP(double *Mmax, double rho_c, int saveRadialProfile, int checkAccuracyMetric, double *Maxerror, struct PP eos);
double solveStarMRPL_PP(double *Mmax, double *LambdaNS, double rho_c, int saveRadialProfile, int checkAccuracyMetric, double *Maxerror,struct PP eos);

double solveStarMRP_TASPP(double *Mmax, double rho_c, int saveRadialProfile, int checkAccuracyMetric, double *Maxerror,struct TASPP eos);
double solveStarMRPL_TASPP(double *Mmax, double *LambdaNS, double rho_c, int saveRadialProfile, int checkAccuracyMetric, double *Maxerror,struct TASPP eos);

double solveStarMRP_TAB(double *Mmax, double rho_c, int saveRadialProfile, int checkAccuracyMetric, double *Maxerror,struct Tab eos);
double solveStarMRPL_TAB(double *Mmax, double *LambdaNS, double rho_c, int saveRadialProfile, int checkAccuracyMetric, double *Maxerror, struct Tab eos);

void starSequence_PP(double rho_ini, double rho_end, double rho_step, int checkAccuracy, struct PP eos);
void starSequence_TASPP(double rho_ini, double rho_end, double rho_step, int checkAccuracy, struct TASPP eos);
void starSequence_TAB(double rho_ini, double rho_end, double rho_step, int checkAccuracy, struct Tab eos);