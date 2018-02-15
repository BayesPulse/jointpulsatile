
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"

void mcmc(Subject_type *sublist,Common_parms *parms_l,Common_parms *parms_f,double **ts_l,double **ts_f,long iter,int N,
          Priors *priors,unsigned long *seed,char *filel,char *file2,double propvar[],Hyper_priors *hyper);

draw_pop_mass_mean(Patient *, PopulationEstimate *, PopulationEstimate *, PopulationPrior *, PopulationProposal *, PopulationPrior *, PopulationProposal *, AssocEstimates *, PopulationProposal *, int , unsigned long *);
draw_pop_width_mean(Patient *, PopulationEstimate *, PopulationPrior *, PopulationProposal *, int , unsigned long *);

double Phi(double y);
void adjust_acceptance(double x,double *X);
void adjust2_acceptance(double x,double **X,double corr);
double gamm(double x) ;

double gammaPdf(double x, double a, double b);
double phi(double, double, double);



