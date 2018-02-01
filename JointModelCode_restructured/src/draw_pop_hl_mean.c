/*******************************************************************/
/*************************draw_bh_mean_trigger.c *****************************/
/*******************************************************************/

/**fix here***/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"



/*********************************************************************/
/*START OF draw_bh_mean SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_bh_mean: this runs the Gibbs sampler draw for mean baseline and half life;
	ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
	Priors *priors; the current values of the prior parameters;
	Common_parms *parms; the current values of the common parameters;
	unsigned long *seed; seed values needed for the randon number generator;
	Hyper_priors *hyper; parameters of the prior distributions;
	RETURNS: None; all updates are made internally
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
 fesum: sum of subject specific mean baselines and half-lifes;
 part of evaluation of gmean and gvar
 gmean: mean of the distribution used in Gibbs sampler
 gvar: variance of distribution used in Gibbs sampler
 *subject: current list of subjects and their characteristics

 SUBROUTINES USED
 rnorm: found in randgen.h; draws from the normal distribution
 kiss: found in randgen.h; draws from U(0,1) distribution
**************************************************************************/

void draw_bh_mean_trigger(Subject_type *sublist, Priors *priors, Common_parms *parms, unsigned long *seed, Hyper_priors *hyper)
{
	/* declare variables */
	int j;
	double fesum, gmean, gvar;
	Subject_type *subject;

	/*declare functions */
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);

	for (j = 0; j < 2; j++){

		fesum = 0;
		subject = sublist->succ;

		/*This while loop goes through the pulses to get the sum involved*/
		while (subject != NULL){
			fesum += subject->basehalf_l[j];
			subject = subject->succ;
		} /*end of loop through pulses*/

		gmean = (hyper->meanvarbh_l[j] * fesum) + (priors->varbh_l[j] * priors->varbh_l[j] * hyper->meanmeanbh_l[j]);
		gmean /= (parms->numsub*hyper->meanvarbh_l[j]) + (priors->varbh_l[j] * priors->varbh_l[j]);

		gvar = hyper->meanvarbh_l[j] * priors->varbh_l[j] * priors->varbh_l[j];
		gvar /= (parms->numsub*hyper->meanvarbh_l[j]) + (priors->varbh_l[j] * priors->varbh_l[j]);

		priors->meanbh_l[j] = rnorm(gmean, sqrt(gvar), seed);
	}

}
