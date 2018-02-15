/*******************************************************************/
/*************************draw_pop_width_mean.c*********************/
/*******************************************************************/

#include "draw_pop_width_mean.h"

void draw_pop_widht_mean(Patient *patientlist, PopulationPriors *priors, Common_parms *parms_f, unsigned long *seed, Hyper_priors *hyper)
{
	/* declare variables */
	int j;
	double fesum[2], gmean[2], gvar[2];
	Subject_type *subject;

	/*declare functions */
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);



	fesum[0] = fesum[1] = 0.0;
	subject = sublist->succ;

	/*This while loop goes through the pulses to get the sum involved*/
	while (subject != NULL){
		fesum[0] += subject->theta_l[1];
		fesum[1] += subject->theta_f[1];

		subject = subject->succ;
	} /*end of loop through pulses*/

	/*Gibbs sampler mean*/

	gmean[0] = (hyper->hvar_l*fesum[0]) + (priors->fe_precision_wl*priors->fe_precision_wl*hyper->hmean_l[1]);
	gmean[0] /= (parms_f->numsub*hyper->hvar_l) + (priors->fe_precision_wl*priors->fe_precision_wl);


	gmean[1] = (hyper->hvar_f*fesum[1]) + (priors->fe_precision_wf*priors->fe_precision_wf*hyper->hmean_f[1]);
	gmean[1] /= (parms_f->numsub*hyper->hvar_f) + (priors->fe_precision_wf*priors->fe_precision_wf);

	/*Gibbs sampler variance*/
	gvar[0] = hyper->hvar_l*priors->fe_precision_wl*priors->fe_precision_wl;
	gvar[0] /= (parms_f->numsub*hyper->hvar_l) + (priors->fe_precision_wl*priors->fe_precision_wl);
	gvar[1] = hyper->hvar_f*priors->fe_precision_wf*priors->fe_precision_wf;
	gvar[1] /= (parms_f->numsub*hyper->hvar_f) + (priors->fe_precision_wf*priors->fe_precision_wf);


	/*This draws the new value for mua/muw using a Gibbs sampler*/
	priors->fe_mean_l[1] = rnorm(gmean[0], sqrt(gvar[0]), seed);
	priors->fe_mean_f[1] = rnorm(gmean[1], sqrt(gvar[1]), seed);



}
