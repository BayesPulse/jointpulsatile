/*******************************************************************/
/*************************draw_bh_var_trigger.c ********************/
/*******************************************************************/

/***fix here***/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"


/*********************************************************************/
/*START OF draw_bh_var SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_bh_var: this runs the M-H draw for standard deviations of subject
				   level baselines and half-lifes;
				   ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
				   Priors *priors; the current values of the prior parameters;
				   double v1; the proposal variance for overall st-dev of baseline;
				   double v2; the proposal variance for overall st-dev of half-life;
				   unsigned long *seed; seed values needed for the randon number generator;
				   Hyper_priors *hyper; parameters of the prior distributions;
				   RETURNS: None; all updates are made internally; this does update the global
				   variables regarding acceptance rate: afebv, afehv, nfebv, and nfehv
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
 *tmp1: current acceptance counters are saved to this vector so we can increment
 inside a loop
 *new_var: proposed values of st. dev. of baseline and half-life
 prop_new: part of evaluation of prop_ratio portion of log(rho)
 prop_old: same as prop_new but with current (not proposed) value
 psum: part of evaluation of prop_ratio portion of log(rho)
 pcomp: part of evaluation of prop_ratio portion of log(rho)
 prior_old: part of evaluation of prior_ratio portion of log(rho)
 prior_new: same as prior_old, but with proposed value
 prior_ratio: comparison between prior distributions under current and proposed values
 prop_ratio: comparison between random effects distributions under current and
 proposed values
 temp: value used in evaluation of log(rho)
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 *subject: current list of subjects and their qualities

 SUBROUTINES USED
 kiss: found in randgen.h; draws from U(0,1) distribution
 rnorm: found in randgen.h; draws from the normal distribution
**************************************************************************/

void draw_bh_var_trigger(Subject_type *sublist, Priors *priors, double v1, double v2, unsigned long *seed, Hyper_priors *hyper)
{
	int j, *tmp1;
	double *new_var, prop_new, prop_old, psum, pcomp, prior_old, prior_new;
	double prior_ratio, prop_ratio, temp, alpha;
	Subject_type *subject;
	double kiss(unsigned long *);
	double rnorm(double, double, unsigned long *);

	/*Allocate Memory*/
	new_var = (double *)calloc(2, sizeof(double));
	tmp1 = (int *)calloc(2, sizeof(int));

	/*Add 1 to the counters for acceptance rates of sigma_b and sigma_h*/
	nfebv_l++;
	nfehv_l++;

	/*Assign current acceptance counts to temporary vector*/
	tmp1[0] = afebv_l;
	tmp1[1] = afehv_l;

	/*Draw proposed values for sigma_b and sigma_h*/
	new_var[0] = rnorm(priors->varbh_l[0], v1, seed);
	new_var[1] = rnorm(priors->varbh_l[1], v2, seed);

	/*Accept or Reject sigma_b, then accept or reject for sigma_h*/
	for (j = 0; j < 2; j++){
		/*Calculate ratio of Cauchy priors*/
		/*prior_new = new_var[j]/hyper->varbh[j];
		prior_old = priors->varbh[j]/hyper->varbh[j];
		prior_new *= prior_new;
		prior_old *= prior_old;
		prior_new += 1.0;
		prior_old += 1.0;
		prior_ratio = prior_old/prior_new;
		*/
		/*Compute the sum included in the "likelihood"*/
		/*Also compute the old value divided by the new value raised to num_node*/
		/*This sum is the same assuming the current value of sigma_j in both*/
		/*the numerator and denominator of rho*/
		psum = 0;
		subject = sublist->succ;
		prior_ratio = 1;
		while (subject != NULL){
			prior_ratio *= priors->varbh_l[j] / new_var[j];
			pcomp = subject->basehalf_l[j] - priors->meanbh_l[j];
			pcomp *= pcomp;
			psum += pcomp;
			subject = subject->succ;
		}

		/*Log prior ratio value*/
		prior_ratio = log(prior_ratio);

		/*Complete the "likelihood" portion of rho*/
		prop_old = 1.0 / priors->varbh_l[j];
		prop_old /= priors->varbh_l[j];
		prop_new = 1.0 / new_var[j];
		prop_new /= new_var[j];
		prop_ratio = psum*0.5*(prop_old - prop_new);

		/*We only can accept the proposed value if it is positive*/
		/*With a uniform(0,2) prior, skip this step if we have a sd less than 2*/
		if (new_var[j] > 0 && new_var[j] < hyper->varbh_l[j]){

			/*Compute log rho, and set alpha equal to min(log rho,0)*/
			/*With a uniform prior, get rid of prior_ratio*/
			alpha = (0 < (temp = (prior_ratio + prop_ratio))) ? 0 : temp;

			/*If log(U) < log rho, accept the proposed value*/
			/*Increase acceptance count by 1*/
			if (log(kiss(seed)) < alpha){
				tmp1[j]++;
				priors->varbh_l[j] = new_var[j];
			}
		}
	}

	/*Set acceptance count equal to temp vector components*/
	afebv_l = tmp1[0];
	afehv_l = tmp1[1];

	free(new_var);
	free(tmp1);

}
