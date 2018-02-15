/*******************************************************************/
/*************************draw_fe_priors_w_var.c *****************************/
/*******************************************************************/

/**fix here**/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"


/*********************************************************************/
/*START OF draw_fe_priors2 SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_fe_priors2: this runs the M-H draw for standard deviations of subject
				   level means of mass and width;
				   ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
				   Priors *priors; the current values of the prior parameters;
				   double v1; the proposal variance for overall st-dev of mean mass;
				   double v2; the proposal variance for overall st-dev of mean width;
				   unsigned long *seed; seed values needed for the randon number generator;
				   Hyper_priors *hyper; parameters of the prior distributions;
				   RETURNS: None; all updates are made internally; this does update the global
				   variables regarding acceptance rate: afepmv, afepwv, nfepmv, and nfepwv
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j,t: generic counters
 *tmp1: current acceptance counters are saved to this vector so we can increment
 inside a loop
 *new_var: proposed values of st. dev. of mass and width
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

void draw_fe_priors_w_var(Subject_type *sublist, Priors *priors, double v1, double v2, unsigned long *seed, Hyper_priors *hyper)
{
	int j, t, *tmp1;
	double *new_var, prop_new[2], prop_old[2], psum[2], pcomp[2], prior_old, prior_new;
	double prior_ratio[2], prop_ratio[2], temp, alpha;
	Subject_type *subject;
	double kiss(unsigned long *);
	double rnorm(double, double, unsigned long *);

	/*Allocate Memory*/
	new_var = (double *)calloc(2, sizeof(double));
	tmp1 = (int *)calloc(2, sizeof(int));

	/*Add 1 to the counters for acceptance rates of sigma_ma and sigma_mw*/
	nfewv_l++;
	nfewv_f++;

	/*Assign current acceptance counts to temporary vector*/
	tmp1[0] = afewv_l;
	tmp1[1] = afewv_f;

	/*Draw proposed values for sigma_ma and sigma_mw*/
	new_var[0] = rnorm(priors->fe_precision_wl, v1, seed);
	new_var[1] = rnorm(priors->fe_precision_wf, v2, seed);

	/*Accept or Reject sigma_ma, then accept or reject for sigma_mw*/

	/*t and j are both used to test what happens if we reverse the
	order of sma and smw*/

	/*Calculate ratio of Cauchy priors*/
	/*      prior_new = new_var[j]/hyper->prec[j];
		  prior_old = priors->fe_precision[j]/hyper->prec[j];
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
	psum[0] = psum[1] = 0;
	subject = sublist->succ;
	prior_ratio[0] = prior_ratio[1] = 1;
	while (subject != NULL){
		prior_ratio[0] *= priors->fe_precision_wl / new_var[0];
		pcomp[0] = subject->theta_l[1] - priors->fe_mean_l[1];
		pcomp[0] *= pcomp[0];
		psum[0] += pcomp[0];

		prior_ratio[1] *= priors->fe_precision_wf / new_var[1];
		pcomp[1] = subject->theta_f[1] - priors->fe_mean_f[1];
		pcomp[1] *= pcomp[1];
		psum[1] += pcomp[1];

		subject = subject->succ;
	}

	/*Log prior ratio value*/
	prior_ratio[0] = log(prior_ratio[0]);
	prior_ratio[1] = log(prior_ratio[1]);


	/*Complete the "likelihood" portion of rho*/
	prop_old[1] = 1.0 / priors->fe_precision_wf;
	prop_old[1] /= priors->fe_precision_wf;
	prop_old[0] = 1.0 / priors->fe_precision_wl;
	prop_old[0] /= priors->fe_precision_wl;


	prop_new[0] = 1.0 / new_var[0];
	prop_new[0] /= new_var[0];

	prop_new[1] = 1.0 / new_var[1];
	prop_new[1] /= new_var[1];
	prop_ratio[0] = psum[0] * 0.5*(prop_old[0] - prop_new[0]);
	prop_ratio[1] = psum[1] * 0.5*(prop_old[1] - prop_new[1]);


	/*We only can accept the proposed value if it is positive*/
	if (new_var[0] > 0 && new_var[0] < hyper->precision_wl){

		/*Compute log rho, and set alpha equal to min(log rho,0)*/
		alpha = (0 < (temp = (prior_ratio[0] + prop_ratio[0]))) ? 0 : temp;

		/*If log(U) < log rho, accept the proposed value*/
		/*Increase acceptance count by 1*/
		if (log(kiss(seed))<alpha){
			tmp1[0]++;
			priors->fe_precision_wl = new_var[0];
		}
	}
	if (new_var[1]>0 && new_var[1] < hyper->precision_wf){

		/*Compute log rho, and set alpha equal to min(log rho,0)*/
		alpha = (0 < (temp = (prior_ratio[1] + prop_ratio[1]))) ? 0 : temp;

		/*If log(U) < log rho, accept the proposed value*/
		/*Increase acceptance count by 1*/
		if (log(kiss(seed)) < alpha){
			tmp1[1]++;
			priors->fe_precision_wf = new_var[1];
		}
	}

	/*Set acceptance count equal to temp vector components*/
	afewv_l = tmp1[0];
	afewv_f = tmp1[1];

	free(new_var);
	free(tmp1);

}
