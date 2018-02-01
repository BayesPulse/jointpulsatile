/*******************************************************************/
/*************************draw_random_effects_trigger.c *****************************/
/*******************************************************************/

/**8fix here***/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"


/*********************************************************************/
/*START OF draw_random_effects SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_random_effects: this runs the M-H draw for individual pulse masses and
					  widths;
					  ARGUMENTS: double **ts; this is the matrix of observed data (a column of
					  times and a column of log(concentration);
					  Subject_type *list; this is the current list of subjects that exist;
					  Common_parms *parms; the current values of the common parameters;
					  int N; the number of observations in **ts;
					  double v1; the proposal variance for individual pulse masses;
					  double v2; the proposal variance for individual pulse widths;
					  unsigned long *seed; seed values needed for the randon number generator;
					  RETURNS: None; all updates are made internally; this does update the global
					  variables regarding acceptance rate: arem, arew, nrem, nrew
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j,k: generic counters
 temp: value used in evaluation of log(rho)
 prior_old: part of evaluation of prior_ratio portion of log(rho)
 prior_new: same as prior_old, but with proposed value
 old_val: current value of the random effect is saved to this variable so we
 can evaluate likelihood under the proposed value
 *tmp1: current acceptance counters are saved to this vector so we can increment
 inside a loop
 *pRE: proposed values of individual mass and width
 prior_ratio: comparison between prior distributions under current and proposed values
 like_ratio: difference between current and proposed log-likelihoods
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 plikelihood: likelihood under proposed value
 *old_contrib: current mean_contrib for a pulse is saved to this vector
 so we can evaluate likelihood under proposed pulse's mean_contrib
 current_like: current likelihood for a given subject; used for comparison
 *subnode: for a given subject, current list of pulses and their qualities
 *subject: current list of subjects and their qualities

 SUBROUTINES USED
 mean_contribution: found in birthdeath.c; evaluates mean_contrib for a
 pulse with specific parameters
 kiss: found in randgen.h; draws from U(0,1) distribution
 rnorm: found in randgen.h; draws from the normal distribution
 likelihood2: found in birthdeath.c; calculates and returns likelihood for a
 given subject under current parameters
**************************************************************************/

void draw_random_effects_trigger(double **ts, Subject_type *sublist, Common_parms *parms, int N, double v1, double v2, unsigned long *seed)
{
	int *tmp1, i, j, k;
	double temp, prior_old, prior_new, old_val, *pRE, prior_ratio, like_ratio;
	double alpha, plikelihood, *old_contrib, current_like, stdold[2], stdnew[2], newint[2], oldint[2];
	Node_type *subnode;
	Subject_type *subject;
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);
	double kiss(unsigned long *);
	double rnorm(double, double, unsigned long *);
	double likelihood2(Node_type *, double **, Common_parms *, int, Node_type *, double);
	double Phi(double);

	/*Set acceptance counts equal to temporary vector*/
	tmp1 = (int *)calloc(2, sizeof(int));
	tmp1[0] = arem_l;
	tmp1[1] = arew_l;

	/*Go to first subject*/
	subject = sublist->succ;

	/*Index of subjects used for likelihood2*/
	parms->subindex = 0;

	/*Loop thru subjects*/
	while (subject != NULL){

		/*Go to start of pulses for this subject*/
		subnode = subject->driver->succ;
		k = 0;

		/*Go through each existing pulse for this subject*/
		while (subnode != NULL) {

			/*Allocate Memory*/
			old_contrib = (double *)calloc(N, sizeof(double));
			pRE = (double *)calloc(2, sizeof(double));

			/*Increase the denominators of the acceptance rates*/
			nrem_l++;
			nrew_l++;

			/*Draw proposed values of current pulse's mass and width*/
			pRE[0] = rnorm(subnode->theta[0], v1, seed);
			pRE[1] = rnorm(subnode->theta[1], v2, seed);
			/*only proceed when we draw positive value*/
			if (pRE[0] > 0.0 && pRE[1] > 0.01 && pRE[1] < 100){

				/*Determine if we accept or reject proposed pulse mass then determine*/
				/*if we accept or reject proposed pulse width*/
				for (j = 0; j < 2; j++){
					current_like = likelihood2(subject->driver, ts, parms, N, subject->driver, subject->basehalf_l[0]);

					/*Compute the log of the ratio of the priors*/
					prior_old = subnode->theta[j] - subject->theta_l[j];
					prior_old *= 0.5*prior_old;
					prior_new = pRE[j] - subject->theta_l[j];
					prior_new *= 0.5*prior_new;
					prior_ratio = (prior_old - prior_new)*subnode->eta[j];
					prior_ratio /= parms->re_precision[j];
					prior_ratio /= parms->re_precision[j];
					stdold[j] = (subnode->theta[j]) / (parms->re_precision[j] / sqrt(subnode->eta[j]));
					stdnew[j] = pRE[j] / (parms->re_precision[j] / sqrt(subnode->eta[j]));

					newint[j] = log(Phi(stdnew[j]));
					oldint[j] = log(Phi(stdold[j]));
					/*prior_ratio = prior_ratio-newint[j]+oldint[j];*/
					/*Save the current value of mass/width*/
					old_val = subnode->theta[j];

					/*Set the pulse's mass/width equal to the proposed value*/
					subnode->theta[j] = pRE[j];

					/*Save the mean_contrib for that pulse*/
					for (i = 0; i < N; i++)
						old_contrib[i] = subnode->mean_contrib[i];

					/*Recalculate that pulse's mean_contrib assuming proposed mass/width */
					mean_contribution(subnode, ts, parms, N, subject->basehalf_l[1]);

					/*Calculate likelihood assuming proposed mass/width */
					plikelihood = likelihood2(subject->driver, ts, parms, N, subject->driver, subject->basehalf_l[0]);

					like_ratio = plikelihood - current_like;

					/*Compute the log of the ratio between the two likelihoods*/

					/*Calculate log rho; set alpha equal to min(0,log rho) */
					alpha = (0 < (temp = (prior_ratio + like_ratio))) ? 0 : temp;

					/*If log U < log rho, accept the proposed value, increase acceptance*/
					/*counter */
					if (log(kiss(seed)) < alpha) {
						tmp1[j]++;
					}

					/*Otherwise, reject the proposed value, set pulse's mass/width back to */
					/*saved current value, and set pulse's mean_contrib back to saved value */
					else {
						subnode->theta[j] = old_val;
						for (i = 0; i < N; i++)
							subnode->mean_contrib[i] = old_contrib[i];
					}
				}
			}		 /*end of loop through mass & width */

			/*Advance to next pulse*/
			subnode = subnode->succ;
			k++;
			/*free memory*/
			free(pRE);
			free(old_contrib);

		} /*end of loop through pulses*/

		/*advance to next subject*/
		subject = subject->succ;
		parms->subindex++;

	} /*end of loop through subjects*/
	/*Set counter equal to temporary vector values*/
	arem_l = tmp1[0];
	arew_l = tmp1[1];
	free(tmp1);

}
