/*******************************************************************/
/*************************draw_bh_trigger.c ************************/
/*******************************************************************/

/***fix hwew***/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"

/*********************************************************************/
/*START OF draw_bh SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_bh: this runs the M-H draw for baseline and halflife (drawn together);
	ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
	Common_parms *parms; the current values of the common parameters;
	Priors *priors; the current values of the prior parameters;
	double **ts; this is the matrix of observed data (a column of
	times and S columns of log(concentration);
	int N; the number of observations in each column of **ts;
	unsigned long *seed; seed values needed for the randon number generator;
	double **var; the proposal variance-covariance matrix for
	baseline and halflife
	RETURNS: None; all updates are made internally; this does update the global
	variables regarding acceptance rate: adelta and ndelta
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j,k: generic counters
 num_node: number of pulses
 *subnode: for a given subject, current list of pulses and their qualities
 *new_node: counter for number of pulses
 *subject: current list of subjects and their characteristics
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 plikelihood: likelihood under proposed value
 *pmd: vector of proposed baseline and halflife
 like_ratio: difference between current and proposed log-likelihoods
 priorb_old: part of evaluation of prior_ratio portion of log(rho)
 priorb_new: same as priorb_old, but with proposed baseline
 priorh_old: part of evaluation of prior_ratio portion of log(rho)
 priorh_new: same as priorb_old, but with proposed halflife
 prior_ratio: comparison between prior distributions under current and proposed values
 currentmd[2]: current baseline and halflife saved to this variable name so we can
 evaluate likelihood under proposed baseline and halflife
 **currentmc: current mean_contrib for each pulse is saved to this matrix
 so we can evaluate likelihood under proposed pulses' mean_contrib
 temp: value used in evaluation of log(rho)
 current_decay: current decay rate is saved to this variable name so we can
 evaluate likelihood under proposed decay rate
 current_like: current likelihood for a given subject; saved so we can compare

 SUBROUTINES USED
 rmvnorm: found in randgen.h; draws from the multivariate normal distribution
 kiss: found in randgen.h; draws from U(0,1) distribh_f
 ution
 mean_contribution: found in birthdeath.c; evaluates mean_contrib for a
 pulse with specific parameters
 likelihood2: found in birthdeath.c; calculates and returns likelihood for a
 given subject under current parameters*/
/**************************************************************************/

void draw_bh_l(Subject_type *sublist, Common_parms *parms, Priors *priors, double **ts,
	int N, unsigned long *seed, double **var)
{
	int i, j, k, num_node;
	Node_type *subnode, *new_node;
	Subject_type *subject;
	double alpha, plikelihood, *pmd, like_ratio, priorb_old, priorb_new, priorh_old;
	double priorh_new, prior_ratio, currentmd[2], **currentmc, temp, current_decay;
	double current_like;
	int rmvnorm(double *, double **, int, double *, unsigned long *, int);
	double kiss(unsigned long *);
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);
	double likelihood2(Node_type *, double **, Common_parms *, int, Node_type *, double);

	/*Go to start of list of subjects*/
	subject = sublist->succ;

	/*this tells likelihood2 which subject's likelihood to calculate*/
	parms->subindex = 0;

	/*Loop thru all subjects*/
	while (subject != NULL){

		/*Calculate the current likelihood for this subject*/
		current_like = likelihood2(subject->driver, ts, parms, N, subject->driver, subject->basehalf_l[0]);

		/*Count number of pulses for this subject*/
		num_node = 0;
		new_node = subject->driver->succ;
		while (new_node != NULL) {
			new_node = new_node->succ;
			num_node++;
		}

		/*Allocate memory*/
		pmd = (double *)calloc(2, sizeof(double));
		currentmc = (double **)calloc(num_node, sizeof(double *));
		for (i = 0; i<num_node; i++)
			currentmc[i] = (double *)calloc(N, sizeof(double));

		/*Increase denominator of acceptance rate for b and hl*/
		ndelta_l++;

		/*Draw proposal values for b and hl*/
		rmvnorm(pmd, var, 2, subject->basehalf_l, seed, 1);

		/*Only proceed if we draw reasonable values*/
		if (pmd[0] > 0 && pmd[1] > 3) {

			/*Compute ratio of prior densities*/
			priorb_old = subject->basehalf_l[0] - priors->meanbh_l[0];
			priorb_old *= 0.5*priorb_old / (priors->varbh_l[0] * priors->varbh_l[0]);
			priorb_new = pmd[0] - priors->meanbh_l[0];
			priorb_new *= 0.5*priorb_new / (priors->varbh_l[0] * priors->varbh_l[0]);
			priorh_old = subject->basehalf_l[1] - priors->meanbh_l[1];
			priorh_old *= 0.5*priorh_old / (priors->varbh_l[1] * priors->varbh_l[1]);
			priorh_new = pmd[1] - priors->meanbh_l[1];
			priorh_new *= 0.5*priorh_new / (priors->varbh_l[1] * priors->varbh_l[1]);

			prior_ratio = priorb_old + priorh_old - priorb_new - priorh_new;

			/*Save current values of b and hl and set new current values equal to
			proposal values */
			for (k = 0; k < 2; k++) {
				currentmd[k] = subject->basehalf_l[k];
				subject->basehalf_l[k] = pmd[k];
			}

			/*Save current decay rate; calculate new current decay rate based on
			  proposal value of hl */
			current_decay = subject->decay_l;
			subject->decay_l = log(2) / subject->basehalf_l[1];

			/*Save current mean_contrib for each pulse; calculate new current
			  mean_contrib for each pulse based on proposed values of b and hl*/
			i = 0;
			subnode = subject->driver->succ;
			while (subnode != NULL) {
				for (j = 0; j < N; j++)
					currentmc[i][j] = subnode->mean_contrib[j];
				mean_contribution(subnode, ts, parms, N, subject->basehalf_l[1]);
				i++;
				subnode = subnode->succ;
			}

			/*Calculate proposed likelihood and then calculate likelihood ratio*/
			plikelihood = likelihood2(subject->driver, ts, parms, N, subject->driver, subject->basehalf_l[0]);
			like_ratio = plikelihood - current_like;

			/*Calculate log rho; set alpha equal to min(0,log rho) */
			alpha = (0 < (temp = (prior_ratio + like_ratio))) ? 0 : temp;

			/*If log U < log rho, increase acceptance rate by 1  */
			if (log(kiss(seed)) < alpha) {
				adelta_l++;
			}

			/*Otherwise, we need to revert values back to their current state */
			else {

				/*Set b and hl back equal to current values*/
				subject->basehalf_l[0] = currentmd[0];
				subject->basehalf_l[1] = currentmd[1];

				/*Set mean_contrib back equal to current values for each pulse*/
				i = 0;
				subnode = subject->driver->succ;
				while (subnode != NULL) {
					for (j = 0; j < N; j++) {
						subnode->mean_contrib[j] = currentmc[i][j];
					}
					i++;
					subnode = subnode->succ;
				}

				/*Set decay rate back equal to current value*/
				subject->decay_l = current_decay;

			} /*End of if else statement*/

		}

		/*Go to next subject*/
		subject = subject->succ;
		parms->subindex++;

		/*Free memory */
		for (i = 0; i < num_node; i++)
			free(currentmc[i]);
		free(currentmc);
		free(pmd);

	} /*End of loop through subjects*/

}
