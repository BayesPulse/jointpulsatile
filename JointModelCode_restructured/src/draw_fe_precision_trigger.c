/*******************************************************************/
/*************************linklistv2.c *****************************/
/*******************************************************************/

/**fixed here***/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"


void draw_fe_precision_trigger(Subject_type *sublist, Priors *priors, Common_parms *parms, double v1, double v2, unsigned long *seed)
{
	int *tmp1, j;
	double *new_var, prop_new, prop_old, psum, pcomp, prior_old, prior_new;
	double prior_ratio, prop_ratio, temp, alpha;
	double stdold[2], stdnew[2], newint[2], oldint[2];
	Node_type *subnode;
	Subject_type *subject;
	double kiss(unsigned long *);
	double Phi(double);
	double rnorm(double, double, unsigned long *);

	/*Allocate Memory*/
	new_var = (double *)calloc(2, sizeof(double));
	tmp1 = (int *)calloc(2, sizeof(int));

	/*Add 1 to the counters for acceptance rates of sigma_a and sigma_w*/
	nfepmv_l++;
	nfepwv_l++;

	/*Assign current acceptance counts to temporary vector*/
	tmp1[0] = afepmv_l;
	tmp1[1] = afepwv_l;

	/*Draw proposed values for sigma_a and sigma_w*/
	new_var[0] = rnorm(parms->re_precision[0], v1, seed);
	new_var[1] = rnorm(parms->re_precision[1], v2, seed);

	/*Accept or Reject sigma_a, then accept or reject for sigma_w*/
	for (j = 0; j < 2; j++){

		/*Calculate ratio of Cauchy priors*/
		/*      prior_new = new_var[j]/priors->re_var[j];
			  prior_old = parms->re_precision[j]/priors->re_var[j];
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
		prior_ratio = 1;
		newint[0] = 0;
		oldint[0] = 0;
		newint[1] = 0;
		oldint[1] = 0;

		subject = sublist->succ;
		while (subject != NULL){
			subnode = subject->driver->succ;
			while (subnode != NULL){
				prior_ratio *= parms->re_precision[j] / new_var[j];
				pcomp = subnode->theta[j] - subject->theta_l[j];
				pcomp *= pcomp;
				psum += pcomp*subnode->eta[j];
				stdold[j] = (subnode->theta[j]) / (parms->re_precision[j] / sqrt(subnode->eta[j]));
				stdnew[j] = (subnode->theta[j]) / (new_var[j] / sqrt(subnode->eta[j]));

				newint[j] += log(Phi(stdnew[j]));
				oldint[j] += log(Phi(stdold[j]));


				subnode = subnode->succ;
			}
			subject = subject->succ;
		}

		/*Log prior ratio value*/
		prior_ratio = log(prior_ratio) - newint[j] + oldint[j];

		/*Complete the "expo " portion of likelihood ratio*/
		prop_old = 1.0 / parms->re_precision[j];
		prop_old /= parms->re_precision[j];
		prop_new = 1.0 / new_var[j];
		prop_new /= new_var[j];
		prop_ratio = psum*0.5*(prop_old - prop_new);

		/*We only can accept the proposed value if it is positive*/
		if (new_var[j] > 0.1 && new_var[j] < priors->re_var_l[j]){

			/*Compute log rho, and set alpha equal to min(log rho,0)*/
			alpha = (0 < (temp = (prior_ratio + prop_ratio))) ? 0 : temp;
			/*If log(U) < log rho, accept the proposed value*/
			/*Increase acceptance count by 1*/
			if (log(kiss(seed)) < alpha){
				tmp1[j]++;
				parms->re_precision[j] = new_var[j];
			}
		}
	} /*end of loop through s_a and s_w */

	/*Set acceptance count equal to temp vector components*/
	afepmv_l = tmp1[0];
	afepwv_l = tmp1[1];

	free(new_var);
	free(tmp1);

}
