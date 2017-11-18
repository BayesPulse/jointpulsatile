/*******************************************************************/
/*************************draw_bh_var_response.c *******************/
/*******************************************************************/

/***fix here***/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"


void draw_bh_var_response(Subject_type *sublist, Priors *priors, double v1, double v2, unsigned long *seed, Hyper_priors *hyper)
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
	nfebv_f++;
	nfehv_f++;

	/*Assign current acceptance counts to temporary vector*/
	tmp1[0] = afebv_f;
	tmp1[1] = afehv_f;

	/*Draw proposed values for sigma_b and sigma_h*/
	new_var[0] = rnorm(priors->varbh_f[0], v1, seed);
	new_var[1] = rnorm(priors->varbh_f[1], v2, seed);

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
			prior_ratio *= priors->varbh_f[j] / new_var[j];
			pcomp = subject->basehalf_f[j] - priors->meanbh_f[j];
			pcomp *= pcomp;
			psum += pcomp;
			subject = subject->succ;
		}

		/*Log prior ratio value*/
		prior_ratio = log(prior_ratio);

		/*Complete the "likelihood" portion of rho*/
		prop_old = 1.0 / priors->varbh_f[j];
		prop_old /= priors->varbh_f[j];
		prop_new = 1.0 / new_var[j];
		prop_new /= new_var[j];
		prop_ratio = psum*0.5*(prop_old - prop_new);

		/*We only can accept the proposed value if it is positive*/
		/*With a uniform(0,2) prior, skip this step if we have a sd less than 2*/
		if (new_var[j] > 0 && new_var[j] < hyper->varbh_f[j]){

			/*Compute log rho, and set alpha equal to min(log rho,0)*/
			/*With a uniform prior, get rid of prior_ratio*/
			alpha = (0 < (temp = (prior_ratio + prop_ratio))) ? 0 : temp;

			/*If log(U) < log rho, accept the proposed value*/
			/*Increase acceptance count by 1*/
			if (log(kiss(seed)) < alpha){
				tmp1[j]++;
				priors->varbh_f[j] = new_var[j];
			}
		}
	}

	/*Set acceptance count equal to temp vector components*/
	afebv_f = tmp1[0];
	afehv_f = tmp1[1];

	free(new_var);
	free(tmp1);

}

