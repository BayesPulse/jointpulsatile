/*******************************************************************/
/*************************draw_fixed_width_trigger.c ***************/
/*******************************************************************/

/**fix here**/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"


void draw_fixed_width_trigger(Subject_type *sublist, Priors *priors, Common_parms *parms, double v, unsigned long *seed)
{
	int j, tmp1;
	double new_mean, prop_new, prop_old, psum_old, psum_new, pcomp, old_prior, new_prior;
	double prior_ratio, prop_ratio, temp, alpha;
	double stdxold, stdxnew, newint, oldint, prior_old, prior_new;
	Node_type *subnode;
	Subject_type *subject;
	double kiss(unsigned long *);
	double Phi(double);
	double rnorm(double, double, unsigned long *);

	/*Allocate Memory*/

	/*Add 1 to the counters for acceptance rates of sigma_a and sigma_w*/
	nfepw_l++;


	/*Assign current acceptance counts to temporary vector*/
	tmp1 = afepw_l;




	/* do the simulation for each subject */

	subject = sublist->succ;

	while (subject != NULL){
		new_mean = rnorm(subject->theta_l[1], v, seed);
		if (new_mean > 0){
			old_prior = (subject->theta_l[1] - priors->fe_mean_l[1])*(subject->theta_l[1] - priors->fe_mean_l[1]);
			old_prior /= priors->fe_precision_wl;
			old_prior /= priors->fe_precision_wl;
			new_prior = (new_mean - priors->fe_mean_l[1])*(new_mean - priors->fe_mean_l[1]);
			new_prior /= priors->fe_precision_wl;
			new_prior /= priors->fe_precision_wl;

			prior_ratio = 0.5*(old_prior - new_prior);
			psum_old = 0.0;
			psum_new = 0.0;
			newint = 0.0;
			oldint = 0.0;
			subnode = subject->driver->succ;
			while (subnode != NULL){
				psum_old += (subnode->theta[1] - subject->theta_l[1])*(subnode->theta[1] - subject->theta_l[1])*subnode->eta[1];
				psum_new += (subnode->theta[1] - new_mean)*(subnode->theta[1] - new_mean)*subnode->eta[1];
				stdxnew = (new_mean) / (parms->re_precision[1] / sqrt(subnode->eta[1]));
				stdxold = (subject->theta_l[1]) / (parms->re_precision[1] / sqrt(subnode->eta[1]));
				newint += log(Phi(stdxnew));
				oldint += log(Phi(stdxold));

				subnode = subnode->succ;
			}

			prop_ratio = (psum_old - psum_new)*0.5 / (parms->re_precision[1] * parms->re_precision[1]) - newint + oldint;

			alpha = (0 < (temp = (prior_ratio + prop_ratio))) ? 0 : temp;
			/*If log(U) < log rho, accept the proposed value*/
			/*Increase acceptance count by 1*/
			if (log(kiss(seed)) < alpha){
				tmp1++;
				subject->theta_l[1] = new_mean;



			}
		}

		subject = subject->succ;
	}





	/*end of loop through s_a and s_w */

	/*Set acceptance count equal to temp vector components*/

	afewv_l = tmp1;



}
