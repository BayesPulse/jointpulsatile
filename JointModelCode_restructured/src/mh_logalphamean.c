/*******************************************************************/
/*************************linklistv2.c *****************************/
/*******************************************************************/

#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"



/***MH algorithm for drawing the mean cluster size for the response pulse location model***/

void mh_logalphamean(Subject_type *subject, Common_parms *parms_cox, Priors *priors_cox, unsigned long *seed, double palphamean_sd) {
	int nparents, npulse;
	double  NewMeanAlpha, NewMeanLogAlpha, lprop_ratio, lprior_ratio, lratio, alpha, integral, pulselocs;
	Node_type *tmpparent, *fnode;
	Node_type *tmplist;
	Subject_type  * subject2;

	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);
	double phi(double, double, double);

	NewMeanLogAlpha = rnorm(parms_cox->lrho, palphamean_sd, seed);
	NewMeanAlpha = exp(NewMeanLogAlpha);

	/*   printf("newmeanalpha %lf palphamean_sd %lf\n",NewMeanAlpha,palphamean_sd);*/
	nrho++;
	lratio = 0;
	lprior_ratio = 0;
	lprop_ratio = 0;

	subject2 = subject->succ;  /****pull information from all the subject  ***/;
	while (subject2 != NULL)

	{
		nparents = subject2->numnode_l;
		npulse = subject2->numnode_f;
		tmpparent = subject2->driver->succ;
		integral = 0;
		while (tmpparent != NULL){
			pulselocs = tmpparent->time;
			integral = phi(fitend, pulselocs, sqrt(parms_cox->nu)) - phi(fitstart, pulselocs, sqrt(parms_cox->nu));
			integral = (-NewMeanAlpha + exp(parms_cox->lrho))*integral;
			lratio += integral;
			tmpparent = tmpparent->succ;
		}
		lratio += npulse*(NewMeanLogAlpha - parms_cox->lrho);

		subject2 = subject2->succ;
	}


	/*    printf("lratio pt1 %lf\n",lratio);*/

	/*log ratio of the prior distributions */
	lprior_ratio = 1 / (2 * priors_cox->rho_prior[1]) * ((parms_cox->lrho - priors_cox->rho_prior[0]) * (parms_cox->lrho - priors_cox->rho_prior[0]) - (NewMeanLogAlpha - priors_cox->rho_prior[0]) * (NewMeanLogAlpha - priors_cox->rho_prior[0]));

	lprop_ratio = lratio + lprior_ratio;

	/*  printf("lprop_ratio %lf lratio %lf lprior_ra  tio %lf\n",lprop_ratio,lratio,lprior_ratio);*/

	alpha = (0 < lprop_ratio) ? 0 : lprop_ratio;

	if (log(kiss(seed)) < alpha) {
		/*  printf("accept\n"); */
		arho++;
		parms_cox->lrho = NewMeanLogAlpha;
		parms_cox->rho = exp(parms_cox->lrho);
		/*update the conditional density function for each response pulse*/
		subject2 = subject->succ;
		while (subject2 != NULL)
		{
			fnode = subject2->response->succ;
			while (fnode != NULL)
			{
				fnode->lambda = skernel_1(fnode, subject2->driver, parms_cox);

				fnode = fnode->succ;
			}
			subject2 = subject2->succ;
		}
	}
}


