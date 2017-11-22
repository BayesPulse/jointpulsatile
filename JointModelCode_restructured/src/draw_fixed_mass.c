/*******************************************************************/
/*************************draw_fixed_mass.c *****************************/
/*******************************************************************/

/**fix here**/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"



void draw_fixed_mass(Subject_type *sublist, Priors *priors, Common_parms *parms_l, Common_parms *parms_f, unsigned long *seed, double **pmean_var)
{
	int i, j, k, numnode_l, numnode_f, tmp1[2];
	double diff_mean_old[2], diff_mean_new[2], log_prior_ratio, temp_old_mean[2], theta[2], temp, alpha, *new_mean, *tmpold, *tmpnew;
	double A_inv[2][2], **re_var_inv, **re_var, stdnew[2], stdold[2], log_RE_ratio, new_sum[2], old_sum[2], **pvartmp;
	double newint, oldint, junk, iden[2][2], stdxold, stdyold, stdxnew, stdynew, corrxy, log_prop_ratio;


	Subject_type *subject;
	Node_type *subnode_f;
	Node_type *subnode_l;
	double MDBNOR(double, double, double);
	int rmvnorm(double *, double **, int, double *, unsigned long *, int);
	double rgamma(double, double, unsigned long *);
	double kiss(unsigned long *);
	void print_list(Node_type *);
	int cholesky_decomp(double **, int);
	double **cholesky_invert(int, double **);
	double Phi(double y);

	double kiss(unsigned long *);
	double rnorm(double, double, unsigned long *);
	tmpold = (double *)calloc(2, sizeof(double));
	tmpnew = (double *)calloc(2, sizeof(double));
	new_mean = (double *)calloc(2, sizeof(double));
	re_var_inv = (double **)calloc(2, sizeof(double *));
	re_var = (double **)calloc(2, sizeof(double *));
	pvartmp = (double **)calloc(2, sizeof(double *));
	/*for (i = 0; i < 2; i++) {
		re_var_inv[i] = (double *)calloc(2, sizeof(double));
		re_var[i] = (double *)calloc(2, sizeof(double));
		pvartmp[i] = (double *)calloc(2, sizeof(double));
	}*/
	/*    printf("iter %d \n",iter);*/




	/*re_var = cholesky_invert(2,priors->fe_precision);*/

	/*Go to start of subject list*/
	subject = sublist->succ;

	/*Go through each existing subject*/
	while (subject != NULL) {


		nfepm++;

		tmp1[0] = afepm;
		/*draw the new pulse pair*/
		theta[0] = subject->theta_l[0];
		theta[1] = subject->theta_f[0];
		old_sum[0] = old_sum[1] = new_sum[0] = new_sum[1] = 0.0;

		/*Go to start of pulse list for this subject*/
		rmvnorm(new_mean, pmean_var, 2, theta, seed, 1);
		if (new_mean[0] > 0 && new_mean[1] > 0) {

			/*calculate the prior ratio*/
			log_prior_ratio = 0;
			diff_mean_old[0] = subject->theta_l[0] - priors->fe_mean_l[0];
			diff_mean_old[1] = subject->theta_f[0] - priors->fe_mean_f[0];
			diff_mean_new[0] = new_mean[0] - priors->fe_mean_l[0];
			diff_mean_new[1] = new_mean[1] - priors->fe_mean_f[0];

			for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++)
				log_prior_ratio += diff_mean_old[i] * priors->fe_precision[i][j] * diff_mean_old[j] - diff_mean_new[i] * priors->fe_precision[i][j] * diff_mean_new[j];

			log_prior_ratio *= 0.5;

			/*calculate the likelihood ratio*/
			/*first the driver's pulse*/
			subnode_l = subject->driver->succ;
			numnode_l = 0;
			newint = 0.;
			oldint = 0.;

			/*Go through each existing pulse and add up log amplitude/width
			 Also count pulses*/
			while (subnode_l != NULL){

				numnode_l++;
				old_sum[0] += (subnode_l->theta[0] - subject->theta_l[0])*(subnode_l->theta[0] - subject->theta_l[0])*subnode_l->eta[0];
				new_sum[0] += (subnode_l->theta[0] - new_mean[0])*(subnode_l->theta[0] - new_mean[0])*subnode_l->eta[0];
				stdxnew = (new_mean[0]) / (parms_l->re_precision[0] / sqrt(subnode_l->eta[0]));
				stdxold = (subject->theta_l[0]) / (parms_l->re_precision[0] / sqrt(subnode_l->eta[0]));
				newint += log(Phi(stdxnew));
				oldint += log(Phi(stdxold));

				subnode_l = subnode_l->succ;

			}
			log_RE_ratio = 0.5/(parms_l->re_precision[0]* parms_l->re_precision[0])*(old_sum[0] - new_sum[0]) + oldint - newint;
			/*start of the fsh loop*/
			subnode_f = subject->response->succ;
			numnode_f = 0;
			newint = 0.;
			oldint = 0.;
			/*Go through each existing pulse and add up log amplitude/width
			 Also count pulses*/
			while (subnode_f != NULL){

				numnode_f++;
				old_sum[1] += (subnode_f->theta[0] - subject->theta_f[0])*(subnode_f->theta[0] - subject->theta_f[0])*subnode_f->eta[0];
				new_sum[1] += (subnode_f->theta[0] - new_mean[1])*(subnode_f->theta[0] - new_mean[1])*subnode_f->eta[0];
				stdxnew = (new_mean[1]) / (parms_f->re_precision[0] / sqrt(subnode_f->eta[0]));
				stdxold = (subject->theta_f[0]) / (parms_f->re_precision[0] /sqrt(subnode_f->eta[0]));
				newint += log(Phi(stdxnew));
				oldint += log(Phi(stdxold));
				subnode_f = subnode_f->succ;


				log_RE_ratio += 0.5/ (parms_f->re_precision[0]* parms_f->re_precision[0])*(old_sum[1] - new_sum[1]) + oldint - newint;

			}

			alpha = (0 < (temp = (log_prior_ratio + log_RE_ratio))) ? 0 : temp;
			/*If log(U) < log rho, accept the proposed value*/
			/*Increase acceptance count by 1*/
			if (log(kiss(seed)) < alpha){
				tmp1[0]++;
				subject->theta_l[0] = new_mean[0];
				subject->theta_f[0] = new_mean[1];

			}
		}
		afepm = tmp1[0];
		/*Go through each existing pulse and add up log amplitude/width
		 Also count pulses*/


		/*Advance to next subject*/
		subject = subject->succ;

	} /*end of loop through subjects*/
	free(re_var_inv);
	free(re_var);
	free(new_mean);
	free(tmpnew);
	free(tmpold);
}

