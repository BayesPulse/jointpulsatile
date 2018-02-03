/*******************************************************************/
/*************************draw_pop_mass_mean.c *******************/
/*******************************************************************/

/***fix here***/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"

/*********************************************************************/
/*START OF draw_pop_mass_mean SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_pop_mass_mean: this is a metropolis-hastings draw of the both population mean pulse masses
	ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
	Priors *priors; the current values of the prior parameters;
	Common_parms *parms; the current values of the common parameters;
	unsigned long *seed; seed values needed for the randon number generator;
	Hyper_priors *hyper; parameters of the prior distributions;
	RETURNS: None; all updates are made internally
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
 fesum: sum of subject specific mean amplitudes or widths;
 part of evaluation of gmean and gvar
 gmean: mean of the distribution used in Gibbs sampler
 gvar: variance of distribution used in Gibbs sampler
 *subject: current list of subjects and their characteristics

 SUBROUTINES USED
 rnorm: found in randgen.h; draws from the normal distribution
 kiss: found in randgen.h; draws from U(0,1) distribution
**************************************************************************/

 void draw_pop_mass_mean(patientlist, popparms, popparms_response, popprior, popprior_response, seed);
void draw_fe_prior_a_mean(Subject_type *sublist, Priors *priors, Common_parms *parms_f, unsigned long *seed, Hyper_priors *hyper)
{
	/* declare variables */

	Subject_type *subject;

	int i, j, nsubj,flag;
	double sum_lmass[2], **fcond_var, **fcond_var_inv, temp_cond_mean[2], *cond_mean, temp[2];
	int rmvnorm(double *, double **, int, double *, unsigned long *, int);


	/*declare functions */
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);
	int cholesky_decomp(double **, int);
	double **cholesky_invert(int, double **);
	/* set up space for dynamically allocated variables */
	fcond_var = (double **)calloc(2, sizeof(double *));
	fcond_var_inv = (double **)calloc(2, sizeof(double *));
	cond_mean = (double *)calloc(2, sizeof(double));
	for (i = 0; i < 2; i++) {
		fcond_var[i] = (double *)calloc(2, sizeof(double));
		fcond_var_inv[i] = (double *)calloc(2, sizeof(double));
	}
	nsubj = parms_f->numsub;
	/* sum the current log masses */
	sum_lmass[0] = 0;
	sum_lmass[1] = 0;


	/*priors->fe_precision[0][0] = 5.263;
	priors->fe_precision[1][1] = 5.263;
	priors->fe_precision[0][1] = priors->fe_precision[1][0] = -4.737;*/



	subject = sublist->succ;

	/*This while loop goes through the pulses to get the sum involved*/
	while (subject != NULL){
		sum_lmass[0] += subject->theta_l[0]; /*sum of lh mean*/
		sum_lmass[1] += subject->theta_f[0];  /*sum of fsh mean*/
		subject = subject->succ;
	} /*end of loop through pulses*/
	for (i = 0; i < 2; i++)
	for (j = 0; j < 2; j++)
		fcond_var_inv[i][j] = (double)nsubj * priors->fe_precision[i][j] + hyper->prec[i][j]; /*the inverse*/

	if (!cholesky_decomp(fcond_var_inv, 2)) {
		printf("fcond_var_inv is not PSD matrix\n");
		exit(0);
	}
	/* invert to get the variance matrix.  This is needed for the mean calculation */
	/* and for the draw of the MVN distribution                                                 */
	fcond_var = cholesky_invert(2, fcond_var_inv);
	/* calculate the means of the full cond. Bivariate normal distn */
	for (i = 0; i < 2; i++) {
		temp_cond_mean[i] = 0;

		temp_cond_mean[i] += hyper->prec[i][0] * hyper->hmean_l[0] + hyper->prec[i][1] * hyper->hmean_f[0];
		temp_cond_mean[i] += priors->fe_precision[i][0] * sum_lmass[0] + priors->fe_precision[i][1] * sum_lmass[1];

	}

	for (i = 0; i < 2; i++) {
		cond_mean[i] = 0;
		for (j = 0; j < 2; j++)
			cond_mean[i] += fcond_var[i][j] * temp_cond_mean[j];
	}

	/*draw the new means*/
	flag =0;
	while (flag == 0){
		rmvnorm(temp, fcond_var, 2, cond_mean, seed, 1);
		if (temp[0]>0 && temp[1]>0)
			flag = 1;
	}
	priors->fe_mean_l[0] = temp[0];
	priors->fe_mean_f[0] = temp[1];
	for (i = 0; i < 2; i++) {
		free(fcond_var[i]);
		free(fcond_var_inv[i]);
	}
	free(fcond_var);
	free(fcond_var_inv);
	free(cond_mean);
}
