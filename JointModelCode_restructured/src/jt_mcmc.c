/*******************************************************************/
/*************************jt_mcmc.c *****************************/
/*******************************************************************/

/**FIX HERE**/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"


/*******************************************************************
*******************GLOBAL VARIABLE DEFINITIONS**********************

mmm: Order statistic used for distribution of pulse locations.
This is inputted by the user and is typically 3
adelta_f: counter that increments when a halflife estimate is accetped in the MH draw
ndelta_f: counter that increments for each loop of the halflife MH draw
adelta_l: counter that increments when a halflife estimate is accetped in the MH draw
ndelta_l: counter that increments for each loop of the halflife MH draw

atime_l: counter that increments when a pulse location estimate is accetped in the MH draw
ntime_l: counter that increments for each loop of the pulse location MH draw

atime_f: counter that increments when a pulse location estimate is accetped in the MH draw
ntime_f: counter that increments for each loop of the pulse location MH draw

arem_f: counter that increments when an individual pulse mass estimate is accetped in the MH draw
nrem_f: counter that increments for each loop of the individual pulse mass MH draw
arem_l: counter that increments when an individual pulse mass estimate is accetped in the MH draw
nrem_l: counter that increments for each loop of the individual pulse mass MH draw


arew_l: counter that increments when an individual pulse width estimate is accetped in the MH draw
nrew_l: counter that increments for each loop of the individual pulse width MH draw
arew_f: counter that increments when an individual pulse width estimate is accetped in the MH draw
nrew_f: counter that increments for each loop of the individual pulse width MH draw

afepmv_l: counter that increments when an overall mean pulse mass standard deviation (sma) is accetped in the MH draw
nfepmv_l: counter that increments for each loop of the overall mean pulse mass standard deviation (sma) MH draw
afepwv_l: counter that increments when an overall mean pulse width standard deviation (smw) is accetped in the MH draw
nfepwv_l: counter that increments for each loop of the overall mean pulse width standard deviation (smw) MH draw

afepmv_f: counter that increments when an overall mean pulse mass standard deviation (sma) is accetped in the MH draw
nfepmv_f: counter that increments for each loop of the overall mean pulse mass standard deviation (sma) MH draw
afepwv_f: counter that increments when an overall mean pulse width standard deviation (smw) is accetped in the MH draw
nfepwv_f: counter that increments for each loop of the overall mean pulse width standard deviation (smw) MH draw

afepm: counter that increments when an overall mean pulse mass standard deviation (sma) is accetped in the MH draw
nfepm: counter that increments for each loop of the overall mean pulse mass standard deviation (sma) MH draw
afepwv_l: counter that increments when an overall mean pulse width standard deviation (smw) is accetped in the MH draw
nfepwv_l: counter that increments for each loop of the overall mean pulse width standard deviation (smw) MH draw

afepwv_f: counter that increments when an overall mean pulse width standard deviation (smw) is accetped in the MH draw
nfepwv_f: counter that increments for each loop of the overall mean pulse width standard deviation (smw) MH draw

afemv: counter that increments when an overall pulse mass standard deviation (sa) is accetped in the MH draw
nfemv: counter that increments for each loop of the overall pulse mass standard deviation (sa) MH draw


afewv_l: counter that increments when an overall pulse width standard deviation (sw) is accetped in the MH draw
nfewv_l: counter that increments for each loop of the overall pulse width standard deviation (sw) MH draw
afewv_f: counter that increments when an overall pulse width standard deviation (sw) is accetped in the MH draw
nfewv_f: counter that increments for each loop of the overall pulse width standard deviation (sw) MH draw

afebv_l: counter that increments when a baseline standard deviation (sb) is accetped in the MH draw
nfebv_l: counter that increments for each loop of the baseline standard deviation (sb) MH draw
afehv_l: counter that increments when a half-life standard deviation (sh) is accetped in the MH draw
nfehv_l: counter that increments for each loop of the half-life standard deviation (sh) MH draw
afebv_f: counter that increments when a baseline standard deviation (sb) is accetped in the MH draw
nfebv_f: counter that increments for each loop of the baseline standard deviation (sb) MH draw
afehv_f: counter that increments when a half-life standard deviation (sh) is accetped in the MH draw
nfehv_f: counter that increments for each loop of the half-life standard deviation (sh) MH draw
ae_f:    for drawing the eta
ne_f:
ae_l:
ne_l:
fitstart: The first time in hours that a pulse may occur
fitend: The last time in hours that a pulse may occur

*********************************************************************/

extern int mmm;
long adelta_l = 0;
long ndelta_l = 0;
long atime_l = 0;
long ntime_l = 0;
long arem_l = 0;
long nrem_l = 0;
long arew_l = 0;
long nrew_l = 0;
long afepw_l = 0;
long nfepw_l = 0;
long ae_ml = 0;
long ne_ml = 0;
long ae_wl = 0;
long ne_wl = 0;

long afepmv_l = 0;
long nfepmv_l = 0;
long afepwv_l = 0;
long nfepwv_l = 0;
long afewv_l = 0;
long nfewv_l = 0;
long afebv_l = 0;
long nfebv_l = 0;
long afehv_l = 0;
long nfehv_l = 0;
long afepm = 0;
long nfepm = 0;
long arho = 0;
long nrho = 0;
long anu = 0;
long nnu = 0;


long adelta_f = 0;
long ndelta_f = 0;
long atime_f = 0;
long ntime_f = 0;
long arem_f = 0;
long nrem_f = 0;
long arew_f = 0;
long nrew_f = 0;
long afepw_f = 0;
long nfepw_f = 0;
long afepmv_f = 0;
long nfepmv_f = 0;
long afepwv_f = 0;
long nfepwv_f = 0;
long afewv_f = 0;
long nfewv_f = 0;
long afebv_f = 0;
long nfebv_f = 0;
long afehv_f = 0;
long nfehv_f = 0;
extern double fitstart;
extern double fitend;
long afemv = 0;
long nfemv = 0;
long ae_mf = 0;
long ne_mf = 0;
long ae_wf = 0;
long ne_wf = 0;

/********************************************************************/
/*SUBROUTINES THAT EXIST IN THIS PROGRAM

 mcmc
 

 **********************************************************************/


/*********************************************************************/
/*START OF MCMC SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*mcmc: this runs the BDMCMC process
	ARGUMENTS: Subject_type *sublist; this is the current list of subjects that exist;
	Common_parms *parms; the current values of the common parameters;
	double **ts; this is the matrix of observed data (a column of
	times and S columns of log(concentration);
	long iter; the number of iterations to run;
	int N; the number of observations in **ts;
	Priors *priors; the parameters of the prior distributions;
	unsigned long *seed; seed values needed for the randon number generator;
	char *file1; the output file name for common parameters
	double propvar[]; vector of proposal variances input by user
	Hyper_priors *hyper; hyperpriors input into algorithm by user
	RETURNS: None; all updates are made internally
*********************************************************************/
/*********************************************************************/
/***NOT UPDATED***/
/*VARIABLE DEFINITIONS
 i,j,k,l: generic counters
 num_node: counter of num_node for loops
 num_node2: number of pulses
 NN=50: Output every NNth iteration to files
 NNN=5000: Output every NNNth interation to screen or output file
 vtime: proposal variance for individual pulse locations
 vrem: proposal variance for individual masses
 vrew: proposal variance for individual widths
 vfepmv: proposal variance for standard deviation of overall mean mass, sma
 vfepwv: proposal variance for standard deviation of overall mean width, smw
 vfemv: proposal variance for standard deviation of pulse mass, sa
 vfewv: proposal variance for standard deviation of pulse width, sw
 vfebv: proposal variance for standard deviation of baseline, sb
 vfehv: proposal variance for standard deviation of half-life, sh
 ssq: sum of squared differences between log(concentration) and expected value
 **pmd_var: proposal variance-covariance matrix for baseline and halflife
 **pmd_vch: cholesky decomposed proposal variance-covariance matrix for
 baseline and halflife
 *new_node: index used to go through nodes in file output
 *subject: index used to go through subjects in file output
 *common: name given to inputted file name for common parameters
 *parm: name given to inputted file name for pulse parameters

 SUBROUTINES USED
 inverse_gamma: found in randgen.h; draws from the inverse gamma distribution
 error_squared: found in this file; calculates sum of squared error
 draw_bh: found in this file; draws baseline and half-life
 draw_times: found in this file; draws individual pulse locations
 birth_death: found in birthdeath.c; runs birthdeath process
 draw_fixed_effects: found in this file; draws subject specific mean mass and width
 draw_fe_precision: found in this file; draws st dev of pulse mass and width
 draw_random_effects: found in this file; draws individual masses and widths
 draw_bh_mean: found in this file; draws mean of baseline and half-life
 draw_bh_var: found in this file; draws st dev of baseline and half-life
 draw_fe_priors: found in this file; draws mean of mean pulse mass and width
 draw_fe_priors2: found in this file; draws st dev of mean pulse mass and width
 adjust_acceptance: found in this file; adjusts proposal variance
 adjust_acceptance2: found in this file; adjusts 2-d proposal variance*/
/***********************************************************************/
void mcmc(Subject_type *sublist, Common_parms *parms_l, Common_parms *parms_f, double **ts_l, double **ts_f, long iter, int N,
	Priors *priors, unsigned long *seed, char *filel, char *file2, double propvar[], Hyper_priors *hyper)
{
	int i, j, k, l, num_node, num_node2, NN = 50, NNN = 5000;
    double vfepmv_l, **sigma_p, vfepwv_l, vrem_l, vrew_l, vfemv_l, vfepw_l, vnu, vrho, vfepw_f, vfewv_l;
    double vfebv_l, vfehv_l, veta_l[2], veta_f[2], vfepmv_f, vfepwv_f, vrem_f, vrew_f, vfemv_f, vfewv_f;
    double vfebv_f, vfehv_f, vtime_l, ssq_l, vtime_f, ssq_f, **pmd_var_l, **pmd_vch_l, **pmd_var_f, **pmd_vch_f, **pmean_var;
	Node_type *new_node, *tdriver, *tresponse;
	Subject_type *subject;
	FILE *common_l, *common_f;
	int cholesky_decomp(double **A, int num_col);
	double *skernel(Node_type *, Node_type *, Common_parms *, int);
	double phi(double, double, double);
	double **cholesky_invert(int, double **);

	Common_parms *tparms_l, *tparms_f;
	Priors *tpriors;
	Subject_type *tsubject, *tsublist, *ttsublist, *ttsubject;

	double inverse_gamma(double alpha, double beta, unsigned long *seed);
	double error_squared_l(double **, Subject_type *, Common_parms *, int);
	double error_squared_f(double **, Subject_type *, Common_parms *, int);
	void draw_fe_prior_a_mean(Subject_type *sublist, Priors *priors, Common_parms *parms_f, unsigned long *seed, Hyper_priors *hyper);
	void draw_fe_prior_w_mean(Subject_type *sublist, Priors *priors, Common_parms *parms_f, unsigned long *seed, Hyper_priors *hyper);
	void draw_fe_prior_a_var(Subject_type *sublist, Priors *priors, Common_parms *parms_f, unsigned long *seed, Hyper_priors *hyper);
	void draw_bh_mean_f(Subject_type *sublist, Priors *priors, Common_parms *parms, unsigned long *seed, Hyper_priors *hyper);
	void draw_bh_mean_l(Subject_type *sublist, Priors *priors, Common_parms *parms, unsigned long *seed, Hyper_priors *hyper);
	void draw_fe_priors_w_var(Subject_type *sublist, Priors *priors, double v1, double v2, unsigned long *seed, Hyper_priors *hyper);

	void draw_bh_l(Subject_type *, Common_parms *, Priors *, double **,
		int, unsigned long *, double **);
	void draw_bh_f(Subject_type *, Common_parms *, Priors *, double **,
		int, unsigned long *, double **);


	void draw_times_f(Subject_type *, Common_parms *, Common_parms *, double **,
		int, unsigned  long *, double);
	void draw_times_l(Subject_type *, Common_parms *, double **,
		int, unsigned long *, double);
	void draw_times_old(Subject_type *, Common_parms *, double **,
		int, unsigned long *, double);
	void draw_times_old_f(Subject_type *, Common_parms *, double **,
		int, unsigned long *, double);

	void draw_bh_var_l(Subject_type *, Priors *, double, double, unsigned long *, Hyper_priors *);
	void draw_bh_var_f(Subject_type *, Priors *, double, double, unsigned long *, Hyper_priors *);


	void birth_death_l(Subject_type *, double **, Common_parms *, int,
		unsigned long *, int);
	void birth_death_f(Subject_type *, double **, Common_parms *, Common_parms *, int,
		unsigned long *, int);

	void draw_fixed_mass(Subject_type *, Priors *, Common_parms *, Common_parms *, unsigned long *, double **);
	void draw_fixed_width_l(Subject_type *, Priors *, Common_parms *, double, unsigned long *);
	void draw_fixed_width_f(Subject_type *, Priors *, Common_parms *, double, unsigned long *);
	void draw_fe_precision_l(Subject_type *, Priors *, Common_parms *, double, double, unsigned long *);
	void draw_fe_precision_f(Subject_type *, Priors *, Common_parms *, double, double, unsigned long *);

	void draw_bh_l(Subject_type *, Common_parms *, Priors *, double **,
		int, unsigned long *, double **);
	void draw_bh_f(Subject_type *, Common_parms *, Priors *, double **,
		int, unsigned long *, double **);
	void draw_random_effects_l(double **, Subject_type *, Common_parms *, int, double, double, unsigned long *);
	void draw_random_effects_f(double **, Subject_type *, Common_parms *, int, double, double, unsigned long *);


	void draw_eta_l(Subject_type *, Common_parms *, unsigned long *, double *);
	void draw_eta_f(Subject_type *, Common_parms *, unsigned long *, double *);
	void destroy_sublist(Subject_type *);
	Subject_type *initialize_subject(void);
	void insert_subject(Subject_type *, Subject_type *);

	double Phi(double y);
	void adjust_acceptance(double x, double *X);
	void adjust2_acceptance(double x, double **X, double corr);
	double gamm(double x);
	double gammalog(double x, double a, double b);
	double gammaPdf(double x, double a, double b);

	void adjust_acceptance(double, double *);
	void adjust2_acceptance(double, double **, double);
	void mh_logalphamean(Subject_type *, Common_parms *, Priors *, unsigned long *, double);
	void mh_logsigmean(Subject_type *, Common_parms *, Priors *, unsigned long *, double);



	
 
strcat(filel, ".out");
	common_l = fopen(filel, "w");

	strcat(file2, ".out");
	common_f = fopen(file2, "w");

	subject = sublist->succ;
	while (subject != NULL){
		subject->csub_l = fopen(subject->common_l, "w");
		subject->psub_l = fopen(subject->pulse_l, "w");
		subject->csub_f = fopen(subject->common_f, "w");
		subject->psub_f = fopen(subject->pulse_f, "w");
		subject = subject->succ;
	}
	sigma_p = (double **)calloc(2, sizeof(double *));
	for (i = 0; i < 2; i++)
		sigma_p[i] = (double *)calloc(2, sizeof(double));
	/*proposal variances for poulation variance of pulse width,sigma_wl,sigma_wf*/
	vfewv_l = propvar[0];
	vfewv_f = propvar[1];

	/*proposal variances for population variance of baseline and halflife,sigma_b, sigma_h)*/
	vfebv_l = propvar[2];
	vfehv_l = propvar[3];
	vfebv_f = propvar[4];
	vfehv_f = propvar[5];
	/*proosal variancee for subject pulse mass pair*/
	/*read in three parameters, lh, fsh, and correlation*/
	/*the output is and 2 by 2 matrix,big_sigma_m*/
	pmean_var = (double **)calloc(2, sizeof(double *));
	for (i = 0; i < 2; i++)
		pmean_var[i] = (double *)calloc(2, sizeof(double));
	pmean_var[0][0] = propvar[6];
	pmean_var[1][1] = propvar[7];
	pmean_var[0][1] = pmean_var[1][0] = sqrt(pmean_var[1][1] * pmean_var[0][0])* propvar[8];
	/*proposal variances for individual width mean,mu_wi*/
	vfepw_l = propvar[9];
	vfepw_f = propvar[10];

	/*proposal variances for common variances of subject level v_pm */
	vfepmv_l = propvar[11];
	vfepwv_l = propvar[12];
	vfepmv_f = propvar[13];
	vfepwv_f = propvar[14];

	/*proposal variances for subject baselines and halflifes*/
	pmd_var_l = (double **)calloc(2, sizeof(double *));
	for (i = 0; i < 2; i++)
		pmd_var_l[i] = (double *)calloc(2, sizeof(double));
	pmd_var_f = (double **)calloc(2, sizeof(double *));
	for (i = 0; i < 2; i++)
		pmd_var_f[i] = (double *)calloc(2, sizeof(double));


	pmd_var_l[0][0] = propvar[15];
	pmd_var_l[1][1] = propvar[16];
	pmd_var_l[0][1] = pmd_var_l[1][0] = propvar[17] * sqrt(pmd_var_l[0][0])*sqrt(pmd_var_l[1][1]);
	pmd_var_f[0][0] = propvar[18];
	pmd_var_f[1][1] = propvar[19];
	pmd_var_f[0][1] = pmd_var_f[1][0] = propvar[20] * sqrt(pmd_var_f[0][0])*sqrt(pmd_var_f[1][1]);

	/*proposa variance for pulse mass and width*/
	vrem_l = propvar[21];
	vrew_l = propvar[22];
	vrem_f = propvar[23];
	vrew_f = propvar[24];


	/*proposal variance for pulse locations*/
	vtime_l = propvar[25];
	vtime_f = propvar[26];

	/*eta*/
	veta_l[0] = propvar[27];
	veta_l[1] = propvar[28];

	veta_f[0] = propvar[29];
	veta_f[1] = propvar[30];
	vrho = propvar[31];
	vnu = propvar[32];

	/***cholesky decompose the proposal var-covar matrix for b and hl***/
	pmd_vch_l = (double **)calloc(2, sizeof(double *));
	for (i = 0; i < 2; i++)
		pmd_vch_l[i] = (double *)calloc(2, sizeof(double));

	for (i = 0; i < 2; i++)
	for (j = 0; j < 2; j++)
		pmd_vch_l[i][j] = pmd_var_l[i][j];

	if (!cholesky_decomp(pmd_vch_l, 2)){
		printf("not PSD matrix A\n");
		exit(0);
	}


	pmd_vch_f = (double **)calloc(2, sizeof(double *));
	for (i = 0; i < 2; i++)
		pmd_vch_f[i] = (double *)calloc(2, sizeof(double));

	for (i = 0; i < 2; i++)
	for (j = 0; j < 2; j++)
		pmd_vch_f[i][j] = pmd_var_f[i][j];

	if (!cholesky_decomp(pmd_vch_f, 2)){
		printf("not PSD matrix A\n");
		exit(0);
	}

	/*********type in the true value  ***********/

	for (i = 0; i < iter; i++) {

		/* run the birth-death algorithm to select pulses */

		parms_f->iter = i;
		parms_l->iter = i;



		/* draw fixed effects priors; mua and muw,the population mean*/

		draw_fe_prior_a_mean(sublist, priors, parms_f, seed, hyper);
		draw_fe_prior_w_mean(sublist, priors, parms_f, seed, hyper);



		/* draw fixed effects prior variances; sma and smw*/

		draw_fe_prior_a_var(sublist, priors, parms_f, seed, hyper);
		draw_fe_priors_w_var(sublist, priors, vfewv_l, vfewv_f, seed, hyper);  /*proposal variance */

		/* draw mean baseline and halflife; mub and muh*/

		draw_bh_mean_f(sublist, priors, parms_f, seed, hyper);
		draw_bh_mean_l(sublist, priors, parms_l, seed, hyper);

		///* draw variance for baseline and halflife; sb and sh*/

		draw_bh_var_l(sublist, priors, vfebv_l, vfehv_l, seed, hyper);
		draw_bh_var_f(sublist, priors, vfebv_f, vfehv_f, seed, hyper);

		/////* draw subject specific means; muak and muwk*/

		draw_fixed_mass(sublist, priors, parms_l, parms_f, seed, pmean_var);
		draw_fixed_width_l(sublist, priors, parms_l, vfepw_l, seed);
		draw_fixed_width_f(sublist, priors, parms_f, vfepw_f, seed);

		/* draw mass and width variances; sa and sw*/

		draw_fe_precision_l(sublist, priors, parms_l, vfepmv_l, vfepwv_l, seed);
		draw_fe_precision_f(sublist, priors, parms_f, vfepmv_f, vfepwv_f, seed);

		birth_death_l(sublist, ts_l, parms_l, N, seed, i);
		//
		
			draw_times_l(sublist, parms_l, ts_l,
				N, seed, vtime_l);
		
			////////
			//draw_times_old(tsublist, tparms_l, ts_l,
			//	N, seed, vtime_l);
			draw_random_effects_l(ts_l, sublist, parms_l, N, vrem_l, vrew_l, seed);
			draw_eta_l(sublist, parms_l, seed, veta_l);
			

		draw_bh_l(sublist, parms_l, priors, ts_l,
			N, seed, pmd_var_l);
		//////

		/////* draw pulse locations; tauki*/
		birth_death_f(sublist, ts_f, parms_f, parms_l, N, seed, i);
		////draw_times_old_f(tsublist, tparms_f, ts_f,
		////	N, seed, vtime_f);
		
		draw_times_f(sublist, parms_f, parms_l, ts_f,
			N, seed, vtime_f);
		draw_eta_f(sublist, parms_f, seed, veta_f);
	
		//
		draw_random_effects_f(ts_f, sublist, parms_f, N, vrem_f, vrew_f, seed);

		draw_bh_f(sublist, parms_f, priors, ts_f,
			N, seed, pmd_var_f);


		/////* draw pulse mass and width; Aki and s2pki*/

		
		mh_logalphamean(sublist, parms_l, priors, seed, vrho);
		mh_logsigmean(sublist, parms_l, priors, seed, vnu);


		/* draw model error; s2e*/
		ssq_l = error_squared_l(ts_l, sublist, parms_l, N);
		parms_l->sigma = inverse_gamma(priors->alpha_l + (parms_l->numsub*N) / 2, priors->beta_l + 0.5*ssq_l, seed);
		parms_l->lsigma = log(parms_l->sigma);
		ssq_f = error_squared_l(ts_f, sublist, parms_f, N);
		parms_f->sigma = inverse_gamma(priors->alpha_f + (parms_f->numsub*N) / 2, priors->beta_f + 0.5*ssq_f, seed);
		parms_f->lsigma = log(parms_f->sigma);

		fflush(stdout);
        
        /**FIX HERE!!**/
/**        if iteration == 50 or 5000 do write out loops
        write_mcmc_output(MAKE A NEW FUNCTION FOUND IN WRITE_MCMC_OUTPUT.C) **/
        if (!(i % 500) && i<25000 && i >0) {
            
            adjust2_acceptance((double)adelta_l / (double)ndelta_l, pmd_var_l, propvar[17]);
            adjust2_acceptance((double)adelta_f / (double)ndelta_f, pmd_var_f, propvar[20]);
            adjust_acceptance((double)atime_l / (double)ntime_l, &vtime_l);
            adjust_acceptance((double)arem_l / (double)nrem_l, &vrem_l);
            adjust_acceptance((double)arew_l / (double)nrew_l, &vrew_l);
            adjust_acceptance((double)atime_f / (double)ntime_f, &vtime_f);
            adjust_acceptance((double)arem_f / (double)nrem_f, &vrem_f);
            adjust_acceptance((double)arew_f / (double)nrew_f, &vrew_f);
            adjust_acceptance((double)afepwv_l / (double)nfepwv_l, &vfepwv_l);
            adjust_acceptance((double)afepmv_l / (double)nfepmv_l, &vfepmv_l);
            adjust_acceptance((double)afepwv_f / (double)nfepwv_f, &vfepwv_f);
            adjust_acceptance((double)afepmv_f / (double)nfepmv_f, &vfepmv_f);
            adjust_acceptance((double)afepwv_l / (double)nfepwv_l, &vfepwv_l);
            adjust_acceptance((double)afepwv_f / (double)nfepwv_f, &vfepwv_f);
            adjust_acceptance((double)ae_ml / (double)ne_ml, &veta_l[0]);
            
            
            
            adjust2_acceptance((double)afemv / (double)nfemv, pmean_var, propvar[8]);
            adjust_acceptance((double)afewv_l / (double)nfewv_l, &vfewv_l);
            adjust_acceptance((double)afebv_l / (double)nfebv_l, &vfebv_l);
            adjust_acceptance((double)afehv_l / (double)nfehv_l, &vfehv_l);
            adjust_acceptance((double)afewv_f / (double)nfewv_f, &vfewv_f);
            adjust_acceptance((double)afebv_f / (double)nfebv_f, &vfebv_f);
            adjust_acceptance((double)afehv_f / (double)nfehv_f, &vfehv_f);
            
            
            /*After adjusting proposal variances, reset the counters for acceptance rates*/
            adelta_l = ndelta_l = 0;
            adelta_f = ndelta_f = 0;
            
            atime_l = ntime_l = 0;
            arem_l = nrem_l = 0;
            arew_l = nrew_l = 0;
            atime_f = ntime_f = 0;
            arem_f = nrem_f = 0;
            arew_f = nrew_f = 0;
            afepmv_l = nfepmv_l = 0;
            afepwv_l = nfepwv_l = 0;
            afepwv_f = nfepwv_f = 0;
            afepmv_f = nfepmv_f = 0;
            
            afemv = nfemv = 0;
            afewv_l = nfewv_l = 0;
            afebv_l = nfebv_l = 0;
            afehv_l = nfehv_l = 0;
            afewv_f = nfewv_f = 0;
            afebv_f = nfebv_f = 0;
            afehv_f = nfehv_f = 0;
            ae_wf = ne_wf = 0;
            ae_wl = ne_wl = 0;
            ae_mf = ne_mf = 0;
            ae_ml = ne_ml = 0;
            
            /***cholesky decompose the proposal var-covar matrix***/
            for (k = 0; k < 2; k++)
                for (l = 0; l < 2; l++)
                    pmd_vch_l[k][l] = pmd_var_l[k][l];
            
            if (!cholesky_decomp(pmd_vch_l, 2)){
                printf("pmd not PSD matrix\n");
                exit(0);
            }
            for (k = 0; k < 2; k++)
                for (l = 0; l < 2; l++)
                    pmd_vch_f[k][l] = pmd_var_f[k][l];
            
            if (!cholesky_decomp(pmd_vch_f, 2)){
                printf("pmd not PSD matrix\n");
                exit(0);
            }
            
        }
    }
    
} /*End of loop through iterations*/

for (i = 0; i < 2; i++)
free(sigma_p[i]);
free(sigma_p);
for (i = 0; i < 2; i++)
free(pmd_var_f[i]);
free(pmd_var_f);
for (i = 0; i < 2; i++)
free(pmd_var_l[i]);
free(pmd_var_l);
for (i = 0; i < 2; i++)
free(pmd_vch_f[i]);
free(pmd_vch_f);
for (i = 0; i < 2; i++)
free(pmd_vch_l[i]);
free(pmd_vch_l);

fclose(common_l);



fclose(common_f);


} /*End of MCMC */

/***************************************************************/

