/*******************************************************************/
/*************************jt_mcmc.c *****************************/
/*******************************************************************/

/**FIX HERE**/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"

extern double M;

//Start of the function: Needs updating
//void mcmc(Subject_type *sublist, Common_parms *parms_l, Common_parms *parms_f, double **ts_l, double **ts_f, int MCMCiter, int N,
	Priors *priors, unsigned long *seed, char *filel, char *file2, double propvar[], Hyper_priors *hyper)
void mcmc(Patient *patientlist, PopulationPriors *popprior, PopulationPriors *popprior_response, PopulationEstimates *popparms, PopulationEstimates *popparms_response,
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

//Run the MCMC for Niter iterations
	for (i = 0; i < MCMCinter; i++) {

		/* run the birth-death algorithm to select pulses */

        //Why do I need iteration number assigned to variables
//		parms_f->iter = i;
//		parms_l->iter = i;


//Draw population mean parameters for pulses: mua and muw
		/* draw fixed effects priors; mua and muw,the population mean*/

		draw_pop_mass_mean(sublist, priors, parms_f, seed, hyper);
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

