/*******************************************************************/
/*************************jt_mcmc.c *****************************/
/*******************************************************************/

#include "jt_mcmc.h"

extern double M;

//Start of the function: Needs updating
//void mcmc(Subject_type *sublist, Common_parms *parms_l, Common_parms *parms_f, double **ts_l, double **ts_f, int MCMCiter, int N,Priors *priors, unsigned long *seed, char *filel, char *file2, double propvar[], Hyper_priors *hyper)
void mcmc(Patient *patientlist, PopulationPriors *popprior, PopulationPriors *popprior_response,  PopulationEstimates *popparms, PopulationProposals *pv_pop, PopulationPopulationEstimates *popparms_response, PopulationProposals *pv_pop_resp, AssocPriors *assocprior, AssocEstimates *assocparms, AssocProposals *pv_assoc, int Nsubj, unsigned long *seed, int MCMCinter,int Nthin, int Nburnin)
{
    int i;
    
//Run the MCMC for Niter iterations
	for (i = 0; i < MCMCinter; i++) {

	
        //Why do I need iteration number assigned to variables
//		parms_f->iter = i;
//		parms_l->iter = i;


//Draw population mean parameters for mean pulse mass and width: mua and muw
		draw_pop_mass_mean(patientlist, popparms, popparms_response, popprior, popprior_response, poprior, popprior_response, pv_assoc, Nubj,seed);   //for trigger and response jointly

        draw_pop_width_mean(patientlist, popparms, popprior, pv_pop, Nsubj, seed); //for trigger: This should be done with Pop Model function

// Draw patient-to-patient variance in pulse mass mean and SD of width
        
		draw_pop_mass_mean_var(sublist, priors, parms_f, seed, hyper);  //this is unique to joint model
        
        draw_pop_width_mean_SD(sublist, priors, vfewv_l, vfewv_f, seed, hyper); //for trigger: This should be done with Pop Model function

// Draw population mean baseline and halflife; mub and muh*/

		draw_pop_base_mean(sublist, priors, parms_f, seed, hyper);  //do for trigger
        draw_pop_hl_mean(sublist, priors, parms_f, seed, hyper);
		draw_pop_base_mean(sublist, priors, parms_l, seed, hyper);  //do for response
        draw_pop_hl_mean(sublist, priors, parms_f, seed, hyper);
        
// draw patient-to-patient SD in baseline and halflife; sb and sh  //why drawn together??

		draw_pop_base_SD(sublist, priors, vfebv_l, vfehv_l, seed, hyper); //do for trigger
        draw_pop_hl_SD(sublist, priors, vfebv_l, vfehv_l, seed, hyper);
		draw_pop_base_SD(sublist, priors, vfebv_f, vfehv_f, seed, hyper);  //do for response
        draw_pop_hl_SD(sublist, priors, vfebv_l, vfehv_l, seed, hyper);

// draw subject specific means; muak and muwk*/

        draw_patient_mass_means(sublist, priors, parms_l, parms_f, seed, pmean_var); //jt draw for trigger&response
        draw_patient_width_mean(sublist, priors, parms_l, vfepw_l, seed);  //for trigger
        draw_patient_width_mean(sublist, priors, parms_f, vfepw_f, seed); //for response
		
// draw pulse-to-pulse mass and width SD's; sa and sw*/

		draw_patient_mass_SD(sublist, priors, parms_l, vfepmv_l, vfepwv_l, seed); //for trigger
		draw_patient_mass_SD(sublist, priors, parms_f, vfepmv_f, vfepwv_f, seed); //for response

		birth_death_trigger(sublist, ts_l, parms_l, N, seed, i);

// draw trigger pulse times
        draw_times(sublist, parms_l, ts_l,N, seed, vtime_l);

// draw individual pulse masses and widths for trigger
        draw_pulse_mass(ts_l, sublist, parms_l, N, vrem_l, vrew_l, seed);
        draw_pulse_width(ts_l, sublist, parms_l, N, vrem_l, vrew_l, seed);)
        draw_tscalemass(sublist, parms_l, seed, veta_l);
        draw_tscalewidth(sublist, parms_l, seed, veta_l);
			
// draw patient level baseline and width for trigger
        draw_patient_base(sublist, parms_l, priors, ts_l,N, seed, pmd_var_l);
		draw_patient_hl(sublist, parms_l, priors, ts_l,N, seed, pmd_var_l);
	

// draw pulse locations for reponse; tauki*/
		birth_death_response(sublist, ts_f, parms_f, parms_l, N, seed, i);

//draw times for response
        
		draw_times(sublist, parms_f, parms_l, ts_f,
			N, seed, vtime_f);
		draw_tscalemass(sublist, parms_f, seed, veta_f);
        draw_tscalewidth(sublist, parms_f, seed, veta_f);
	
// draw individual pulse masses and widths for the response
		draw_pulse_mass(ts_f, sublist, parms_f, N, vrem_f, vrew_f, seed);
        draw_pulse_width(ts_f, sublist, parms_f, N, vrem_f, vrew_f, seed);
        
//draw patient level baseline and halflife
		draw_patient_base(sublist, parms_f, priors, ts_f,
			N, seed, pmd_var_f);
        draw_patient_hl(sublist, parms_f, priors, ts_f,
                          N, seed, pmd_var_f);

//draw cox model parameters

		
		mh_logclustersize(sublist, parms_l, priors, seed, vrho);
		mh_logclusterwidth(sublist, parms_l, priors, seed, vnu);


		/* draw model error; s2e*/
		ssq_l = error_squared(ts_l, sublist, parms_l, N);  //for trigger
		parms_l->sigma = inverse_gamma(priors->alpha_l + (parms_l->numsub*N) / 2, priors->beta_l + 0.5*ssq_l, seed);
		parms_l->lsigma = log(parms_l->sigma);
		ssq_f = error_squared(ts_f, sublist, parms_f, N); //for response
		parms_f->sigma = inverse_gamma(priors->alpha_f + (parms_f->numsub*N) / 2, priors->beta_f + 0.5*ssq_f, seed);
		parms_f->lsigma = log(parms_f->sigma);

		fflush(stdout);
        
//STOPPED HERE
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

