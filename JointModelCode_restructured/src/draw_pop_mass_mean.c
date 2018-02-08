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

void draw_pop_mass_mean(Patient *patientlist, PopulationEstimate *popparms, PopulationEstimate *popparms_response, PopulationPrior *popprior, PopulationPrior *popprior_response, AssocEstimates *assocparms, PopulationProposal *pv_assoc, int Nsubj, unsigned long *seed);
//void draw_fe_prior_a_mean(Subject_type *sublist, Priors *priors, Common_parms *parms_f, unsigned long *seed, Hyper_priors *hyper)
{
	/* declare variables */

	Patient *patient;

	int i, j;
	double sum_mass[2], *current_means, proposalmeans[2], **patient_mass_var, **patient_mass_var_inv;
    double sum_norm_int_numerator, sum_norm_int_denominator, **prior_mass_var, **prior_mass_var_inv;
 
	patient_mass_var = (double **)calloc(2, sizeof(double *));
	patient_mass_var_inv = (double **)calloc(2, sizeof(double *));
	current_means = (double *)calloc(2, sizeof(double));
	for (i = 0; i < 2; i++) {
		patient_mass_var[i] = (double *)calloc(2, sizeof(double));
        patient_mass_var_inv[i] = (double *)calloc(2, sizeof(double));
	}
    
//Draw the proposal
    current_means[0] = popparms->mass_mean;
    current_means[1] = popparms_response->mass_mean;
    rmvnorm(proposalmeans, pv_assoc->massmatrix, 2, current_means, seed, 1);

    if (proposalmeans[0]>0 & proposalmeans[1]>0) {
        
//Compute the ratio of "likelihoods"
    
        //Step 1: set current patient level mass matrix: it is stored as SD's and a correlation for estimation
        patient_mass_var[0][0] = popparms->mass_SD * popparms->mass_SD;
        patient_mass_var[1][1] = popparms_response->mass_SD * popparms_response->mass_SD;
        patient_mass_var[0][1] = patient_mass_var[1][0] = assocparms->mass_corr * popparms->mass_SD * popparms_response->mass_SD;
   
        //Step 2: invert the patient level mass matrix.
        if (!cholesky_decomp(patient_mass_var,2)) {
            printf("prior variance matrix not PSD matrix\n");
            exit(0);
        }
        patient_mass_var_inv = cholesky_invert(2, patient_mass_var);
    
        //sum exponent in the bivariate normal likelihood: distribution of the patient level pulse mass means
        sum_numerator = 0;
        sum_denominator = 0;
        sum_norm_int_numerator = 0;
        sum_norm_int_denominator = 0;

        patient = patientlist->succ;

        /*This while loop goes through the pulses to get the sum involved*/
        while (patient != NULL){
            sum_numerator += (patient->patient_parms->mass_mean - proposalmeans[0]) * patient_mass_var_inv[0][0] * (patient->patient_parms->mass_mean - proposalmeans[0]);
            sum_numerator += (patient->resp_patient_parms->mass_mean - proposalmeans[1]) * patient_mass_var_inv[1][0] * (patient->patient_parms->mass_mean - proposalmeans[0]);
            sum_numerator += (patient->resp_patient_parms->mass_mean - proposalmeans[1]) * patient_mass_var_inv[0][1] * (patient->patient_parms->mass_mean - proposalmeans[0]);
            sum_numerator += (patient->resp_patient_parms->mass_mean - proposalmeans[1]) * patient_mass_var_inv[1][1] * (patient->resp_patient_parms->mass_mean - proposalmeans[1]);
            sum_norm_int_numerator += NEED log MVN CDF;
        
            sum_denominator += (patient->patient_parms->mass_mean - popparms->mass_mean) * patient_mass_var_inv[0][0] * (patient->patient_parms->mass_mean - popparms->mass_mean);
            sum_denominator += (patient->resp_patient_parms->mass_mean - popparms_response->mass_mean) * patient_mass_var_inv[1][0] * (patient->patient_parms->mass_mean - popparms->mass_mean);
            sum_denominator += (patient->resp_patient_parms->mass_mean - popparms_response->mass_mean) * patient_mass_var_inv[0][1] * (patient->patient_parms->mass_mean - popparms->mass_mean);
            sum_denominator += (patient->resp_patient_parms->mass_mean - popparms_response->mass_mean) * patient_mass_var_inv[1][1] * (patient->resp_patient_parms->mass_mean - popparms_response->mass_mean);
            sum_norm_int_denominator += NEED log MNV CDF;
        
            patient = patientlist->succ;
        } /*end of loop through pulses*/
	
        likelihood_ratio = 0.5 * sum_numerator - sum_norm_int_numerator - 0.5 * sum_denominator + sum_norm_int_denominator;
        
//Compute the prior_ratio.
        //Set current population level mass matrix: it is stored as SD's and a correlation for estimation
        //NOT FIXED HERE. WHAT ABOUT CORRELATED PRIOR.
        patient_mass_var[0][0] = popparms->mass_SD * popparms->mass_SD;
        patient_mass_var[1][1] = popparms_response->mass_SD * popparms_response->mass_SD;
        patient_mass_var[0][1] = patient_mass_var[1][0] = assocparms->mass_corr * popparms->mass_SD * popparms_response->mass_SD;
        
        //Step 2: invert the patient level mass matrix.
        if (!cholesky_decomp(patient_mass_var,2)) {
            printf("prior variance matrix not PSD matrix\n");
            exit(0);
        }
        patient_mass_var_inv = cholesky_invert(2, patient_mass_var);
        prior_numerator = 0;
        prior_denominator = 0;
        prior_int_numerator = 0;
        prior_int_denominator = 0;
        
        prior_numerator += (proposalmeans[0] - popprior->mass_mean) * patient_mass_var_inv[0][0] * (patient->patient_parms->mass_mean - proposalmeans[0]);
        prior_numerator += (patient->resp_patient_parms->mass_mean - proposalmeans[1]) * patient_mass_var_inv[1][0] * (patient->patient_parms->mass_mean - proposalmeans[0]);
        sum_numerator += (patient->resp_patient_parms->mass_mean - proposalmeans[1]) * patient_mass_var_inv[0][1] * (patient->patient_parms->mass_mean - proposalmeans[0]);
        sum_numerator += (patient->resp_patient_parms->mass_mean - proposalmeans[1]) * patient_mass_var_inv[1][1] * (patient->resp_patient_parms->mass_mean - proposalmeans[1]);
        sum_norm_int_numerator += NEED log MVN CDF;
        
        sum_denominator += (patient->patient_parms->mass_mean - popparms->mass_mean) * patient_mass_var_inv[0][0] * (patient->patient_parms->mass_mean - popparms->mass_mean);
        sum_denominator += (patient->resp_patient_parms->mass_mean - popparms_response->mass_mean) * patient_mass_var_inv[1][0] * (patient->patient_parms->mass_mean - popparms->mass_mean);
        sum_denominator += (patient->resp_patient_parms->mass_mean - popparms_response->mass_mean) * patient_mass_var_inv[0][1] * (patient->patient_parms->mass_mean - popparms->mass_mean);
        sum_denominator += (patient->resp_patient_parms->mass_mean - popparms_response->mass_mean) * patient_mass_var_inv[1][1] * (patient->resp_patient_parms->mass_mean - popparms_response->mass_mean);
        prior_ratio = PUT IN;

//NEED TO ADD MH DECISION LOOP HERE
        
    }
    for (i = 0; i < 2; i++) {
        free(fcond_var[i]);
        free(fcond_var_inv[i]);
    }
    free(fcond_var);
	free(fcond_var_inv);
	free(cond_mean);
}
