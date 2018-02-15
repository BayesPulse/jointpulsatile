/*******************************************************************/
/*************************draw_pop_width_mean.c*********************/
/*******************************************************************/

#include "draw_pop_width_mean.h"

/*********************************************************************/
/*START OF draw_pop_width_mean SUBROUTINE                             */
/*********************************************************************/
/*********************************************************************/
/*draw_pop_widht_mean: this is a metropolis-hastings draw of the both population mean pulse masses
 /* The subroutine only exists for the joint hormone models
 ARGUMENTS:  Patient *patientlist,
 PopulationEstimate *popparms,
 PopulationPrior *popprior,
 PopulationProposal *pv_pop
 
 PopulationProposal *pv_assoc,
 int Nsubj,
 unsigned long *seed
 RETURNS: None; all updates are made internally
 Internal update will be for Population Mean Pulse for both hormones.
 *********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 TO BE FILLED IN
 
 SUBROUTINES USED
 
 TO BE FILLED IN
 **************************************************************************/
void draw_pop_width_mean(Patient *patientlist, PopulationEstimate *popparms, PopulationPrior *popprior, PopulationProposal *pv_pop, int Nsubj, unsigned long *seed)
    {
        /* declare variables */
        
        Patient *patient;
        
        int i, j;
        double sum_mass, proposalmean, sum_numerator, sum_denominator;
        double sum_norm_int_numerator, sum_norm_int_denominator, alpha, temp;
        
        //Draw the proposed new pulse width mean
        proposalmean = rnorm(popparms->width_mean, pv_pop->pv_width_mean, seed);
        
        //Increment the denominator of the acceptance probability;
        pv_pop->N_width_mean++;
        
        //Calculate the proposal ratio for the MH.  If the draw is <0 for either parameter the likelihood would be zero
        // so reject the proposal and do not waste computation time
        
        if (proposalmean>0) {
            
            //Compute the ratio of "likelihoods"
            
            //sum exponent in the normal likelihood: distribution of the patient level pulse mass means
            sum_numerator = 0;
            sum_denominator = 0;
            sum_norm_int_numerator = 0;
            sum_norm_int_denominator = 0;
            
            patient = patientlist->succ;
            
            /*This while loop goes through the pulses to get the sum involved*/
            while (patient != NULL){
                sum_numerator += -0.5 * (patient->patient_parms->width_mean - proposalmean) * (patient->patient_parms->width_mean - proposalmean)/ (popparms->width_SD * popparms->width_SD);
                sum_norm_int_numerator += log(1 - phi(0,proposalmean,popparms->width_SD));
                
                sum_denominator += -0.5 * (patient->patient_parms->width_mean - popparms->width_mean) * (patient->patient_parms->width_mean - popparms->width_mean)/ (popparms->width_SD * popparms->width_SD);
                
                sum_norm_int_denominator += log(1 - phi(0,popparms->width_mean,popparms->width_SD));
                
                patient = patientlist->succ;
            } /*end of loop through pulses*/
            
            likelihood_ratio = (sum_numerator - sum_denominator) - sum_norm_int_numerator + sum_norm_int_denominator;
            
            //Compute the prior_ratio.
            
            prior_numerator = =0.5 *(proposalmean - popprior->width_mean) * (proposalmean - popprior->width_mean)/(popprior->width_variance);
            

            prior_denomonator = -0.5 * (popparms->width_mean - popprior->width_mean) * (popparms->width_mean - popprior->width_mean)/ (popprior->width_variance);
            
            prior_ratio = prior_numerator - prior_denomonator;  //note the integrals in the denomonators cancel
            //also note that the proposal ratios cancel because we are doing a random walk with a Gaussian proposal distribution
            
            
            //MH decision here.
            alpha = (0 < (temp = (likelihood_ratio + prior_ratio))) ? 0 : temp;
            
            /*If log(U) < log rho, accept the proposed value*/
            /*Increase acceptance count by 1*/
            if (log(kiss(seed)) < alpha){
                popparms->width_mean = proposalmean;
                pv_pop->Naccept_width_mean++;
                
            }
        }
    }

