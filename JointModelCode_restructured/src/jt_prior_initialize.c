//
//  jt_prior_initialize.c
//  
//
//  Created by Nichole on 1/22/18.
//
//

#include "jt_prior_initialize.h"

/**********************************************************************
 MAIN PROGRAM STARTS HERE
 
 VARIABLE DEFINITIONS
 int *nobs: Number of observations for each subject
 int MCMCiter: Number of MCMC interations to perform
 int nsubj: Number of subjects in the analysis.  This is the number of columns in the dataset minus 2 (for obs number and time of observation).
 int Nthin: how to thin the output to reduce file storage requirements and autocorrelation in the observations used to summarize the results.
 int Nburnin: how much burnin prior to outputting the results.
 
 
 char datafile: Name of the rectangular file with the trigger hormone series.  One patient for each column.
 char datafile_response: Name of the rectangular file with the response hormone series. One patient for each column. Columns must line up between the trigger and response.
 char pulse_fnameroot: root of the filename for output for pulse specific parameters.
 char patient_fnameroot: root of the filename for output for patient specific parameters.
 char pop_fnameroot: root of the filename for output for population level parameters.
 char assoc_fname: filename for output for parameters with the association components.
 
 double PopMeanMassPrior_input: User defined mean in the prior on the population level mean pulse mass for the hormone. On natural scale since using truncated normal distributions for this prior
 double PopMeanMassPriorVar_input: User defined variance on the prior for the population level mean pulse mass for hormone.
 double PulseMassSDMax_input: User defined maximum on the prior distribution for the pulse-to-pulse variation in the pulse masses of the hormone. Pulse masses are on the natural scale.
 double PatientMeanMassSDMax_input: User defined maximum on the maximum on the prior
 
 double PopMeanWidthPrior_input: User defined mean in the prior on the population level mean pulse Width for the hormone. On natural scale since using truncated normal distributions for this prior
 double PopMeanWidthPriorVar_input: User defined variance on the prior for the population level mean pulse Width for the hormone.
 double PulseWidthSDMax_input: User defined maximum on the prior distribution for the pulse-to-pulse variation in the pulse Widthes of the hormone. Pulse Widthes are on the natural scale.
 double PatientMeanWidthSDMax_input: User defined maximum on the maximum on the prior
 
 double strauss_range1_input: The prior on the trigger pulse locations is a multilevel Stauss.  The first component is a Hard Core, this is the range of the hard core (strict repulsion) in hours.
 double strauss_range2_input: The 2nd component is a strauss with repulsion.  This is the range of the repulsion in hours.
 double straussrate_input: Rate of the number of pulses per hour on the Strauss.
 double straussrepulsion_input; Strength of repulsion.  Needs to be between 0 and 1.  0 is strict repulsion (hard core), 1 is no repulsion.
 
 double PopMeanBasePrior_input: User defined mean on the prior for the population mean baseline
 double PopMeanBasePriorVar_input: User defined variance on the prior for the population mean baseline
 double PatientBaseSDMax_input: user defined maximum of the SD in the Unif prior on the patient-to-patient variation in baseline
 
 double PopMeanHLPrior_input: User defined mean on the prior for the population mean halflife
 double PopMeanHLPriorVar_input: User defined variance on the prior for the population mean halflife
 double PatientHLSDMax_input: User defined maximum on the SD in the Unif prior on the patient-to-patient variance in half-life
 
 double alpha_input: User defined input in into the Gamma for the inverse of the model error variance for each subject
 dobule beta_input: User defined input into the Gamma prior for the invese of the model error variance for each subject
 
 
 ************************************************************************/


void prior_initialize(File *finput,PopulationPriors *popprior,PopulationPriors *popprior_response) {
    
/*****************************************/
/***Variable defining***/
/*****************************************/
    double PopMeanMassPrior_input, PopMeanMassPriorVar_input, PulseMassSDMax_input, PatientMeanMassSDMax_input; //user inputs to define trigger prior info for population mean mass and corresponding SD's
    double PopMeanWidthPrior_input, PopMeanWidthPriorVar_input, PulseWidthSDMax_input, PatientMeanWidthSDMax_input; //user inputs to define trigger prior info for population mean Width and corresponding SD's
    double strauss_range1_input, strauss_range2_input, straussrate_input, straussrepulsion_input; //user inputs to define the prior on the trigger pulse locations
    double PopMeanBasePrior_input, PopMeanBasePriorVar_input, PatientBaseSDMax_input; //user inputs to define trigger prior into for the popultion mean baseline and corresponding SD's
    double PopMeanHLPrior_input, PopMeanHLPriorVar_input, PatientHLSDMax_input; //user inputs to define trigger prior into for the popultion mean halflife and corresponding SD's
    double alpha_input, beta_input; //user inputs into the Gamma prior on the inverse of the model error for each subject

    
/*Receive user input for the priors for the pulse masses (mean and variance of prior on mean and max on priors on SDs) and then set them in the Population Priors data structure*/
    fscanf(finput,"%lf %lf %lf %lf\n", &PopMeanMassPrior_input, &PopMeanMassPriorVar_input, &PulseMassSDMax_input, &PatientMeanMassSDMax_input);
    popprior->mass_mean = PopMeanMassPrior_input;
    popprior->mass_variance = PopMeanMassPriorVar_input;
    popprior->mass_mean_SD_max = PatientMeanMassSDMax_input;
    popprior->mass_SD_max = PulseMassSDMax_input;
    
    
    fscanf(finput,"%lf %lf %lf %lf\n", &PopMeanMassPrior_input, &PopMeanMassPriorVar_input, &PulseMassSDMax_input, &PatientMeanMassSDMax_input);
    popprior_response->mass_mean = PopMeanMassPrior_input;
    popprior_response->mass_variance = PopMeanMassPriorVar_input;
    popprior_response->mass_mean_SD_max = PatientMeanMassSDMax_input;
    popprior_response->mass_SD_max = PulseMassSDMax_input;
    
    
/*Receive user input for the priors for the pulse widths (mean and variance of prior on mean and max on priors on SDs) and then set them in the Population Priors data structures */
    fscanf(finput,"%lf %lf %lf %lf\n", &PopMeanWidthPrior_input, &PopMeanWidthPriorVar_input, &PulseWidthSDMax_input, &PatientWidthWidthSDMax_input);
    popprior->width_mean = PopMeanWidthPrior_input;
    popprior->width_variance = PopMeanWidthPriorVar_input;
    popprior->width_mean_SD_max = PatientMeanWidthSDMax_input;
    popprior->width_SD_max = PulseWidthSDMax_input;
    
    
    fscanf(finput,"%lf %lf %lf %lf\n", &PopMeanWidthPrior_input, &PopMeanWidthPriorVar_input, &PulseWidthSDMax_input, &PatientMeanWdithSDMax_input);
    popprior_response->width_mean = PopMeanWidthPrior_input;
    popprior_response->width_variance = PopMeanWidthPriorVar_input;
    popprior_response->width_mean_SD_max = PatientMeanWidthSDMax_input;
    popprior_response->width_SD_max = PulseWidthSDMax_input;

    
/*Receive user input for the priors for the pulse locations of the trigger pulse location model*/
    fscanf(finput, "%lf %lf %lf %lf\n", &strauss_range1_input, &strauss_range2_input, &straussrate_input, &straussrepulsion_input);
    popprior->strauss_range1 = strauss_range1_input;
    popprior->strauss_range2 = strauss_range2_input;
    popprior->StraussRate = straussrate_input
    popprior->StraussRepulsion = straussrepulsion_input;
    
/*For good practice set response parameters to zero.  These parameters are not used for the response hormone*/
    popprior_response->strauss_range1 = strauss_range1_input;
    popprior_response->strauss_range2 = strauss_range2_input;
    popprior_response->StraussRate = straussrate_input
    popprior_response->StraussRepulsion = straussrepulsion_input;
    
    
/*Receive user input for the priors for baseline (mean and variance of prior on mean and max on priors on SD) and then set them in the Population Priors data structures */
    fscanf(finput,"%lf %lf %lf\n",PopMeanBasePrior_input, PopMeanBasePriorVar_input, PatientBaseSDMax_input);
    popprior->baseline_mean = PopMeanBasePrior_input;
    popprior->baseline_variance = PopMeanBasePriorVar_input;
    popprior->baseline_SD_max = PatientBaseSDMax_input;
    
    fscanf(finput,"%lf %lf %lf\n",PopMeanBasePrior_input, PopMeanBasePriorVar_input, PatientBaseSDMax_input);
    popprior_response->baseline_mean = PopMeanBasePrior_input;
    popprior_response->baseline_variance = PopMeanBasePriorVar_input;
    popprior_response->baseline_SD_max = PatientBaseSDMax_input;
    

/*Receive user input for the priors for half-life (mean and variance of prior on mean and max on priors on SD) and then set them in the Population Priors data structures */
    fscanf(finput,"%lf %lf %lf\n",PopMeanHLPrior_input, PopMeanHLPriorVar_input, PatientHLSDMax_input);
    popprior->HLline_mean = PopMeanHLPrior_input;
    popprior->HLline_variance = PopMeanHLPriorVar_input;
    popprior->HLline_SD_max = PatientHLSDMax_input;
    
    fscanf(finput,"%lf %lf %lf\n",PopMeanHLPrior_input, PopMeanHLPriorVar_input, PatientHLSDMax_input);
    popprior_response->HLline_mean = PopMeanHLPrior_input;
    popprior_response->HLline_variance = PopMeanHLPriorVar_input;
    popprior_response->HLline_SD_max = PatientHLSDMax_input;
    
/*Receive user input for the priors for the model error and set them in the Popultion Priors data structures */
    fscanf(finput,"%lf %lf\n",&alpha_input, &beta_input);
    popprior->alpha = alpha_input;
    popprior->beta = beta_input;
    
    fscanf(finput,"%lf %lf\n",&alpha_input, &beta_input);
    popprior_response->alpha = alpha_input;
    popprior_response->beta = beta_input;
    
}
