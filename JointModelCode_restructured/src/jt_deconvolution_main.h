#pragma once
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <imsls.h>


/*****This structure contains the information of each pulse********/
typedef struct node_tag{
	struct node_tag *succ;
	struct node_tag *pred;
	double time;  /*pulse location in minutes*/
	double mass;  /*pulse mass on natural scale*/
    double width; /*pulse width: variance of the normal distribution of the pulse shape*/
    double tvarscalemass;  /*variance scale for mass t-dist*/
    double tvarscalewidth;   /*variance scal for width t-dist*/
	double *mean_contrib; /*pulse level contribution to the integral of mean concentration */
	//double lambda; /*for fsh pulse only, lambda, denomsum*/  /***NOT SURE WHAT THIS TERM IS FOR***/
} PulseEstimates;

//Proposal variances and MH acceptance probability components for pulse specific parameter.
//These will go away with integration with pop model package.
typedef struct {
    double pv_time;
    double pv_mass;
    double pv_width;
    double pv_tscalemass;
    double pv_tscalewidth;
    
    int Naccept_pulse_mass;
    int N_pulse_mass;
    int Naccept_pulse_width;
    int N_pulse_width;
    int Naccept_pulse_location;
    int N_pulse_location;
    int Naccept_tscalemass;
    int N_tscalemass;
    int Naccept_tscalewidth;
    int N_tscalewidth;
    
} PulseProposals;

/****This structure contains the patient level parameters***/
typedef struct {
    
    PulseEstimates *pulselist;
    double baseline;
    double halflife;
    double decay;     // decay rate converted from above half-life
    double errorsq;    // model error (variance)
    double logerrorsq; // log of model error (may not be used)
    double mass_mean;
    double width_mean;
    int pulse_count;
    
} PatientEstimates;

//Proposal variances and MH acceptance probability components for patient specific parameter.
//These will go away with integration with pop model package.
typedef struct {
    double pv_baseline;
    double pv_halflife;
    double pv_mass_mean;
    double pv_width_mean;
    
    int Naccept_baseline;
    int N_baseline;
    int Naccept_halflife;
    int N_halflife;
    int Naccept_mass_mean;
    int N_mass_mean;
    int Naccept_width_mean;
    int N_width_mean;
    
} PatientProposals;

typedef struct {
    
    double *concentration;
    double *time;
    double *response_concentration;
    int *number_of_obs;
    double avg_period_of_obs;  /**in minutes**/
    double duration_of_obs;  /**in minutes**/
    
    char *common_filename; /*patient level filename for output for baseline, half-life, mass_mean, width_mean*/
    FILE *csubfile;   /*patient level file for baseline, half-life, mass_mean, width_mean*/
    char *pulse_filename; /*patient level filename for output for pulse level chains, mass, width, locations */
    FILE *psubfile; /*patient level file for above*/
    char *resp_common_filename; /*response hormone filename for output for baseline, half-life, mass_mean, width_mean*/
    FILE *resp_csub_file;  /*response hormone patient level file for output for baseline, half-life, mass_mean, etc*/
    char *resp_pulse_filename; /* resp hormone filename for pulse level output: mass, width, locations*/
    FILE *resp_psub_file; /*resp hormone file for pulse level output*/

} PatientData;

/****Thus structure contains the patient level data and parameters for both trigger and response hormones***/
typedef struct patient_tag {
    struct patient_tag *succ;
    struct patient_tag *pred;
   
    PatientEstimates *patient_parms;
    PatientEstimates *resp_patient_parms;
    PatientData *patient_data;
    PatientProposals *patient_pv;
    PatientProposals *resp_patient_pv;
    PulseProposals *pulse_pv;
    PulseProposals *resp_pulse_pv;
    
} Patient;

/*******the population level parameters.  These are also in a Bayesian framework, priors on the patient level parameters*******/
typedef struct{
    double baseline_mean;  /*theta_b*/
    double halflife_mean;  /*theta_h*/
    
    double baseline_SD;  /*patient to patient variation of baselines*/
    double halflife_SD;  /*patient to patient variation of half-lifes*/
    
	double mass_mean;     /*population mean mass: mean of the patient level means*/
    double width_mean;    /*population mean width: mean of the patient level means*/

    double mass_SD;  /* population level pulse-to-pulse variation in pulse mass*/
    double width_SD; /* population level pulse-to-pulse variation in pulse width*/

    double mass_mean_SD; /*patient-to-patient SD in mean pulse mass*/
    double width_mean_SD; /*patient-to-patient SD in mean pulse width*/
    
    char *pop_filename; /*filename of output for population level parameters*/
    FILE *popfile; /*The file containing the population level parameters*/
    
    
} PopulationEstimates;

//Proposal variances and MH acceptance probability components.  These will go away with integration with pop model package.
typedef struct {
    double pv_baseline_mean;
    double pv_HL_mean;
    double pv_baseline_SD;
    double pv_HL_SD;
    double pv_mass_mean;
    double pv_width_mean;
    double pv_mass_SD;
    double pv_width_SD;
    double pv_mass_mean_SD;
    double pv_width_mean_SD;
    
    int Naccept_baseline_mean;
    int N_baseline_mean;
    int Naccept_HL_mean;
    int N_HL_mean;
    int Naccept_baseline_SD;
    int N_baseline_SD;
    int Naccept_HL_SD;
    int N_HL_SD;
    int Naccept_mass_mean;
    int N_mass_mean;
    int Naccept_width_mean;
    int N_width_mean;
    int Naccept_mass_SD;
    int N_mass_SD;
    int Naccept_width_SD;
    int N_width_SD;
    int Naccept_mass_mean_SD;
    int N_mass_mean_SD;
    int Naccept_width_mean_SD;
    int N_width_mean_SD;
    
} PopulationProposal;

/****This structure contains the priors on the patient level info: population parameters***/
typedef struct{
    double cluster_size; /*cluster size (rho) in the cox process for the response hormone intensity*/
    double log_cluster_size; /*log scale cluster size*/
    double cluster_width; /*cluster width (nu) in the cox process for the respose hormone intensity*/
    double log_cluster_width;
    double mass_corr;
    double popprior_mass_corr;

    char *popassoc_filename; /* output filename of the association parameters*/
    FILE *popassocfile; /* file of the association parameters*/
    
} AssocEstimates;

typedef struct{
    double pv_clustersize;
    double pv_clusterwidth;
    double pv_masscorr;
    double pv_massmatrix[2][2]; //This is proposal for the population pulse masses
    
    int Naccept_clustersize;
    int N_clustersize;
    int Naccept_clusterwidth;
    int N_clusterwidth;
    int Naccept_masscorr;
    int N_masscorr;
    
} AssocProposals;

/*The user defined values for the priors on the population level parameters*/
/*This structure contains all the variables that the user sets when setting the priors*/
typedef struct {
    double mass_mean;  // mean of the prior on the population mean pulse mass
    double mass_variance; // variance of the prior on the population mean pulse mass
    double mass_SD_max;  // maximum of the Unif distribution for the pulse-to-pulse SD of pulse masses within a patient
    double mass_mean_SD_max; // maximum of the Unif distributio for the patient-to-patient SD of patient level mena pulse masses
    
    double width_mean;
    double width_variance;
    double width_SD_max;
    double width_mean_SD_max;
    
    double strauss_range1; /*(hard core): range of interaction*/
    double strauss_range2; /*strauss range of interactio*/
    double StraussRate; /*intensity function for strauss process: beta*/
    double StraussRepulsion; /*repulsion function for struass: gamma */

    double baseline_mean;
    double baseline_variance;
    double baseline_SD_max;
    
    double halflife_mean;
    double halflife_variance;
    double halflife_SD_max;
    
    double alpha;  /* prior parameter in gamma for model error */
    double beta;
 

} PopulationPriors;

/****This structure contains the priors on the patient level info: population parameters***/
typedef struct{
    double mean_log_cluster_size; /*mean of prior on cluster size (rho) in the cox process for the response hormone intensity*/
    double variance_log_cluster_size; /*variance in prior on the cluster size*/
    double mean_log_cluster_width; /*mean on prior of cluster width (nu) in the cox process for the respose hormone intensity*/
    double variance_log_cluster_width; /*variance in prior on the cluster width*/
    double var_popmeanmass_df;  /*df in the wishart prior on the mean pulse masses*/
    double var_popmeanmass[2][2];   /*matrix in the prior on the wishart for the mean pulse masses*/
} AssocPriors;

//Function declaration
void mcmc(Patient *patientlist, PopulationPriors *popprior, PopulationPriors *popprior_response,  PopulationEstimates *popparms, PopulationProposals *pv_pop, PopulationPopulationEstimates *popparms_response, PopulationProposals *pv_pop_resp, AssocPriors *assocprior, AssocEstimates *assocparms, AssocProposals *pv_assoc, int Nsubj, unsigned long *seed, int MCMCinter,int Nthin, int Nburnin);

