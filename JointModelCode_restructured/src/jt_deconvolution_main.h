#pragma once
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>


/*****This structure contains the information of each pulse********/
typedef struct {
    double pv_time;
    double pv_mass;
    double pv_width;
    double pv_tscalemass;
    double pv_tscalewidth;
    
} PulseProposals;

typedef struct node_tag{
	struct node_tag *succ;
	struct node_tag *pred;
	double time;  /*tau_ki*/
	double mass;  /*theta^s_ki*/
    double width; /*theta^w_ki*/
    double tvarscalemass;  /*variance scale for mass t-dist*/
    double tvarscalewidth;   /*variance scal for width t-dist*/
	double *mean_contrib;
	double lambda; /*for fsh pulse only, lambda, denomsum
	                */  /***NOT SURE WHAT THIS TERM IS FOR***/
    PulseProposals *pulsepv;
} PulseEstimate;

/****This structure contains the patient level parameters***/
typedef struct {
    
    Node_type *pulselist;
    double baseline;
    double halflife;
    double decay;     // decay rate converted from above half-life
    double errorsq;    // model error (variance)
    double logerrorsq; // log of model error (may not be used)
    double mass_mean;
    double width_mean;
    int pulse_count;
    
} PatientEstimates;

/****This structure contains the patient level proposal variances***/
typedef struct {
    double pv_baseline;
    double pv_halflife;
    double pv_mass_mean;
    double pv_width_mean;
} PatientProposals;

typedef struct {
    
    double *concentration;
    double *time;
    double *response_concentration;
    double *response_time;
    int number_of_obs;
    double avg_period_of_obs;  /**in minutes**/
    double duration_of_obs = number_of_obs * period_of_obs;  /**in minutes**/
    
    char *common_filename; /*patient level filename for output for baseline, half-life, mass_mean, width_mean*/
    FILE *csubfile;   /*patient level file for baseline, half-life, mass_mean, width_mean*/
    char *pulse_filename; /*patient level filename for output for pulse level chains, mass, width, locations */
    FILE *psubfile; /*patient level file for above*/
    char *resp_common_filename; /*response hormone filename for output for baseline, half-life, mass_mean, width_mean*/
    FILE *resp_csub_filename;  /*response hormone patient level file for output for baseline, half-life, mass_mean, etc*/
    char *resp_pulse_filename; /* resp hormone filename for pulse level output: mass, width, locations*/
    FILE *resp_psub_filename; /*resp hormone file for pulse level output*/

} PatientData;

/****Thus structure contains the patient level data and parameters for both trigger and response hormones***/
typedef struct patient_tag {
    struct patient_tag *succ;
    struct patient_tag *pred;
    PatientEstimates *patient_estimates;
    PatientEstimates *resp_patient_estimates;
    PatientData *patient_data;
    PatientProposals *patient_pv;
    
} Patient;

/*******the population level parameters.  These are also in a Bayesian framework, priors on the patient level parameters*******/
typedef struct {
    double pv_baseline_mean;
    double pv_HL_mean;
    double pv_baseline_variance;
    double pv_HL_variance;
    double pv_mass_mean;
    double pv_width_mean;
    double pv_mass_SD;
    double pv_width_SD;
    double pv_mass_mean_SD;
    double pv_width_mean_SD;
} PopulationProposal;


typedef struct{
    double baseline_mean;  /*theta_b*/
    double halflife_mean;  /*theta_h*/
    
    double baseline_variance;  /*patient to patient variation of baselines*/
    double halflife_variance;  /*patient to patient variation of half-lifes*/
    
	double mass_mean;     /*population mean mass: mean of the patient level means*/
    double width_mean;    /*population mean width: mean of the patient level means*/

    double mass_SD;  /* population level pulse-to-pulse variation in pulse mass*/
    double width_SD; /* population level pulse-to-pulse variation in pulse width*/

    double mass_mean_SD; /*patient-to-patient SD in mean pulse mass*/
    double mass_mean_variance; /**?  **/

    double width_mean_SD; /*patient-to-patient SD in mean pulse width*/
    
    PopulationProposal *pv_population;
    char *pop_filename; /*filename of output for population level parameters*/
    FILE *popfile; /*The file containing the population level parameters*/
    
    
} PopulationEstimates;

typedef struct{
    double pv_clustersize;
    double pv_clusterwidth;
    double pv_masscorr;
} AssocProposals;

/****This structure contains the priors on the patient level info: population parameters***/
typedef struct{
    double cluster_size; /*cluster size (rho) in the cox process for the response hormone intensity*/
    double log_cluster_size; /*log scale cluster size*/
    double cluster_width; /*cluster width (nu) in the cox process for the respose hormone intensity*/
    double log_cluster_width;
    
    AssocProposals *pv_assoc;
    
    char *popassoc_filename; /* output filename of the association parameters*/
    FILE *popassocfile; /* file of the association parameters*/
    
} AssocEstimates;

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
    double corr_alpha;  /*alpha parameter in the beta prior on the correlation between the pulse masses*/
    double corr_beta;   /*beta parameter in the beta prior on the correlation between the pulse masses*/
} AssocPriors;



