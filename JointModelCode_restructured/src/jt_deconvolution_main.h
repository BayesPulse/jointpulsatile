#pragma once
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

/*****This structure contains the information of each pulse********/
typedef struct node_tag{
	struct node_tag *succ;
	struct node_tag *pred;
	double time;  /*tau_ki*/
	double mass;  /*theta^s_ki*/
    double width; /*theta^w_ki*/
    double kappamass;
    double kappawidth;
	double *mean_contrib;
	double lambda; /*for fsh pulse only, lambda, denomsum
	                */  /***NOT SURE WHAT THIS TERM IS FOR***/
} PulseEstimate;

/****This structure contains the patient level parameters***/
typedef struct subject_tag {
    struct subject_tag *succ;
    struct subject_tag *pred;
    Node_type *pulselist;
    double baseline;
    double halflife;
    double decay;     // decay rate converted from above half-life
    double errorsq;    // model error (variance)
    double logerrorsq; // log of model error (may not be used)
    double mass_mean;
    double width_mean;
    double mass_sd;
    double width_sd;
    int pulse_count;
    
    /*Question: should we have the data in this structure or just the parameters*/
    /*Question: if no data should we also move the output files to the data structure*/

    char *common; /*patient level output for baseline, half-life, mass_mean, width_mean*/
    char *pulse; /*patient level output for pulse level chains, mass, width, locations */
    FILE *csub;
    FILE *psub;
} PatientEstimates;

/*******the prior distribution is collasping all information of fsh and lh*******/

typedef struct{
    double baseline_mean;  /*theta_b*/
    double halflife_mean;  /*theta_h*/
    
    double baseline_variance;  /*patient to patient variation of baselines*/
    double halflife_variance;  /*patient to patient variation of half-lifes*/
    
	double mass_mean;     /*population mean mass: mean of the patient level means*/
    double width_mean;    /*population mean width: mean of the patient level means*/

	double **fe_precision; /*the variance matrix of correlated pulse mass on the population level */

	double **re_var; /*inverse of fe_precision*/

	double fe_precision_wf;/*sigma_w，drawe_fe_priors2*/
	double fe_precision_wl;/*sigma_w，drawe_fe_priors2*/

	double *re_var_f; /*specified, the half-cauchy or uniform prior for subject-level variance */
	double *re_var_l;
	double alpha_f;  /* prior parameters for model error,fsh */
	double beta_f;
	double alpha_l;  /* prior parameters for model error,lh */
	double beta_l;
    double rho_prior[2]; /*log scale*/
	double nu_prior[2];  /*log scale*/

} Priors;

/*the hyper is also defined together for both the driver and response */
typedef struct {
	double hmean_f[2]; /*ma,mw*/
	double hmean_l[2]; /*ma,mw*/

	double hvar_f; /*vw*/
	double hvar_l; /*vw*/
	
	double precision_wf;
	double precision_wl;
	double **sig_a_inv;

	double **prec;  /* prior variance-covariance of pulse mass,,F-1*/
	double meanmeanbh_l[2];  /*mhn,prior mean fo baseline and half life*/
	double meanvarbh_l[2]; /*prior variance for baseline and halflife */
	double varbh_l[2]; /*prior for variance for baseline and halflife*/
	double meanmeanbh_f[2];  /*mhn,prior mean fo baseline and half life*/
	double meanvarbh_f[2]; /*prior variance for baseline and halflife */
	double varbh_f[2]; /*prior for variance for baseline and halflife*/

} Hyper_priors;


/**INCORRECT---INCOMPLETE EDITING
/****This structure contains the priors on the patient level info: population parameters***/
typedef struct{
    
    double Rinter1; /*(hard core)*/
    double Rinter2; /*struass*/
    double *re_precision;  /*nu_a,nu_w,draw_fe_precision same for each subject*/
    int numsub;
    int subindex;
    int iter;
    /*	double nprior;*/
    double rho; /*cluster size, used for lh only,on the log scale*/
    double lrho;
    double lnu;
    double nu; /*cluster width, used for lh only,it is the square not standard deviation*/
    double beta; /*intensity function for strauss process*/
    double gamma; /*repulsion function for struass */
} PopulationEstimates;
