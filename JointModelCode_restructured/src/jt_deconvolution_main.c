/*******************************************************************/
/********************MAIN BDMCMC PROGRAM***************************/
/*******************************************************************/

#include "jt_deconvolution_main.h"

/*******************************************************************
*******************GLOBAL VARIABLE DEFINITIONS**********************

 * M: Sets the environment precision and used for all random number generation
    * Needed specfically for the KISS subroutine
 * mmm: Order statistic used for distribution of pulse locations.
    * This is inputted by the user and is typically 3
 * fitstart: The first time in hours that a pulse may occur
 * fitend: The last time in hours that a pulse may occur

*********************************************************************/

double M;
//int mmm;
double fitstart;
double fitend;

#if _WIN32 || _WIN64
#if _WIN64
#define ENVIRONMENT 64
#else
#define ENVIRONMENT 32
#endif
#endif

// Check GCC
#if __GNUC__
#if __x86_64__ || __ppc64__
#define ENVIRONMENT 64
#else
#define ENVIRONMENT 32
#endif
#endif

/********************************************************************/
/*SUBROUTINES THAT EXIST IN THIS PROGRAM

 **********************************************************************/

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
    char fnameroot: root of the filename for output.
 
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
    double beta_input: User defined input into the Gamma prior for the invese of the model error variance for each subject
 
 ***cox process/association parameters****
    double mean_cluster_input: mean on the prior of the cluster size.  On log scale.
    double var_cluster_input: variance on the prior of the cluster size. On log scale squared.
    double mean_clusterwidth_input: mean on the prior of the cluster width. on log scale.
    double var_clusterwidth_input: variance on the prior of the cluster width. on log scale.
 double mass_corr_beta_mean_input: mean of the beta prior on the correlation of the pulse masses.
 double mass_corr_beta_var_input: variance of the beta prior on the correlation of the pulse masses.
 double corr_prior_alpha: alpha parameter in the beta prior on the correlation.
 double corr_prior_beta: beta parameter in the beta prior on the correlation.
 
 
 
 ***Parameters at population level and starting values****
   double svmean, svptsd, svsd;
 PopulationEstimates *popparms;
 PopulationEstimates *popparms_response;
 PopulationProposal *pv_pop;
 PopulationProposal *pv_pop_response;
 double pvmean, pvptsd, pvsd;
 double svmean, svptsd, svsd;
 
 double sv_ptmass_input:
 double sv_ptwidth_input:
 double sv_ptbase_input:
 double sv_pthl_input:
 double sv_pt_modelerrorvar:
 double pv_pt_mass_input:
 double pv_pt_width_input:
 double pv_pt_base_input:
 double pv_pt_hl_input:
 
 AssocPriors *assocprior;
 AssocEstimates *assocparms;
 AssocProposals *pv_assoc;
 
 double sv_cluster_input:
 double sv_width_input:
 double sv_mass_corr_input:
 
 double pv_time_input,pv_mass_input,pv_width_input,pv_tscalemass_input, pv_tscalewidth_input;
************************************************************************/

int main(int argc,char *argv[])
{
    
/*****************************************/
/***Variable defining***/
/*****************************************/
    int *Nobs,Nsubj, MCMCiter,Nthin, Nburnin, i;
    
    char datafile[90];
    char datafile_response[100];
    char fnameroot[100];
    
    FILE *fseed, *finput;
    
    unsigned long *seed;
    
    PopulationPriors *popprior, *popprior_response;  //Data structures that hold the population prior settings from the user for trigger and response, respectively
  
    double PopMeanMassPrior_input, PopMeanMassPriorVar_input, PulseMassSDMax_input, PatientMeanMassSDMax_input; //user inputs to define trigger prior info for population mean mass and corresponding SD's
    double PopMeanWidthPrior_input, PopMeanWidthPriorVar_input, PulseWidthSDMax_input, PatientMeanWidthSDMax_input; //user inputs to define trigger prior info for population mean Width and corresponding SD's
    double strauss_range1_input, strauss_range2_input, straussrate_input, straussrepulsion_input; //user inputs to define the prior on the trigger pulse locations
    double PopMeanBasePrior_input, PopMeanBasePriorVar_input, PatientBaseSDMax_input; //user inputs to define trigger prior into for the popultion mean baseline and corresponding SD's
    double PopMeanHLPrior_input, PopMeanHLPriorVar_input, PatientHLSDMax_input; //user inputs to define trigger prior into for the popultion mean halflife and corresponding SD's
    double alpha_input, beta_input; //user inputs into the Gamma prior on the inverse of the model error for each subject
    
    double mean_cluster_input, var_cluster_input; //user inputs for mean and variance on the cluster size
    double mean_clusterwidth_input, var_clusterwidth_input; //user inputs for mean and variance on the cluster width
    double mass_corr_beta_mean_input, mass_corr_beta_var_input, corr_prior_alpha, corr_prior_beta; //user inputs for beta prior on correlation of the patient level mean pulse masses
    double svmean, svptsd, svsd; //user inputs for starting values of population level parameters
    
    AssocPriors *assocprior;
    AssocEstimates *assocparms;
    AssocProposals *pv_assoc;
    
    double sv_cluster_input, sv_width_input, sv_mass_corr_input; //user input starting value for association parameters
    double pv_cluster_size, pv_cluster_width, pv_mass_corr; //user input starting value for proposal variances for association parameters
    
    PopulationEstimates *popparms;
    PopulationEstimates *popparms_response;
    PopulationProposal *pv_pop;
    PopulationProposal *pv_pop_response;
    double pvmean, pvptsd, pvsd;
    double svmean, svptsd, svsd;
    
    Patient *patientlist, *patient;
    double sv_ptmass_input,sv_ptwidth_input,sv_ptbase_input,sv_pthl_input,sv_pt_modelerrorvar;
    double pv_pt_mass_input,pv_pt_width_input,pv_pt_base_input,pv_pt_hl_input;
    double sv_ptmassresp_input,sv_ptwidthresp_input,sv_ptbaseresp_input,sv_pthlresp_input,sv_pt_respmodelerrorvar;
    double pv_pt_massresp_input,pv_pt_widthresp_input,pv_pt_baseresp_input,pv_pt_hlresp_input;
    
    double pv_time_input,pv_mass_input,pv_width_input,pv_tscalemass_input, pv_tscalewidth_input; //user input for starting values of the proposal variance for the individual pulses.
    double pv_resptime_input,pv_respmass_input,pv_respwidth_input,pv_resptscalemass_input, pv_resptscalewidth_input; //user input for starting values of the proposal variance for the individual resp pulses.
    
/*This assesses if there is an input filename with the program call*/

    if (argc != 2 ) {
        printf("There is no inputfile given.\n");
        exit(0);
    }

/************************************/
/* M is used in the KISS random number generator */
    M = exp(-ENVIRONMENT*log(2.0));
/*************************************************************/
    
/*Read in the random number generator seed*/
/*This might go away with R interface*/
    seed = (unsigned long *)calloc(3,sizeof(unsigned long));
    fseed = fopen("seed.dat","r");
    fscanf(fseed,"%lu %lu %lu\n",&seed[0],&seed[1],&seed[2]);
    fclose(fseed);

/*Open the input file and read in the information*/
    finput = fopen(argv[1],"r");

    fscanf(finput,"%s %s %s\n", datafile, datafile_response, fnameroot);  /*First line: datafile names for trigger and respones*/
    fscanf(finput,"%d %d %d %d\n", &Nsubj, &MCMCiter, &Nthin, &Nburnin);  /*third line: number of subject, number of iteration*/

/* read in the hormonal time series */
    Nobs = (int *)calloc(1,sizeof(int)); // define a dynamic variable for the number of observations on a subject
    
    ts = read_data_file(datafile,Nobs,nsubj); // time series of the trigger hormone concentration.
    ts_response = read_data_file(datafile_response,Nobs,nsubj); //time series of the response hormone concentration.
    fitend = ts[*Nobs-1][0]+ ts[1][0] * 4;  /*search 2 units farther in time: these set the boundaries for looking for pulse locations.*/
    fitstart = -ts[1][0] * 4;  /*search 4 units in the past*/

/**    mmm = 3;  Not needed in this program.  We don't use an order statistic for the trigger hormone pulse location model**/
    
/*Receive user input for the priors for the pulse masses (mean and variance of prior on mean and max on priors on SDs) and then set them in the Population Priors data structure*/
    
    popprior = calloc(1,sizeof(PopulationPriors));
    popprior_response = calloc(1,sizeof(PopulationPriors));
    
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
    popprior->StraussRate = straussrate_input;
    popprior->StraussRepulsion = straussrepulsion_input;
    
/*For good practice set response parameters to zero.  These parameters are not used for the response hormone*/
    popprior_response->strauss_range1 = 0;
    popprior_response->strauss_range2 = 0;
    popprior_response->StraussRate = 0;
    popprior_response->StraussRepulsion = 0;
    
    
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
    
    
/**Read in the association parameters ***/
    
    assocprior = calloc(1,sizeof(AssocPriors));
    
    fscanf(finput,"%lf %lf\n",&mean_cluster_input, &var_cluster_input);
    fscanf(finput,"%lf %lf\n",&mean_clusterwidth_input, &var_clusterwidth_input);
    fscanf(finput,"%lf %lf\n",&mass_corr_alpha_input, &mass_corr_beta_input);
    
    assocprior->mean_log_cluster_size = log(mean_cluster_input);
    assocprior->variance_log_cluster_size = log(var_cluster_input);
    assocprior->mean_log_cluster_width = log(mean_clusterwidth_input);
    assocprior->variance_log_cluster_width = log(var_clusterwidth_input);
    assocprior->corr_alpha = mass_corr_alpha_input;
    associprior->corr_beta = mass_corr_beta_input;
    
/**Read in the starting values for population parameters**/
    
    popparms = calloc(1,sizeof(PopulationPriors));
    popparms_response = calloc(1,sizeof(PopulationPriors));
    
    fscanf(finput,"%lf %lf %lf\n",&svmean,&svptsd,&svsd); // read in the mass starting values
    popparms->mass_mean = svmean;
    popparms->mass_SD = svptsd;
    popparms->mass_mean_SD = svsd;
    
    fscanf(finput,"%lf %lf %lf\n",&svmean,&svptsd,&svsd); // read in the mass starting values for the response
    popparms_response->mass_mean = svmean;
    popparms_response->mass_SD = svptsd;
    popparms_response->mass_mean_SD = svsd;
    
    fscanf(finput,"%lf %lf %lf\n",&svmean,&svptsd,&svsd); // read in the width starting values
    popparms->width_mean = svmean;
    popparms->width_SD = svptsd;
    popparms->width_mean_SD = svsd;
    
    fscanf(finput,"%lf %lf %lf\n",&svmean,&svptsd,&svsd); // read in the width starting values for the response
    popparms_response->width_mean = svmean;
    popparms_response->width_SD = svptsd;
    popparms_response->width_mean_SD = svsd;
    
    fscanf(finput,"%lf %lf\n",&svmean,&svsd); //read in the baseline starting values
    popparms->baseline_mean = svmean;
    popparms->baseline_SD = svsd;
    
    fscanf(finput,"%lf %lf\n",&svmean,&svsd); // read in the baseline starting values for the response
    popparms_response->baseline_mean = svmean;
    popparms_response->baseline_SD = svsd;
    
    fscanf(finput,"%lf %lf\n",&svmean,&svsd); //read in the half-life starting values
    popparms->halflife_mean = svmean;
    popparms->halflife_SD = svsd;
    
    fscanf(finput,"%lf %lf\n",&svmean,&svsd); // read in the half-life starting values for the response
    popparms_response->halflife_mean = svmean;
    popparms_response->halflife_SD = svsd;
    
    //name the output file for trigger
    popparms->pop_filename = (char *)calloc(40,sizeof(char *));
    strcpy(popparms->pop_filename,fnameroot);
    strcat(popparms->pop_filename,"_pop_trig.out");
    popparms->popfile = fopen(popparms->pop_filename,"w");
    
    //name the output file for response
    popparms_response->pop_filename = (char *)calloc(40,sizeof(char *));
    strcpy(popparms_response->pop_filename,fnameroot);
    strcat(popparms_response->pop_filename,"_pop_resp.out");
    popparms_response->popfile = fopen(popparms_response->pop_filename,"w");
    
/***Read in the starting values for the association parameters***/
    
    assocparms = calloc(1,sizeof(AssocEstimates));
    
    fscanf(finput,"%lf %lf %lf\n",sv_cluster_input,sv_width_input,sv_mass_corr_input);
    assocparms->cluster_size = sv_cluster_input;
    assocparms->log_cluster_size = log(sv_cluster_input);
    assocparms->cluster_width = sv_width_input;
    assocparms->log_cluster_width = log(sv_width_input);
    assocparms->mass_corr = sv_mass_corr_input;
    
    pv_assoc = calloc(1,sizeof(AssocProposals));
    
    fscanf(finput,"%lf %lf %lf\n",pv_cluster_size, pv_cluster_width, pv_mass_corr);
    //name the output file for association parameters
    assocparms->popassoc_filename = (char *)calloc(40,sizeof(char *));
    strcpy(assocparms->popassoc_filename,fnameroot);
    strcat(assocparms->popassoc_filename,"_pop_assoc.out");
    assocparms->popfile = fopen(assocparms->popassoc_filename,"w");
    

/**Read in the starting values for patient level parameters**/
    //we need to create the patient level data info
    //we need to create the patient level estimates for trigger and response
    
    fscanf(finput,"%lf %lf %lf %lf %lf\n",sv_ptmass_input,sv_ptwidth_input,sv_ptbase_input,sv_pthl_input,sv_pt_modelerrorvar);
    fscanf(finput "%lf %lf %lf %lf\n",pv_pt_mass_input,pv_pt_width_input,pv_pt_base_input,pv_pt_hl_input);
    fscanf(finput "%lf %lf %lf %lf %lf\n",pv_time_input,pv_mass_input,pv_width_input,pv_tscalemass_input, pv_tscalewidth_input);
    
    fscanf(finput,"%lf %lf %lf %lf %lf\n",sv_ptmassresp_input,sv_ptwidthresp_input,sv_ptbaseresp_input,sv_pthlresp_input,sv_pt_respmodelerrorvar);
    fscanf(finput "%lf %lf %lf %lf\n",pv_pt_massresp_input,pv_pt_widthresp_input,pv_pt_baseresp_input,pv_pt_hlresp_input);
    fscanf(finput "%lf %lf %lf %lf %lf\n",pv_resptime_input,pv_respmass_input,pv_respwidth_input,pv_resptscalemass_input, pv_resptscalewidth_input);
    
    patientlist = initialize_subject();
    
    i=0;
    for (i=0;i<Nsubj;i++) {
        patient=initialize_subject();
        
        //fill patient data: create filenames, set data
        patient->patient_data->common_filename = (char *)calloc(40,sizeof(char *));
        patient->patient_data->pulse_filename = (char *)calloc(40,sizeof(char *));
        patient->patient_data->resp_common_filename = (char *)calloc(40,sizeof(char *));
        patient->patient_data->resp_pulse_filename = (char *)calloc(40,sizeof(char *));

        strcpy(patient->patient_data->common_filename,fnameroot);
        strcat(patient->patient_data->common_filename,"_ptparms_trig");
        strcat(patient->patient_data->common_filename,".out");
        
        strcpy(patient->patient_data->resp_common_filename,fnameroot);
        strcat(patient->patient_data->resp_common_filename,"_ptparms_resp");
        strcat(patient->patient_data->resp_common_filename,".out");

        sprintf(patientnumb,"%d",Nsubj-i);
        
        strcpy(subject->common_l,commonl);
        strcat(subject->common_l,"s");
        strcat(subject->common_l,tmp);
        strcat(subject->common_l,".out");
        strcpy(patient->patient_data->pulse_filename,fnameroot);
        strcat(patient->patient_data->pulse_filename,"_pulse_trig");
        strcat(patient->patient_data->pulse_filename,patientnumb);
        strcat(patient->patient_data->pulse_filename,".out");
        
        strcpy(patient->patient_data->resp_pulse_filename,fnameroot);
        strcat(patient->patient_data->resp_pulse_filename,"_pulse_resp");
        strcat(patient->patient_data->resp_pulse_filename,".out");

        patient->patient_data->number_of_obs = Nobs;
        patient->patient_data->avg_period_of_obs = ts[1][0];
        
        patient->patient_data->concentration = calloc(1,sizeof(Nobs));
        patient->patient_data->time = calloc(1,sizeof(Nobs));
        patient->patient_data->resp_concentration = calloc(1,sizeof(Nobs));
        
        for (j=0;j<Nobs;j++) {
            patient->patient_data->concentration[j] = ts[(i+1),j];
            patient->patient_data->time[j] = ts[0,j];
            patient->patient_data->resp_concentration = ts_response[(i+1),j];
        }
        
        //set starting values for patient estimates
        patient->patient_parms->mass_mean = sv_ptmass_input;
        patient->patient_parms->width_mean = sv_ptwidth_input;
        patient->patient_parms->baseline = sv_ptbase_input;
        patient->patient_parms->halflife = sv_pthl_input;
        patient->patient_parms->errorsq = sv_pt_modelerrorvar;
        
        patient->resp_patient_parms->mass_mean = sv_ptmassresp_input;
        patient->resp_patient_parms->width_mean = sv_ptwidthresp_input;
        patient->resp_patient_parms->baseline = sv_ptbaseresp_input;
        patient->resp_patient_parms->halflife = sv_pthlresp_input;
        patient->resp_ patient_parms->errorsq = sv_pt_respmodelerrorvar;
        
        //set initial proposal variances
        patient->patient_pv->pv_mass_mean = pv_pt_mass_input;
        patient->patient_pv->pv_width_mean = pv_pt_width_input;
        patient->patient_pv->pv_baseline = pv_pt_base_input;
        patient->patient_pv->pv_halflife = pv_pt_hl_input;
        
        patient->resp_patient_pv->pv_mass_mean = pv_pt_massresp_input;
        patient->resp_patient_pv->pv_width_mean = pv_pt_widthresp_input;
        patient->resp_patient_pv->pv_baseline = pv_pt_baseresp_input;
        patient->resp_patient_pv->pv_halflife = pv_pt_hlresp_input;
        
        //set initial proposal variances for pulse parameters
        patient->pulse_pv->pv_time = pv_time_input;
        patient->pulse_pv->pv_mass = pv_mass_input;
        patient->pulse_pv->pv_width = pv_width_input;
        patient->pulse_pv->pv_tscalemass = pv_tscalemass_input;
        patient->pulse_pv->pv_tscalewidth = pv_tscalewidth_input;
        
        patient->pulse_pv_response->pv_time = pv_resptime_input;
        patient->pulse_pv_response->pv_mass = pv_respmass_input;
        patient->pulse_pv_response->pv_width = pv_respwidth_input;
        patient->pulse_pv_response->pv_tscalemass = pv_resptscalemass_input;
        patient->pulse_pv_response->pv_tscalewidth = pv_resptscalewidth_input;
        
        insert_subject(patient,patientlist);
    }

    
    pv_pop = calloc(1,sizeof(PopulationPriors));
    pv_pop_response = calloc(1,sizeof(PopulationPriors));
    
    fscanf(finput,"%lf %lf %lf\n",&pvmean,&pvptsd,&pvsd); // read in the mass starting proposal variances
    pv_pop->pv_mass_mean = pvmean;
    pv_pop->pv_mass_SD = pvptsd;
    pv_pop->pv_mass_mean_SD = pvsd;
    
    fscanf(finput,"%lf %lf %lf\n",&pvmean,&pvptsd,&pvsd); // read in the mass starting proposal variances for the response
    pv_pop_response->pv_mass_mean = pvmean;
    pv_pop_response->pv_mass_SD = pvptsd;
    pv_pop_response->pv_mass_mean_SD = pvsd;
    
    fscanf(finput,"%lf %lf %lf\n",&pvmean,&pvptsd,&pvsd); // read in the width starting proposal variance
    pv_pop->pv_width_mean = pvmean;
    pv_pop->pv_width_SD = pvptsd;
    pv_pop->pv_width_mean_SD = pvsd;
    
    fscanf(finput,"%lf %lf %lf\n",&pvmean,&pvptsd,&pvsd); // read in the width starting proposal variance for the response
    pv_pop_response->pv_width_mean = pvmean;
    pv_pop_response->pv_width_SD = pvptsd;
    pv_pop_response->pv_width_mean_SD = pvsd;
    
    fscanf(finput,"%lf %lf\n",&pvmean,&pvsd); //read in the baseline starting values
    pv_pop->pv_baseline_mean = pvmean;
    pv_pop->pv_baseline_SD = pvsd;
    
    fscanf(finput,"%lf %lf\n",&pvmean,&pvsd); // read in the baseline starting values for the response
    pv_pop_response->pv_baseline_mean = pvmean;
    pv_pop_response->pv_baseline_SD = pvsd;
    
    fscanf(finput,"%lf %lf\n",&pvmean,&pvsd); //read in the half-life starting values
    pv_pop->pv_HL_mean = pvmean;
    pv_pop->pv_HL_SD = pvsd;
    
    fscanf(finput,"%lf %lf\n",&pvmean,&pvsd); // read in the half-life starting values for the response
    pv_pop_response->pv_HL_mean = pvmean;
    pv_pop_response->pv_HL_SD = pvsd;
    
    
 /***********saved for future use **************/
/*********************************************/ 
  /*input prior for rho and nu here*/
/*************************************************/
    fscanf(finput,"%lf %lf %lf %lf\n", &mstart1, &mstart2, &mstart3,&mstart4); /*10th line, starting value for mean pulse features*/
    priors->fe_mean_l[0] = mstart1;
	priors->fe_mean_f[0]=mstart2;
	priors->fe_mean_l[1] = mstart3;
	priors->fe_mean_f[1]=mstart4;
  
    
    fscanf(finput,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &pstart1, &pstart2, &pstart3,&pstart4,&pstart5,&pstart6,&pstart7,&pstart8,&pstart9 );/*11th line starting value for variance pulse features*/
    /*first three is the variance-covariance for population variance and cor of pulse mass*/
	/*next two (4-5) is the population variance for pulse width */
	/*next four (6-9) is the subject-level variance for pulse mass and width*/
	priors->re_var[0][0] = pstart1;
	priors->re_var[1][1] = pstart2;
	priors->re_var[0][1]=priors->re_var[1][0] = sqrt(pstart2)*sqrt(pstart1)*pstart3;
	  if (!cholesky_decomp(priors->re_var,2)) {
     printf("proposal variance matrix not PSD matrix\n");
     exit(0);
   }
   priors->fe_precision = cholesky_invert(2,priors->re_var);
   priors->fe_precision_wl = pstart4;
   priors->fe_precision_wf = pstart5;  
   parms_l->re_precision[0]= pstart6;
   parms_l->re_precision[1]= pstart7;
   parms_f->re_precision[0]= pstart8;
   parms_f->re_precision[1]= pstart9;

   
   

   
    fscanf(finput,"%lf %lf %lf %lf\n", &bstart1, &bstart2,&bstart3, &bstart4); /*12th line baseline*/
    priors->meanbh_l[0] = bstart1;
    priors->varbh_l[0] = bstart2;
	priors->meanbh_f[0] = bstart3;
    priors->varbh_f[0] = bstart4;
    
    fscanf(finput,"%lf %lf %lf %lf\n", &hstart1, &hstart2,&hstart3, &hstart4); /*13th line half life */
    priors->meanbh_l[1] = hstart1;
    priors->varbh_l[1] = hstart2;
	priors->meanbh_f[1] = hstart3;
    priors->varbh_f[1] = hstart4;

    fscanf(finput,"%lf %lf\n", &sstart1, &sstart2);  /*14th line error*/
    parms_l->sigma = sstart1;
    parms_l->lsigma = log(parms_l->sigma);
	parms_f->sigma = sstart2;
    parms_f->lsigma = log(parms_f->sigma);
	parms_l -> rho= 1;
    parms_l ->lrho = log(parms_l->rho);
	
	parms_l ->nu= 4;
	parms_l -> lnu=log(parms_l ->nu);
	parms_l->gamma = 0.1;
	parms_l->beta =4;
	


    parms_l->numsub =parms_f->numsub = subj;

    fscanf(finput,"%lf %lf \n", &propvar[0], &propvar[1]);  /*15th lineï¼?prop for population width variance */
    fscanf(finput,"%lf %lf %lf %lf\n", &propvar[2],&propvar[3], &propvar[4], &propvar[5]); /*16 th line: prop for four variance term of population md*/
    fscanf(finput,"%lf %lf %lf\n", &propvar[6], &propvar[7],&propvar[8]); /*17 th line: prop for subject mass, third is the cor*/
    fscanf(finput,"%lf %lf\n", &propvar[9], &propvar[10]); /*18th line: uwi*/
    fscanf(finput,"%lf %lf %lf %lf \n", &propvar[11],&propvar[12],&propvar[13],&propvar[14]); /*19th line v*/
	fscanf(finput,"%lf %lf %lf %lf %lf %lf\n", &propvar[15],&propvar[16],&propvar[17],&propvar[18],&propvar[19],&propvar[20]); /*20th line bi, hi,3rd and 6th is the cor*/
    fscanf(finput,"%lf %lf %lf %lf\n", &propvar[21],&propvar[22],&propvar[23],&propvar[24]); /*21th line v, random effect, lh first and fsh*/
    fscanf(finput,"%lf %lf \n", &propvar[25],&propvar[26]); /*22th line v, time, lh first and fsh*/
	fscanf(finput, "%lf %lf %lf %lf  \n", &propvar[27], &propvar[28], &propvar[29], &propvar[30]); /*23th line v, eta, lh first and fsh*/
	fscanf(finput, "%lf %lf \n", &propvar[31], &propvar[32]); /*24th line rho and niu*/

    fclose(finput);

     /*Create the list of subjects*/
      sublist = initialize_subject();
          
      i=0;
      for(i=0;i<subj;i++){
          subject=initialize_subject();
          subject->common_l = (char *)calloc(30,sizeof(char *));
          subject->pulse_l = (char *)calloc(30,sizeof(char *));
          subject->common_f = (char *)calloc(30,sizeof(char *));
          subject->pulse_f = (char *)calloc(30,sizeof(char *));

          sprintf(tmp,"%d",subj-i);

          strcpy(subject->common_l,commonl);
          strcat(subject->common_l,"s");
          strcat(subject->common_l,tmp);
          strcat(subject->common_l,".out");
		  
		  strcpy(subject->common_f,commonf);
          strcat(subject->common_f,"s");
          strcat(subject->common_f,tmp);
          strcat(subject->common_f,".out");
          
          strcpy(subject->pulse_l,parml);
          strcat(subject->pulse_l,"s");
          strcat(subject->pulse_l,tmp);
          strcat(subject->pulse_l,".out");

		  strcpy(subject->pulse_f,parmf);
          strcat(subject->pulse_f,"s");
          strcat(subject->pulse_f,tmp);
          strcat(subject->pulse_f,".out");
		  
          subject->basehalf_l[0]=bstart1;
		  subject->basehalf_f[0]=bstart3;
          subject->basehalf_l[1]=hstart1;
		  subject->basehalf_l[1]=hstart3;


          subject->theta_l[0]=mstart1;
		  subject->theta_f[0]=mstart2;
          subject->theta_l[1]=mstart3;
		  subject->theta_f[1]=mstart4;

          insert_subject(subject,sublist);
      }
  
      fflush(stdout);

   //   mcmc(sublist,parms_l,parms_f,ts_l,ts_f,iter,*N,priors,seed,commonl,commonf,propvar,hyper);  //not correct now
  
    destroy_sublist(sublist);
    /**************************/
    
    /* save the current random number as the seed for the next simulation */
    fseed = fopen("seed.dat","w");
    fprintf(fseed,"%lu %lu %lu\n",seed[0],seed[1],seed[2]);
    fclose(fseed);
    /**********************************************************************/

    /* deallocate resources */

    free(seed);
    for (i=0;i<*N;i++)
      free(ts_f[i]);
    free(ts_f);
	 for (i=0;i<*N;i++)
      free(ts_l[i]);
    free(ts_l);
	
    free(N);
    free(priors->re_var_f);
    free(priors->re_var_l);
	free(priors->fe_mean_l);
    free(priors->fe_mean_f);
	free(parms_f->re_precision);
    free(parms_l->re_precision);
	
	 for (i=0;i<2;i++)
      free(hyper->prec[i]);
    free(hyper->prec);

	 for (i=0;i<2;i++)
      free(priors->fe_precision[i]);
    free(priors->fe_precision);
	
	for (i=0;i<2;i++)
    free(priors->re_var[i]);
 free(priors->re_var);
	

   

    free(hyper);
    free(priors);
    free(parms_l);
	free(parms_f);
	free(temp2);
/************************/

  return 0;
}


