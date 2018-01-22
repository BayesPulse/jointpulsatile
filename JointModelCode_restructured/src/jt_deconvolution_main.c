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

int main(int argc,char *argv[])
{
    
/*****************************************/
/***Variable defining***/
/*****************************************/
    int *nobs,Nsubj, MCMCiter,Nthin, Nburnin;
    
    char datafile[90];
    char datafile_response[100];

    char pulse_fnameroot[100];
    char patient_fnameroot[100];
    char pop_fnameroot[100];
    char assoc_fname[100];

    double PopMeanMassPrior_input, PopMeanMassPriorVar_input, PulseMassSDMax_input, PatientMeanMassSDMax_input; //user inputs to define trigger prior info for population mean mass and corresponding SD's
    double PopMeanWidthPrior_input, PopMeanWidthPriorVar_input, PulseWidthSDMax_input, PatientMeanWidthSDMax_input; //user inputs to define trigger prior info for population mean Width and corresponding SD's
    double strauss_range1_input, strauss_range2_input, straussrate_input, straussrepulsion_input; //user inputs to define the prior on the trigger pulse locations
    double PopMeanBasePrior_input, PopMeanBasePriorVar_input, PatientBaseSDMax_input; //user inputs to define trigger prior into for the popultion mean baseline and corresponding SD's
    double PopMeanHLPrior_input, PopMeanHLPriorVar_input, PatientHLSDMax_input; //user inputs to define trigger prior into for the popultion mean halflife and corresponding SD's
    double alpha_input, beta_input; //user inputs into the Gamma prior on the inverse of the model error for each subject

    
    char tmp[5];
  unsigned long *seed;
  double **ts_l,**ts_f,propvar[33] ,**temp2;
 double mprior1, mprior2, mprior3, mprior4;
  double psprior1, psprior2, psprior3, psprior4,psprior5, psprior6, psprior7, psprior8,psprior9;
  double bprior1, bprior2, bprior3, bprior4, bprior5, bprior6;
  double hprior1, hprior2, hprior3,hprior4, hprior5, hprior6;
  double sprior1, sprior2,sprior3,sprior4;
  /*int nprior1, nprior2;*/
  double gprior,betaprior;
  double rhoprior1,rhoprior2;
  double nuprior1,nuprior2;
  
  double mstart1, mstart2, mstart3, mstart4;
  double pstart1, pstart2, pstart3,pstart4, pstart5, pstart6,pstart7, pstart8, pstart9,psprior10, psprior11;
  double bstart1, bstart2,bstart3, bstart4;
  double hstart1, hstart2,hstart3, hstart4;
  double sstart1,sstart2;
  double rhostart,nustart;

  FILE *fseed, *finput;
  Common_parms *parms_l,*parms_f;
  Priors *priors;
  Hyper_priors *hyper;
  Subject_type *subject, *sublist;

  double **read_data_file(char *,int *,int);
  double rnorm(double,double,unsigned long *);
  double rgamma(double,double,unsigned long *);
  void destroy_sublist(Subject_type *);
  Subject_type *initialize_subject(void);
  void insert_subject(Subject_type *,Subject_type *);
  void mcmc(Subject_type *, Common_parms *, Common_parms *, double **, double **,long, int, Priors *, unsigned long *,
	  char *,  char *, double[], Hyper_priors *);
  char *itoa(int , char *, int );
  int cholesky_decomp(double **,int);
  double **cholesky_invert(int, double **);



  Subject_type *initialize_subject(void);
    
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

    fscanf(finput,"%s %s \n", datafile, datafile_response);  /*First line: data for lh and fsh*/
    fscanf(finput,"%s %s %s %s \n", pulse_fnameroot, patient_fnameroot, pop_fnameroot, assoc_fname); /*second line: filename roots for the output files-pulse level, patient level, population level and association */
 
    fscanf(finput,"%d %d \n", &Nsubj, &MCMCiter, &Nthin, &Nburnin);  /*third line: number of subject, number of iteration*/

/* read in the hormonal time series */
    Nobs = (int *)calloc(1,sizeof(int)); // define a dynamic variable for the number of observations on a subject
    
    ts = read_data_file(datafile,Nobs,nsubj); // time series of the trigger hormone concentration.
    ts_response = read_data_file(datafile_response,Nobs,nsubj); //time series of the response hormone concentration.
    fitend = ts[*Nobs-1][0]+ ts[1][0] * 4;  /*search 2 units farther in time: these set the boundaries for looking for pulse locations.*/
    fitstart = -ts[1][0] * 4;  /*search 4 units in the past*/

/**    mmm = 3;  Not needed in this program.  We don't use an order statistic for the trigger hormone pulse location model**/
    
/*Set up prior structure for the population parameters*/
    
    popprior = (PopulationPriors *)calloc(1,sizeof(PopulationPriors));
    popprior_response = (PopulationPriors *)calloc(1,sizeof(PopulationPriors));
    
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

    

    


    fscanf(finput,"%lf %lf %lf %lf\n", &mprior1, &mprior2, &mprior3, &mprior4); /*fourth line: lh pulse Width mean, fsh pulse mass mean, lh pulse width mean, fsh pulse width mean*/
    hyper->hmean_l[0] = mprior1;
    hyper->hmean_f[0] = mprior2;
    hyper->hmean_l[1] = mprior3;
    hyper->hmean_f[1] = mprior4;

	fscanf(finput, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &psprior1, &psprior2, &psprior3, &psprior4, &psprior5, &psprior6, &psprior7, &psprior8, &psprior9,&psprior10, &psprior11); /*5th line: first 3 is var-cor for pulse mass variance, 4-5 is for pulse width */
    
	temp2[0][0] = psprior1;
    temp2[1][1]= psprior2;
	temp2[1][0]= temp2[0][1]= sqrt(psprior2)*sqrt(psprior1)*psprior3;
		  if (!cholesky_decomp(temp2,2)) {
     printf("prior variance matrix not PSD matrix\n");
     exit(0);
   }
      hyper->prec = cholesky_invert(2, temp2);

    hyper->precision_wl=psprior4; /*prior for population variance of pulse width */
    hyper->precision_wf=psprior5;

	hyper->hvar_l=psprior6;    /*variance term for population mean of puse width*/
	hyper->hvar_f=psprior7;
	priors->re_var_l[0]=psprior8;  /*prioor for subject variance */
	priors->re_var_l[1]=psprior9;
	priors->re_var_f[0]=psprior10;
	priors->re_var_f[1]=psprior11;


	fscanf(finput, "%lf %lf %lf %lf %lf %lf\n", &bprior1, &bprior2, &bprior3, &bprior4, &bprior5, &bprior6); /*6th baseline*/
    hyper->meanmeanbh_l[0] = bprior1;
    hyper->meanvarbh_l[0] = bprior2;
    hyper->varbh_l[0] = bprior3;
	hyper->meanmeanbh_f[0] = bprior4;
    hyper->meanvarbh_f[0] = bprior5;
    hyper->varbh_f[0] = bprior6;

    fscanf(finput,"%lf %lf %lf %lf %lf %lf\n", &hprior1, &hprior2, &hprior3,&hprior4, &hprior5, &hprior6); /*7th halflife*/
    hyper->meanmeanbh_l[1] = hprior1;
    hyper->meanvarbh_l[1] = hprior2;
    hyper->varbh_l[1] = hprior3;
	 hyper->meanmeanbh_f[1] = hprior4;
    hyper->meanvarbh_f[1] = hprior5;
    hyper->varbh_f[1] = hprior6;

    fscanf(finput,"%lf %lf %lf %lf\n", &sprior1, &sprior2, &sprior3, &sprior4); /*8th error*/
    priors->alpha_l = sprior1; 
    priors->beta_l= sprior2;
	 priors->alpha_f = sprior3;
    priors->beta_f = sprior4;
	
	/*fscanf(finput,"%d %d\n", &nprior1, &nprior2);*/ /*9th number of pulse*/
/*	parms_l->nprior = nprior1;
	parms_f->nprior = nprior2;*/
	priors->rho_prior[0] = -0.105;/*cluster size*/
	priors->rho_prior[1] = 1;
	priors->nu_prior[0] = 1.5;
	priors->nu_prior[1] = 10;
	parms_l->Rinter1 = 36;
	parms_l->Rinter2=96;



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

      mcmc(sublist,parms_l,parms_f,ts_l,ts_f,iter,*N,priors,seed,commonl,commonf,propvar,hyper);
  
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


