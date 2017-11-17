/*******************************************************************/
/************************death_rate_trigger.c ******************************/
/*******************************************************************/

/**FIX HERE**/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"

/*******************************************************************
*******************GLOBAL VARIABLE DEFINITIONS**********************

fitstart: The first time in hours that a pulse may occur
fitend: The last time in hours that a pulse may occur
mmm: Order statistic used for distribution of pulse locations.
This is inputted by the user and is typically 3

*********************************************************************/

extern double fitstart;
extern double fitend;
extern int mmm;

/********************************************************************/
/*SUBROUTINES THAT EXIST IN THIS PROGRAM

 *calc_death_rate_trigger

 **********************************************************************/



/*********************************************************************/
/*START OF calc_death_rate_trigger SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*calc_death_rate_trigger: calculates a vector of death rates, one for each existing pulse
	ARGUMENTS: 
       Node_type *parent: link list of the cluster centers and cluster specific parameters
        Node_type *list: link list of pulse locations for a subject
      int num_node: current number of pulses on this subject
    , double *partial_likelihood, double full_likelihood, double Birth_rate, Common_parms *parms_cox
 Node_type *list; this is the current list of pulses that exist;
	int num_node; current number of pulses
	double *partial_likelihood; vector of likelihoods, where the ith
	element represents the likelihood with the ith pulse removed
	double full_likelihood; value of the full likelihood
	double Birth_rate; value of the birth rate
    par
	RETURNS: death_rate; a vector where the ith element represents the death
	rate of the ith pulse
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j: generic counters
 x: variable for death rate
 *death_rate: vector of death rates
 coef_denom, coef_num: part of the calculation of death rate
 *node: counter for pulses

 SUBROUTINES USED
 None
***********************************************************************/

/* calculate the log of death rate*/
double *calc_death_rate_trigger(Node_type *parent, Node_type *list, int num_node, double *partial_likelihood, double full_likelihood, double Birth_rate, Common_parms *parms_cox)
{
	/* This function calculates the death rate vector */
	int i, j, ninter1, ninter2;
	double x, *death_rate, itsum, ratio, kdeath, lkdeath, denom;
	double coef_denom, coef_num;
	Node_type *node, *tmpcnode, *deathnode;
	double phi(double, double, double);

	/*calculate the actually value when number of lh pulse is greater than 1*/
	if (num_node > 1) {
		death_rate = (double *)calloc(num_node, sizeof(double));
		i = 0;
		deathnode = parent->succ;

		while (deathnode != NULL){
			/* Calculate the coefficient of the distribution of the taus conditional
			 on number of pulses. In this portion, I have an extra num_node in the
			 numerator */
			node = list->succ;
			ratio = 0;

			while (node != NULL){
				lkdeath = -1 / (2 * parms_cox->nu) * (node->time - deathnode->time) * (node->time - deathnode->time);
				kdeath = exp(parms_cox->lrho) / sqrt(2 * 3.14159 *parms_cox->nu) * exp(lkdeath);
				/*printf("k %d denomsum[k] %le kdeath %le diff %le ratio %le lr %le\n",k,denomsum[k],kdeath,(denomsum[k] - kdeath),denomsum[k]/(denomsum[k] - kdeath),log(denomsum[k]/(denomsum[k] - kdeath)));*/
				if (node->lambda > kdeath) {
					/* printf("loop \n");*/
					/*for model1: no epsilon*/
					/*ratio += log(parms_cox->epsilon + denomsum[k]) - log(parms_cox->epsilon + denomsum[k] - kdeath); */;
					ratio += log(node->lambda) - log(node->lambda - kdeath);
					/* printf("k %d ratio %le ldenomsum %lf lkdeath %lf  diff %le  \n",k,ratio,log(denomsum[k]),lkdeath,denomsum[k] - kdeath);*/
				}
				else if (node->lambda != 0){
					tmpcnode = parent->succ;
					denom = 0;
					while (tmpcnode != NULL) {
						if (tmpcnode->time != deathnode->time) {
							denom += exp(parms_cox->lrho) / sqrt(2 * 3.14159 *parms_cox->nu) * exp(-1 / (2 * parms_cox->nu) * (node->time - tmpcnode->time) * (node->time - tmpcnode->time));
						}
						tmpcnode = tmpcnode->succ;
					}
					if (denom != 0){

						ratio += log(1 + kdeath / (denom));
					}
					else
					{
						ratio += 10;
					}

				}
				else ratio += 0;
				node = node->succ;

				/* printf("k %d ratio %lf\n",k,ratio);   */
			}


			death_rate[i] = exp(parms_cox->lrho)*(phi(fitend, deathnode->time, sqrt(parms_cox->nu)) - phi(fitstart, deathnode->time, sqrt(parms_cox->nu))) - ratio + partial_likelihood[i] - full_likelihood;
			deathnode = deathnode->succ;
			i++;


		}
	}
	else{
		if (num_node == 0)
			return NULL;
		else {
			death_rate = (double *)calloc(1, sizeof(double));
			death_rate[0] = -1e300;
		}
	}

	return death_rate;
}

