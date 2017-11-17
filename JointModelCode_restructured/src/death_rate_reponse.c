/*******************************************************************/
/************************death_rate_response.c ******************************/
/*******************************************************************/

/***FIX HERE***/
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

 *calc_death_rate_response

 **********************************************************************/


/*********************************************************************/
/*START OF calc_death_rate_response SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*calc_death_rate: calculates a vector of death rates, one for each existing pulse
	ARGUMENTS: Node_tye *parent: this is a link list of the trigger pulse location (the cluster centers in the Cox process) for a subject.
     Node_type *list; this is the current list of pulses that exist for a subject;
	int num_node; current number of pulses
	double *partial_likelihood; vector of likelihoods, where the ith
	element represents the likelihood with the ith pulse removed for this subject
	double full_likelihood; value of the full likelihood for this subject
	double Birth_rate; value of the birth rate

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
double *calc_death_rate_response(Node_type *parent, Node_type *list, int num_node, double *partial_likelihood, double full_likelihood, double Birth_rate, Common_parms *parms_cox)
{
	/* This function calculates the death rate vector */
	int i, j;
	double x, *death_rate, intsum;
	double coef_denom, coef_num;
	Node_type *node, *pnode;
	double skernel_1(Node_type *, Node_type *, Common_parms *);
	if (num_node > 1) {
		death_rate = (double *)calloc(num_node, sizeof(double));
		node = list->succ;
		i = 0;

		/* Calculate the coefficient of the distribution of the taus conditional
		 on number of pulses. In this portion, I have an extra num_node in the
		 numerator */

		while (node != NULL) {
			/*		pnode = parent->succ; */
					/*calculate the conditional intensity of pulse l*/
			/*intsum = 0;
			death_rate[i] = 0.01;
			while (pnode != NULL) {
			intsum += exp(parms_cox->lrho) / sqrt(2 * 3.14159*parms_cox->nu) * exp(-1 / (2 * parms_cox->nu) * (node->time - pnode->time) * (node->time - pnode->time));
			pnode = pnode->succ;
			}  */
			if (node->lambda != 0)
			{

				x = partial_likelihood[i] - full_likelihood + log(Birth_rate) - log(node->lambda);
				death_rate[i] = x;
			}
			else{
				x = partial_likelihood[i] - full_likelihood + log(Birth_rate) - log(1e-300);
				death_rate[i] = x;
			}
			node = node->succ;
			i++;

		}


		/*Save to death rate vector */




		/*Advance to next pulse*/

		/*Advance counter 1*/

	}
	else {

		/*if we have 0 or 1 pulses*/
		if (num_node == 0) /*if we don't have any pulses, no death rate*/
			return NULL;

		else {

			/*If we have 1 pulse, don't kill it*/
			death_rate = (double *)calloc(num_node, sizeof(double));
			death_rate[0] = -1e300;
		}
	}
	return death_rate;
}




