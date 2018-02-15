/*******************************************************************/
/************************subject_likelihood.c ******************************/
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


 subj_likelihood

/*********************************************************************/
/*START OF likelihood SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*subj_likelihood: computes the current likelihood for the subject because subjects are independent
                   using the observed log-concentrations and mean concentration
			  ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
			  double **ts; this is the matrix of observed data (a column of
			  times and a column of log(concentration);
			  Common_parms *parms; the current values of the common parameters;
			  int N; the number of observations in **ts;
			  Node_type *node_out; if we want, we can ignore a pulse, this is critical for bd algorithm;

			  RETURNS: x, the likelihood computed based on the inputs
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i: generic counter
 x=0: likelihood to be calculated; initialized at zero
 *mean: vector of sums of mean_contribs calculated using mean_concentration

 SUBROUTINES USED
 mean_concentration: found in this file; computes a vector of sums of mean_contribs
 ***********************************************************************/


double subj_likelihood(Node_type *list, double **ts, Common_parms *parms, int N,
	Node_type *node_out, double baseline)
{
	/* This function computes the likelihood under inputted parameters */
	/* The output is a scalar */
	int i, j;
	double x = 0, *mean;

	double *mean_concentration(Node_type *, Common_parms *, int, Node_type *, double);

	j = parms->subindex;

	/* Sum across mean_contribs */
	mean = mean_concentration(list, parms, N, node_out, baseline);
	for (i = 0; i < N; i++)
		x += (ts[i][j + 1] - mean[i])*(ts[i][j + 1] - mean[i]);

	x /= (-2.0*parms->sigma);
	x += -0.5*N*(1.8378771 + parms->lsigma);
	free(mean);
	return x;
}





