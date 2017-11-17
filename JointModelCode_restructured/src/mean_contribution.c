/*******************************************************************/
/************************mean_contribution.c ******************************/
/*******************************************************************/

/**Fix here**/
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

 *mean_concentration

 **********************************************************************/


/*********************************************************************/
/*START OF mean_concentration SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*mean_concentration: this takes each pulse's mean_contrib vector and sums
					  across them
					  ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
					  Common_parms *parms; the current values of the common parameters;
					  int N; the number of observations in **ts;
					  Node_type *node_out; if we want, we can ignore a pulse--this is a critical use in bd algorithm;
					  double **ts; this is the matrix of observed data (a column of
					  times and a column of log(concentration);
					  RETURNS: x, the vector of sums
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i: generic counter
 *x: the vector of sums
 *node: the counter of pulses

 SUBROUTINES USED
 rnorm: found in randgen.h; draws from the normal distribution
***********************************************************************/

double *mean_concentration(Node_type *list, Common_parms *parms, int N, Node_type *node_out, double baseline)
{
	/* This function sums the mean_contrib vector across pulses */
	/* The result is a vector of sums */
	int i;
	double *x;
	Node_type *node;
	double rnorm(double, double, unsigned long *);

	x = (double *)calloc(N, sizeof(double));

	/* add the contribution to the mean from each pulse */
	node = list->succ;
	while (node != NULL) {
		if (node != node_out) {
			for (i = 0; i < N; i++)
				x[i] += node->mean_contrib[i];
		}
		node = node->succ;
	}

	/* add the baseline contribution */
	for (i = 0; i < N; i++) {
		x[i] += baseline;
		x[i] = log(x[i]);
	}

	return x;
}





