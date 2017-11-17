/*******************************************************************/
/************************birthdeath.c ******************************/
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

 mean_contribution
 **********************************************************************/


/*********************************************************************/
/*START OF mean_contribution SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*mean_contribution: this updates a pulse's mean_contrib vector based on current
					 values of parameters
					 ARGUMENTS: Node_type *node; what pulse are we updating;
					 double **ts; this is the matrix of observed data (a column of
					 times and a column of log(concentration);
					 Common_parms *parms; the current values of the common parameters;
					 int N; the number of observations in **ts;
					 RETURNS: None; all updates are made internally
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i: generic counter
 x,y,z,tmp: these are all variables that are part of the arithmetic needed in
 calculating the mean contribution

 SUBROUTINES USED
 erf: this is a subroutine used to help in integration
 ***********************************************************************/

void mean_contribution(Node_type *node, double **ts, Common_parms *parms, int N, double halflife)
{
	/* This function updates a pulse's mean_contrib vector based on inputted parms*/
	int i;
	double x, y, z, tmp;
	double erf(double);

	z = node->theta[1] * 0.6931472 / halflife;
	y = 0.6931472*(0.5*z / halflife + node->time / halflife);
	z += node->time;
	tmp = sqrt(2.*node->theta[1]);
	for (i = 0; i < N; i++) {
		x = (ts[i][0] - z) / tmp;
		x = 0.5*(1. + erf(x));
		if (x == 0) node->mean_contrib[i] = 0;
		else node->mean_contrib[i] = node->theta[0] * x*exp(y - ts[i][0] * 0.6931472 / halflife);

	}

}




