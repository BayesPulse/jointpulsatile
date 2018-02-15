/*******************************************************************/
/************************death_rate_kernel_1.c ******************************/
/**This is a denoninator sum in the Cox death rate.  It is stored to speed computation time***/
/*******************************************************************/

/***FIX***/
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

skernel_1

 **********************************************************************/


double skernel_1(Node_type *fnode, Node_type *parents, Common_parms *parms_cox) {

	int k;
	double *sumkernel, tempsum, temp;
	Node_type *node;

	tempsum = 0;
	node = parents->succ;
	while (node != NULL) {
		temp = -1 / (2 * parms_cox->nu) * (fnode->time - node->time) * (fnode->time - node->time);
		if (temp > -685)
			tempsum += parms_cox->rho / sqrt(2 * 3.14159 * parms_cox->nu) * exp(temp);

		node = node->succ;
	}
	if (tempsum == 0)
		tempsum = 1e-300;

	return tempsum;
}






