/*******************************************************************/
/************************death_rate_kernel.c ******************************/
/**To speed calcuation of the death rates for a cox process,
/*******************************************************************/

/***FIX NOW***/
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

 skernel

 **********************************************************************/


double *skernel(Node_type *list, Node_type *parents, Common_parms *parms_cox, int npulse) {

	int k;
	double *sumkernel, tempsum, temp;
	Node_type *node, *fnode;

	sumkernel = (double *)calloc(npulse, sizeof(double)); /*npulse is number of fsh pulse*/
	fnode = list->succ;
	k = 0;
	while (fnode != NULL){
		tempsum = 0;
		node = parents->succ;
		while (node != NULL) {
			temp = -1 / (2 * parms_cox->nu) * (fnode->time - node->time) * (fnode->time - node->time);
			if (temp > -685)
				tempsum += exp(parms_cox->lrho) / sqrt(2 * 3.14159 * parms_cox->nu) * exp(temp);
			node = node->succ;
		}
		sumkernel[k] = tempsum;
		k++;
		fnode = fnode->succ;
	}

	return sumkernel;
}



