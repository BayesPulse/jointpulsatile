/*******************************************************************/
/************************birthdeath.c ******************************/
/*******************************************************************/

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

bd_response
 **********************************************************************/


/*********************************************************************/
/*bd_response: this runs the birth-death step of the BDMCMC process
	ARGUMENTS: 
    Subject_type *sublist: This is the linklist of the subjects information.
   double **ts: Matrix of observed data from all subjects.
   Common_parms *parms: The current values of the population parameters.
   Common_parms *parms_cox: The current values in the cox process, which are common across all subjects.
    int N; the number of observations in **ts;
	unsigned long *seed; seed values needed for the randon number generator;
	int iter; which iteration are we on;
 
 Subject_type *sublist, double **ts, Common_parms *parms, Common_parms *parms_cox, int N,
	unsigned long *seed, int iter
	RETURNS: None; all updates are made internally
 *********************************************************************/


/*********************************************************************/
/***NOT UPDATED!!!***/
/*VARIABLE DEFINITIONS
 i,j,k: generic counters
 remove: if a death occurs, remove is a draw from the multinomial distribution
 and represents which pulse we will remove
 num_node: number of pulses
 flag: if a birth occurs, we must draw parameters for this new pulse; flag will
 be set to 0 until a positive mass and a positive width are drawn
 aaa: this counter ensures we do not run the birthdeath loop too many times
 max_num_node=60: this variable is set to 60 and ensures a death if we exceed
 60 pulses
 S: the sum of exponential draws; the birthdeath loop stops when S exceeds T
 Birth_rate: a constant birth rate needed for the process
 T=1: stop birthdeath loop when S exceeds T=1
 full_likelihood: current likelihood before a birth or death occurs
 max: part of the calculation for Death_rate
 Death_rate: total death rate
 *death_rate: vector of death rates for each pulse
 position: if a birth occurs, the pulse's location is position, drawn from a
 uniform distribution
 *partial_likelihood: a vector of likelihoods where the ith element is the
 likelihood if pulse i were removed
 *tmp: used when drawing new pulse's mass and width
 *node: counter through pulses
 *new_node: if a birth occurs, we give the new pulse this name

 SUBROUTINES USED
 rnorm: found in randgen.h; draws from the normal distribution
 rmultinomial: found in randgen.h; draws from the multinomial distribution
 kiss: found in randgen.h; draws from U(0,1) distribution
 rexp: found in randgen.h; draws from the exponential distribution
 runif_atob: found in randgen.h; draws from U(a,b) distribution
 mean_contribution: found in this file; evaluates mean_contrib for a
 pulse with specific parameters
 likelihood: found in this file; calculates and returns likelihood under
 current parameters Node_type *initialize_node(void);
 *initialize_node: found in hash.h; allocates memory for a new pulse and creates it;
 insert_node: found in hash.h; inserts a newly created node in the appropriate
 spot in the linked list
 delete_node: found in hash.h; eliminates a node and links the neighbors
 *calc_death_rate: found in this file; creates a vector where the ith element
 represents the death rate for pulse i
 ***********************************************************************/

void birth_death_response(Subject_type *sublist, double **ts, Common_parms *parms, Common_parms *parms_cox, int N,
	unsigned long *seed, int iter)
{
	int i, j, k, remove, num_node, flag, aaa, ninter1,max_num_node = 60;
	double S, Birth_rate, tot_birth_rate, T = 1., full_likelihood, full_likelihood2, max;
	double Death_rate, *death_rate, position, *partial_likelihood, *tmp,tmp2[2];
	Node_type *node, *new_node, *list;
	Subject_type *subject;
	double rnorm(double, double, unsigned long *); void print_list(Node_type *);
	long rmultinomial(double *, long, unsigned long *);
	double kiss(unsigned long *);
	double rgamma(double, double, unsigned long*);

	double rexp(double, unsigned long *);
	double runif_atob(unsigned long *, double, double);
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);
	double likelihood(Subject_type *list, double **ts, Common_parms *parms, int N,
		Node_type *node_out);
	double likelihood2(Node_type *list, double **ts, Common_parms *parms, int N,
		Node_type *node_out, double baseline);
	Node_type *initialize_node(void);
	void insert_node(Node_type *new_node, Node_type *list);
	/* void delete_node(Node_type *node, Node_type *list);*/
	double *calc_death_rate_f(Node_type *, Node_type *, int, double *, double, double, Common_parms *);
	double skernel_1(Node_type *, Node_type *, Common_parms *);


	/*For first 100 iterations, birth rate is 10, then it becomes 1;
	This is to force more births early on.*/
	if (iter < 100) {
		tot_birth_rate = 10.;
	}
	else {
		tot_birth_rate = 24 * 4.;
	}
	Birth_rate = tot_birth_rate / 1520;

	subject = sublist->succ;

	parms->subindex = 0;


	while (subject != NULL){

		list = subject->response;
		aaa = 0;

		S = 0.0;
		tmp = (double *)calloc(2, sizeof(double));

		/*Save Likelihood*/
		/*full_likelihood = *likeli;
		*/
		/*Go until the loop is broken*/

		while (1) {

			/*This counter keeps this loop from running too many times*/
			aaa++;

			/*Count number of pulses*/
			num_node = 0;
			node = list->succ;
			while (node != NULL) {
				num_node++;
				node = node->succ;
			}


			/*Allocate memory for partial likelihood vector*/
			partial_likelihood = (double *)calloc(num_node, sizeof(double));

			/*Calculate the likelihood if pulse i is removed*/
			i = 0;
			node = list->succ;

			/*  printf("num_node = %d \n",num_node);i*/

			while (node != NULL) {
				partial_likelihood[i] = likelihood2(subject->response, ts, parms, N, node, subject->basehalf_f[0]);
				/*printf("partial likelihood = %lf \n",partial_likelihood[i]);*/
				i++;
				node = node->succ;
			}
			subject->numnode_f = num_node;
			full_likelihood2 = likelihood2(subject->response, ts, parms, N, subject->response, subject->basehalf_f[0]);

			/* CALCULATE DEATH RATE FOR EACH COMPONENT */
			death_rate = NULL;
			if (iter < 1){
				death_rate = calc_death_rate(list, num_node, partial_likelihood, full_likelihood2, tot_birth_rate, 10.0);
			}
			else{
				death_rate = calc_death_rate_f(subject->driver, list, num_node, partial_likelihood, full_likelihood2, Birth_rate, parms_cox);
			}
			/*    node = list->succ;
			i=0;
			while (node != NULL) {
			printf("death_rate = %lf \n",death_rate[i]);
			i++;
			node = node->succ;
			}
			*/
			/*The next portion computes D = sum(d_i); This is a little more complicated
			than summing them up because of precision issues.*/
			if (death_rate != NULL) {
				Death_rate = death_rate[0];
				for (i = 1; i<num_node; i++) {
					max = (Death_rate > death_rate[i]) ? Death_rate : death_rate[i];
					Death_rate = log(exp(Death_rate - max) + exp(death_rate[i] - max)) + max;
				}

				for (i = 0; i < num_node; i++) {
					death_rate[i] -= Death_rate;
					death_rate[i] = exp(death_rate[i]);
				}

				for (i = 1; i < num_node; i++)
					death_rate[i] += death_rate[i - 1];

				Death_rate = exp(Death_rate);
			}
			else
				Death_rate = 0;


			free(partial_likelihood);


			if (num_node <= 1) Death_rate = 0;

			/* Draw from exp(B+D) and add to current S */
			S += rexp(tot_birth_rate + Death_rate, seed);

			/*If S exceeds T or if we've run this too many times, break*/
			if (S > T)  break;
			if (aaa > 5000) break;

			/* SIMULATE JUMP TYPE (BIRTH OR DEATH) */
			/* If we only have 0 or 1 pulses, we definitely have a birth*/
			if (num_node <= 1)
				max = 1.1;

			/* If we have too many pulses (60), we definitely have a death*/
			else if (num_node >= max_num_node)
				max = -0.1;

			/* Otherwise, set max = B/(B+D) */
			else
				max = tot_birth_rate / (tot_birth_rate + Death_rate);

			/*    printf("Birth_rate = %lf \n",Birth_rate);
			printf("Death_rate = %lf \n",Death_rate);
			*/
			/*If U < B/(B+D), a birth occurs */
			if (kiss(seed) < max) {
				ninter1 = 0;
				flag=1;
				/*simulate a pulse position uniformly */
				/*if the new pulse location is too close to the old one, try it again*/
				/*like an hard core process*/
				while (flag==1)
				{
					
					position = runif_atob(seed, fitstart, fitend);
					node = list->succ;
					ninter1 =0;
					while (node != NULL)
					{

						if (fabs(node->time - position) <10) {
							ninter1++;
						}
						node = node->succ;
					}
					if (ninter1 == 0) {
						flag = 0;
					
						
					}

				}
				subject->numnode_f++;
				node = list->succ;
				flag = 1;

				/* simulate random effects from the prior */
				if (flag) {
					for (k = 0; k<100; k++) {
						for (j = 0; j<2; j++){
							tmp2[j] = rgamma(2., 2., seed);
						tmp[j] = rnorm(subject->theta_f[j], parms->re_precision[j]/sqrt(tmp2[j]), seed);
					}
						if ((tmp[0] > 0) && (tmp[1] > 0)) {
							flag = 0;
							break;
						}

					}
					/*Only continue once we have positive values of mass and width*/
					if (!flag) {
						/* initialize a new node */
						new_node = initialize_node();
						new_node->time = position;
						new_node->theta[0] = tmp[0];
						new_node->theta[1] = tmp[1];
						new_node->eta[0] = tmp2[0];
						new_node->eta[1] = tmp2[1];


						/*          new_node->width = exp(new_node->theta[1]);*/
						new_node->mean_contrib = (double *)calloc(N, sizeof(double));
						mean_contribution(new_node, ts, parms, N, subject->basehalf_f[1]);
						new_node->lambda = skernel_1(new_node, subject->driver, parms_cox);
						/* insert the new node; this routine also mskes sure it's
						located properly */
						insert_node(new_node, list);
					}
				}
			}

			/*Otherwise, a death occurs */
			else {
				/* Pick a node to remove based on the multinomial */

				remove = (int)rmultinomial(death_rate, (long)num_node, seed) + 1;
				node = list;
				for (i = 0; i < remove; i++)
					node = node->succ;
				delete_node(node, list);
				subject->numnode_f--;

			}
			free(death_rate);
			/*  fflush(stdout);*/
		} /*End of While Loop*/

		/*Update likelihood before leaving Birthdeath.c*/
		/**likeli = full_likelihood;
		*/

		free(tmp);
		free(death_rate);

		subject = subject->succ;
		parms->subindex++;

	}

}





