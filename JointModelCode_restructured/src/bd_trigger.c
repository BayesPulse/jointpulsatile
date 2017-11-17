/*******************************************************************/
/************************bd_driver.c ******************************/
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

 birth_death_driver

 **********************************************************************/

/*********************************************************************/
/*START OF birth_death SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*birth_death_driver: this runs the birth-death step of the BDMCMC process
	ARGUMENTS: 
 
    Subject_type *sublist, is a link list of the different subject in the analysis.
    double **ts, data matrix.  One column for each subject's data.  Column of observations.
	Common_parms *parms the current values of the population parameters, which are common across subjects;
	int N the number of observations in **ts;
	unsigned long *seed: seed values needed for the randon number generator;
	int iter: which iteration are we on;
	RETURNS: None all updates are made internally
 *********************************************************************/
/*********************************************************************/

/*NOT UPDATED FOR BD_DRIVER SUBROUTINE*/
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

void birth_death_driver(Subject_type *sublist, double **ts, Common_parms *parms, int N,
	unsigned long *seed, int iter)
{
	int i, j, k, remove, num_node, flag, aaa, max_num_node = 24, ninter1, ninter2;
	double S, Birth_rate, T = 1., full_likelihood, full_likelihood2, max, lambda, knew, r;
	double Death_rate, *death_rate, position, *partial_likelihood, *tmp, br_parent, tmp2[2];
	Node_type *node, *new_node, *list, *fnode, *cnode;
	Subject_type *subject;
	double rnorm(double, double, unsigned long *);
	double rgamma(double, double, unsigned long  *);
	void print_list(Node_type *);
	long rmultinomial(double *, long, unsigned long *);
	double kiss(unsigned long *);
	double rexp(double, unsigned long *);
	double runif_atob(unsigned long *, double, double);
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);

	double likelihood2(Node_type *list, double **ts, Common_parms *parms, int N,
		Node_type *node_out, double baseline);
	Node_type *initialize_node(void);
	void insert_node(Node_type *new_node, Node_type *list);
	double *calc_death_rate_l(Node_type *, Node_type *, int, double *, double, double, Common_parms *);
	double *calc_death_rate(Node_type *, int, double *, double, double, double);


/**THIS IS PROBLEM SPECIFIC CODE AND NEEDS TO BE CHANGED***/

	/*br_parent is the birth rate in a day*/
	br_parent = parms->beta;
	if (iter < 10) {
		Birth_rate = 15 * br_parent;
	}
	else {
		Birth_rate = 24 * br_parent * 2;
	}
/***CHANGE WILL END HERE***/

	subject = sublist->succ;

	parms->subindex = 0;

	while (subject != NULL){

		list = subject->driver;
		aaa = 0;  /**THIS IS BREAK CODE AND SHOULD BE CHANGED TO BE A CHECK AND LOOP EXIT WITHOUT A BREAK COMMAND--BREAK COMMAND IS GOING TO BE AN ISSUE WITH R INTEGRATION I BET**/

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

			subject->numnode_l = num_node;
			/*Allocate memory for partial likelihood vector*/
			partial_likelihood = (double *)calloc(num_node, sizeof(double));

			/*Calculate the likelihood if pulse i is removed*/
			i = 0;
			node = list->succ;

			/*  printf("num_node = %d \n",num_node);*/

			while (node != NULL) {
				partial_likelihood[i] = likelihood2(subject->driver, ts, parms, N, node, subject->basehalf_l[0]);
				/*printf("partial likelihood = %lf \n",partial_likelihood[i]);*/
				i++;
				node = node->succ;
			}

			full_likelihood2 = likelihood2(subject->driver, ts, parms, N, subject->driver, subject->basehalf_l[0]);

			/* CALCULATE DEATH RATE FOR EACH COMPONENT */
			death_rate = NULL;
			r = 10.0;
			if (iter < 1){
				death_rate = calc_death_rate(list, num_node, partial_likelihood, full_likelihood2, Birth_rate, r);
			}

			else {

				death_rate = calc_death_rate_l(subject->driver, subject->response, subject->numnode_l, partial_likelihood, full_likelihood2, Birth_rate, parms);
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
			S += rexp(Birth_rate + Death_rate, seed);

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
				max = Birth_rate / (Birth_rate + Death_rate);

			/*    printf("Birth_rate = %lf \n",Birth_rate);
			   printf("Death_rate = %lf \n",Death_rate);
			   */
			/*If U < B/(B+D), a birth occurs */
			if (kiss(seed) < max) {

				/*simulate a pulse position uniformly */
				position = runif_atob(seed, fitstart, fitend);
				/***part 1: calculate interaction sums ***/
				ninter1 = 0;
				ninter2 = 0;
				node = list->succ;

				cnode = list->succ;
				while (cnode != NULL && ninter1 == 0) {
					if (fabs(cnode->time - position) <= parms->Rinter1) {
						ninter1++;
					}
					else {
						if (fabs(cnode->time - position) <= parms->Rinter2) {
							ninter2++;
						}
					}
					cnode = cnode->succ;
				}

				if (ninter1 == 0) {
					lambda = parms->beta * exp((double)ninter2 * log(parms->gamma));
				}
				else
					lambda = 0;

				if (kiss(seed) <= lambda / br_parent) {

					flag = 1;
					subject->numnode_l++;
					/* simulate random effects from the prior */
					if (flag) {
						for (k = 0; k<100; k++) {
							for (j = 0; j<2; j++){
								tmp2[j] = rgamma(2.0, 2.0, seed);

								tmp[j] = rnorm(subject->theta_l[j], parms->re_precision[j]/sqrt(tmp2[j]), seed);
							}
							if ((tmp[0] > 0) && (tmp[1] > 0)) {
								flag = 0;
								break;
							}

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
						mean_contribution(new_node, ts, parms, N, subject->basehalf_l[1]);
						/* insert the new node; this routine also mskes sure it's
						 located properly */
						insert_node(new_node, list);
						fnode = subject->response->succ;
						while (fnode != NULL){
							knew = -1 / (2 * parms->nu)* (fnode->time - position) * (fnode->time - position);

							if (fnode->lambda == 0) {
								fnode->lambda = exp(parms->lrho) / sqrt(2 * 3.14159 * parms->nu) * exp(knew);
							}
							else {
								if (fabs(log(fnode->lambda) - knew) >= 22) {     /*adding the new location is different from the old*/
									if (knew > log(fnode->lambda))  /*new location adds more to the weight */
										fnode->lambda = exp(parms->lrho) / sqrt(2 * 3.14159 * parms->nu) * exp(knew);  /*just use the information of the new location  */
								}
								else {
									fnode->lambda += exp(parms->lrho) / sqrt(2 * 3.14159 * parms->nu) * exp(knew);
								}
							}
							fnode = fnode->succ;
						} /*end of caculating lambda*/

					} /*end of flag=0*/
				}/*end of pulse location*/

			} /*end of birth occurs*/
			/*Otherwise, a death occurs */
			else {
				/* Pick a node to remove based on the multinomial */
				//subject->numnode_l--;
				remove = (int)rmultinomial(death_rate, (long)num_node, seed) + 1;
				node = list;
				for (i = 0; i < remove; i++)
					node = node->succ;
				fnode = subject->response->succ;
				while (fnode != NULL){
					knew = -1 / (2 * parms->nu)* (fnode->time - node->time) * (fnode->time - node->time);
					if (log(fnode->lambda) - knew >= 0) {
						/* printf("k %d ldenom %lf tempterm %lf \n",k,log(denomsum[k]),tempterm);*/
						fnode->lambda -= parms->rho / sqrt(2 * 3.14159 *parms->nu) * exp(knew);
					}
					else
						fnode->lambda = 0;
					fnode = fnode->succ;
					/*printf("k %d denomsum[k] %le\n",k,denomsum[k]);*/
				}
				delete_node(node, list);


			}
			free(death_rate);
			/*  fflush(stdout);*/
		} /*End of While Loop*/

		/*Update likelihood before leaving Birthdeath.c*/
		/**likeli = full_likelihood;
	  */
		/* if (num_node>1)
			  free(death_rate);*/

		free(tmp);
		free(death_rate);
		/*		free(partial_likelihood);*/

		subject = subject->succ;
		parms->subindex++;

	}
}

