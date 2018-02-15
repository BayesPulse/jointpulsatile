//
//  mh_logsigmean.c
//  
//
//  Created by Nichole on 11/17/17.
//
//

#include "mh_logsigmean.h"

void mh_logsigmean(Subject_type *subject, Common_parms *parms_cox, Priors *priors_cox, unsigned long *seed, double psigmean_sd) {
    /*declare variables */
    int i, k, npulse2, npulses, nparents;
    double tmp;
    Node_type *tmpparent, *parent, *node, *fnode;
    Node_type *ploclist;
    Subject_type  * subject2;
    
    double NewMeanLogOmega, NewMeanOmega, lprop_ratio, *sumlogk, sumk, temp, diff_alpha, alpha, integral, tempsum, lprior_ratio, logexpterm, kold, knew, *denomsum, lratio, *pulselocs;
    
    /*declare functions */
    double rnorm(double, double, unsigned long *);
    double kiss(unsigned long *);
    double phi(double, double, double);
    double *skernel(Node_type *, Node_type *, Common_parms *, int);
    double skernel_1(Node_type *, Node_type *, Common_parms *);
    
    
    /*draw the new log(c_omega)   ***********/
    NewMeanLogOmega = rnorm(parms_cox->lnu, psigmean_sd, seed);
    
    NewMeanOmega = exp(NewMeanLogOmega);
    
    lratio = 0.;
    integral = 0;
    nnu++;
    knew = 0;
    kold = 0;
    tmp = parms_cox->nu;
    parms_cox->nu = NewMeanOmega;
    subject2 = subject->succ;  /****pull information from all the subject  ***/;
    
    while (subject2 != NULL){
        
        nparents = subject2->numnode_l;
        npulses = subject2->numnode_f;
        parent = subject2->driver;
        node = subject2->response;
        tmpparent = subject2->driver->succ;
        
        
        while (tmpparent != NULL){
            integral += exp(parms_cox->lrho)*(phi(fitend, tmpparent->time, sqrt(exp(NewMeanLogOmega))) - phi(fitstart, tmpparent->time, sqrt(exp(NewMeanLogOmega))) - phi(fitend, tmpparent->time, sqrt(exp(parms_cox->lnu))) + phi(fitstart, tmpparent->time, sqrt(exp(parms_cox->lnu))));
            tmpparent = tmpparent->succ;
        }
        fnode = node->succ;
        /*calculate the sum of old log lambda */
        while (fnode != NULL){
            kold += log(fnode->lambda);
            
            fnode = fnode->succ;
        }
        
        
        
        fnode = node->succ;
        while (fnode != NULL)
        {
            fnode->lambda = skernel_1(fnode, parent, parms_cox);
            knew += log(fnode->lambda);
            
            fnode = fnode->succ;
        }
        /*calculate the sum of old log lambda */
        
        
        
        
        
        
        
        subject2 = subject2->succ;
    }
    /*printf("endloop\n");*/
    /*calculate the acceptance ratio */
    lratio = knew - kold - integral;
    lprior_ratio = 1 / (2 * priors_cox->nu_prior[1]) * ((parms_cox->lnu - priors_cox->nu_prior[0]) * (parms_cox->lnu - priors_cox->nu_prior[0]) - (NewMeanLogOmega - priors_cox->nu_prior[0]) * (NewMeanLogOmega - priors_cox->nu_prior[0]));
    lprop_ratio = lratio + lprior_ratio;
    /*lprop_ratio = knew-kold;*/
    
    /*   printf("lprop_ratio %lf\n",lprop_ratio);*/
    /*make the decision on whether to keep the new value */
    alpha = (0 < lprop_ratio) ? 0 : lprop_ratio;
    
    if (log(kiss(seed)) < alpha) {
        /*         printf("accept\n");*/
        anu++;
        parms_cox->lnu = NewMeanLogOmega;
    }
    else{
        parms_cox->nu = tmp;
        subject2 = subject->succ;
        while (subject2 != NULL)
        {
            fnode = subject2->response->succ;
            parent = subject2->driver;
            
            while (fnode != NULL)
            {
                fnode->lambda = skernel_1(fnode, parent, parms_cox);
                
                fnode = fnode->succ;
            }
            subject2 = subject2->succ;
        }
    }
    
    
}

