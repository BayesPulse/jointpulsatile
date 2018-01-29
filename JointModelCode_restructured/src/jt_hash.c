
/*******************************************************************/
/****************************hash.h ********************************/
/******************************************************************************/
/* this module implements the linked list.  in this implementation, the first */
/* node in the linked list is a dummy node.  one must initialize the list with*/
/* the dummy node as the first.                                               */
/******************************************************************************/
/*******************************************************************/

#include "jt_deconvolution_main.h"

/*******************************************************************
*******************GLOBAL VARIABLE DEFINITIONS**********************

 fitstart: The first time in hours that a pulse may occur

*********************************************************************/

extern double fitstart;

/********************************************************************/
/*SUBROUTINES THAT EXIST IN THIS PROGRAM

 *initialize_node
 *sequential_search
 insert_node
 delete_node
 print_list
 destroy_list

**********************************************************************/

/*********************************************************************/
               /*START OF initialize_node SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*initialize_node: this allocates memory for a new pulse and gives it some
                   initial parameters
    ARGUMENTS: None
    RETURNS: p, the created pulse
**********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i: generic counter
 *p: initialized pulse

 SUBROUTINES USED
  None
************************************************************************/

Patient *initialize_subject(void)
{
  int i;
  Patient *s;
  PulseEstimates *initialize_node(void);

/* initialize memory for a subject */
  if ((s = (Patient *)malloc(sizeof(Patient))) == NULL) {
    printf("initialize_subject: out of memory, can't allocate node\n");
    exit(1);
  }
       
/* all nodes initialized with a pointer to the NULL pointer     */
/* and time = 0, the parameters all initialized to zero        */
/* it is up to the programmer to insure that the correct values */
/* are inserted into the node at the appropriate time           */

  s->succ = s->pred = NULL;
  

    s->patient_parms = calloc(1,sizeof(PatientEstimates));
    s->resp_patient_parms = calloc(1,sizeof(PatientEstimates));
    s->patient_parms->pulselist = initialize_node();
    s->resp_patient_parms->pulselist = initialize_node();
    
    s->patient_data = calloc(1,sizeof(PatientData));
    s->patient_pv = calloc(1,sizeof(PatientProposals));
    s->patient_pv_response = calloc(1,sizeof(PatientProposals));
    s->pulse_pv = calloc(1,sizeof(PulseProposals));
    s->resp_pulse_pv = calloc(1,sizeof(PulseProposals));
    
  return s;
}

PulseEstimates *initialize_node(void)
{
  int i;
  PulseEstimates *p;

/* initialize memory for a node */
    if ((p = (PulseEstimates *)malloc(sizeof(PulseEstimates))) == NULL) {
    printf("initialize_node: out of memory, can't allocate node\n");
    exit(1);
  }

/* all nodes initialized with a pointer to the NULL pointer     */
/* and time = 0, the parameters all initialized to zero        */
/* it is up to the programmer to insure that the correct values */
/* are inserted into the node at the appropriate time           */

    p->succ = p->pred = NULL;
    p->mean_contrib = NULL;
    p->time = fitstart;
    p->mass = 0.1;
    p->width = 0.1;
    p->tvarscalemass = 1;
    p->tvarscalewidth = 1;
//    p->lambda =0;
  return p;
}

/*********************************************************************/
               /*START OF sequential_search SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*sequential_search: finds where a pulse's location fits within an established
                     linked list
    ARGUMENTS: Pulse_Estimates *list; this is the current list of pulses that exist;
               double time; this is the location of a new pulse; note that we
                  do not actually input a new pulse, we just input its location
    RETURNS: *loc, the newly created pulse's prececessor
**********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 *loc: the newly created pulse's predecessor
 *locm1: the newly created pulse's successor

 SUBROUTINES USED
  None
************************************************************************/

PulseEstimates *sequential_search(PulseEstimates *list,double time)
{
  PulseEstimates *loc,*locm1;

/* searches through the linked list "list" for the first occassion of */
/* "list->time" that is less than "time".                             */
/* the new node will be inserted between "locm1" and "loc".  "locm1"  */
/* preceeds "loc" in "list" by one position                           */

  locm1 = list;
  for (loc=list;loc;loc=loc->succ){
    if (loc->time >= time)
      return locm1;
    else
      locm1 = loc;
  }
  return locm1;
}

/*********************************************************************/
               /*START OF insert_node SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*insert_node: integrates a newly created pulse into the linked list
    ARGUMENTS: Pulse_Estimates *new_node; the newly created pulse;
               Pulse_Estimates *list; this is the current list of pulses that exist;
    RETURNS: None; all updates are made internally
**********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 *node: new_node's predecessor

 SUBROUTINES USED
  sequential_search: found in this file; identifies inputted pulse's predecessor
***********************************************************************/

void insert_node(PulseEstimates *new_node,PulseEstimates *list)
{
  Pulse_Estimates *node;

/* insert the node "new_node" in the linked list "list" */
/* find the position in which to insert "new_node"      */

  node = sequential_search(list,new_node->time);
/********************************************************/

/* reassign pointers so that everyone is pointing to the correct node */

  new_node->pred = node;
  new_node->succ = node->succ;
  if (node->succ != NULL)
    (node->succ)->pred = new_node;
  node->succ = new_node;
}


void insert_subject(Patient *new_subj,Patient *sublist)
{

/* insert the node "new_node" in the linked list "list" */
/* find the position in which to insert "new_node"      */

/********************************************************/

/* reassign pointers so that everyone is pointing to the correct node */
  if (sublist->succ != NULL)
    (sublist->succ)->pred = new_subj;
  new_subj->pred = sublist;
  new_subj->succ = sublist->succ;
  sublist->succ = new_subj;

}

/*********************************************************************/
               /*START OF delete_node SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*delete_node: frees memory associated with inputted pulse and reassigns
               pointers of remaining pulses
    ARGUMENTS: Pulse_Estimates *node; the pulse to be deleted;
               Pulse_Estimates *list; this is the current list of pulses that exist;
    RETURNS: None; all updates are made internally
**********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
  None

 SUBROUTINES USED
  None
************************************************************************/

void delete_node(PulseEstimates *node,PulseEstimates *list)
{
  free(node->mean_contrib);
  
  if (node->succ == NULL)
    (node->pred)->succ = node->succ;
  else {
    (node->pred)->succ = node->succ;
    (node->succ)->pred = node->pred;
  }
  free(node);
}

/*********************************************************************/
               /*START OF print_list SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*print_list: prints information about all pulses in the linked list
    ARGUMENTS: Pulse_Estimates *list; this is the current list of pulses that exist;
    RETURNS: None; all updates are made internally
**********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i: generic counter
 *node: counter through pulses

 SUBROUTINES USED
  None
************************************************************************/

void print_list(PulseEstimates *list)
{
  int i;
  PulseEstimates *node;

/* traverses the linked list and prints the contents of each node */
  i = 1;
  node = list->succ;
  while (node != NULL) {
    printf("%2d %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf\n",i,node->time,node->mass,node->width,node->tvarscalemass,node->tvarscalewidth);
    node = node->succ;
    i++;
  }
}

/*********************************************************************/
               /*START OF destroy_list SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*destroy_list: frees memory throughout the linked list
    ARGUMENTS: Pulse_Estimates *list; this is the current list of pulses that exist;
    RETURNS: None; all updates are made internally
**********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 *loc: counter through pulses

 SUBROUTINES USED
  None
************************************************************************/

void destroy_list(PulseEstimates *list)
{
  PulseEstimates *loc;

  while (list->succ != NULL) {
    loc = list->succ;
    free(list->mean_contrib);
    free(list);
    list = loc;
  }
  free(list->mean_contrib);
  free(list);
}

void destroy_sublist(Patient *sublist)
{
  Patient *loc;
  void destroy_list(PulseEstimates *);

  while (sublist->succ != NULL) {
      loc = sublist->succ;
      destroy_list(sublist->patient_parms->pulselist);
      destroy_list(sublist->resp_patient_parms->pulselist);
      free(patient_parms);
      free(resp_patient_parms);
      free(patient_data->concentration);
      free(patient_data->time);
      free(patient_data->response_concentration);
      free(patient_data->common_filename);
      free(patient_data->pulse_filename);
      free(patient_data->resp_common_filename);
      free(patient_data->resp_pulse_filename);
      free(patient_data);
      free(patient_pv);
      free(patient_pv_response);
      free(pulse_pv);
      free(resp_pulse_pv);
      free(sublist);
    sublist = loc;
  }
    destroy_list(sublist->patient_parms->pulselist);
    destroy_list(sublist->resp_patient_parms->pulselist);
    free(patient_parms);
    free(resp_patient_parms);
    free(patient_data->concentration);
    free(patient_data->time);
    free(patient_data->response_concentration);
    free(patient_data->common_filename);
    free(patient_data->pulse_filename);
    free(patient_data->resp_common_filename);
    free(patient_data->resp_pulse_filename);
    free(patient_data);
    free(patient_pv);
    free(patient_pv_response);
    free(pulse_pv);
    free(resp_pulse_pv);

  free(sublist);
}
