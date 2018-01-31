#pragma once
#include "jt_deconvolution_main.h"
Patient *initialize_subject(void);
PulseEstimates *initialize_node(void);
PulseEstimates *sequential_search(PulseEstimates *list, double time);
void insert_node(PulseEstimates *new_node, PulseEstimates *list);
void insert_subject(Patient *new_subj, Patient *sublist);
void delete_node(PulseEstimates *node, PulseEstimates *list);
void print_list(PulseEstimates *list);
void destroy_list(PulseEstimates *list);
void destroy_sublist(Patient *sublist);
