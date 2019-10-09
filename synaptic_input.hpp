//
//  synaptic_input.hpp
//
//
#ifndef _synaptic_input_hpp
#define _synaptic_input_hpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <algorithm>

using namespace std;

#define NCOL      3
#define NROW      3

#define DURATION  1        /* of a binary spike, in msec */
#define TAU       2      /* in alpha function */
#define TMAX      26       /* duration for a response (precision 0.1%) */

/* allocates the binary array (0 or 1) and the double one */

const double tstep = 0.5;

void allocation (int size, int **pbin, double **ptraited) {
    int i;
    
    if (size <=0 ) {
        printf ("error : length should be positive\n");
        return;
    }
    
    *pbin = (int *) malloc (size * sizeof (int));
    *ptraited = (double *) malloc (size * sizeof (double));
    
    if ((*pbin == NULL) || (*ptraited == NULL)) {
        printf ("error: cannot allocate memory\n");
        return;
    }
    for (i=0; i<size; i++) {
        (*pbin)[i] = 0;
        (*ptraited)[i] = 0.;
    }
}

void allocate_train (int size, double **train) {
    int i;
    
    (*train) = (double *) malloc (size * sizeof (double));
    if ((*train) == NULL) {
        printf ("error : cannot allocate memory\n");
        return;
    }
    for (i=0; i<size; i++)
        (*train)[i] = 0.;
}

void sum (int size, double *a, double *b) {
    int i;
    
    for (i=0; i<size; i++)
        a[i] = a[i] + b[i];
}

void multiply (int size, double *a, double lambda) {
    int i;
    
    for (i=0; i<size; i++)
        a[i] = lambda * a[i];
}

/* generates the spikes according to a Poisson process of frequence 'freq' */

void gen_bin (int size, int *bin, int freq, double tstep) {
    int i, j;
    double duration = DURATION;
    int npoint;
    double p;
    
    npoint = (int) (duration / tstep);
    p = exp (-(double)freq/1000.0);
    
    for (i=0; i<size; i++) {
        if (drand48() > p)
            for (j=0; j<npoint; j++)
                bin [i+j] = 1;
        else
            for (j=0; j<npoint ;j++)
                bin [i+j] = 0;
    }
}


/* alpha function (modelizes a synapse) */

double alpha (double t_local, double tau) {
    double result;
    double Tau;
    
    Tau=tau;
    result = t_local / Tau * exp (1-t_local/Tau);
    return result;
}

void gen_traited (int size, int *bin, double *traited, double tstep) {
    int i, j;
    int duration = DURATION;
    double tau = TAU;
    double tmax = TMAX;
    
    for (i=0; i<size; i += (int) (duration/tstep)) {
        if (bin[i] == 1) {
            for (j=i; (j<size) && ((j-i) < (int)(tmax/tstep)); j++)
                traited [j] = alpha ((double)((j-i)*tstep), tau);
        }
    }
}


void generate_synaptic_input(vector<double> &synaptic, double freq, int nsyn_specific, int size) {
    int i;
    int c = 0;
    int row;
    int col;
    int *bin;
    double *traited;
    double *train_spe;
    
    allocation (size, &bin, &traited);
    allocate_train(size, &train_spe);
    
    if(nsyn_specific>0){
        for (i=0;i<nsyn_specific;i++){
            gen_bin(size, bin, freq, tstep);
            gen_traited(size, bin, traited, tstep);
            sum(size, train_spe, traited);
        }
        multiply(size, train_spe, 1./(nsyn_specific));
    }
    else if(nsyn_specific==0){
        for (i=0;i<size;i++) *(train_spe+i)=0;
    }
    
    synaptic.clear();
    synaptic.insert(synaptic.begin() , train_spe, train_spe + size);
    
    free(train_spe);
    free (bin);
    free(traited);
}

#endif
