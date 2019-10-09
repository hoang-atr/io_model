/*
 * io_spike.cpp
 *
 * This program generates spike trains of 3x3 IO neurons
 * Using: ./io_spike [$trial] [$gc] [$gi] [$synaptic_input] [$threshold]
 *
 * Input:
 *    + [$trial]: the trial index used to name the output files
 *    + [$nstep]: number of integrations, in 0.5 ms step size (e.g, nstep = 100000 ~ 50 s)
 *    + [$gc]: gap-junctional conductance
 *    + [$gi]: inhibitory conductance
 *    + [$synaptic_input]: synaptic input constant
 *    + [$threshold]: threshold for spike detection
 * Output: <io_spike_$trial_$neuron.dat> ($neuron = 1..9)
 * Note: for other parameters see <io_model.hpp>
 * History
 *    + March 29, 2017 - Huu Hoang - version 1.0
 */


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>

#include "io_model.hpp"

typedef boost::array< double , NCELL > ncell_type;
typedef boost::array< double , NCELL*CELLDIM > state_type2;

int main( int argc , char **argv ){
    if (argc < 7) {
        cout << "Error: do not enough required parameters" << endl;
        cout << "Usage: ./io_spike [$trial] [$nstep] [$gc] [$gi] [$synaptic_input] [$threshold]" << endl;
        cout << "\t" << "+ [$trial]: the trial index used to name the output files" << endl;
        cout << "\t" << "+ [$nstep]: number of integrations, in 0.5 ms step size (e.g, nstep = 100000 ~ 50 s)" << endl;
        cout << "\t" << "+ [$gc]: gap-junctional conductanc" << endl;
        cout << "\t" << "+ [$gi]: inhibitory conductance" << endl;
        cout << "\t" << "+ [$synaptic_input]: synaptic input" << endl;
        cout << "\t" << "+ [$threshold]: threshold for spike detection" << endl;
        
        return -1;
    }
    
    // parse input parameters
    size_t trial = atoi(argv[1]);
    size_t nstep = atoi(argv[2]);
    gc = atof(argv[3]);
    gsyn_in = atof(argv[4]);
    synaptic_input = atof(argv[5]);
    double threshold = atof(argv[6]);
    
    state_type2 v;
    ofstream ofs[NCELL];
    
    ncell_type prev_vs;
    fill( prev_vs.begin() , prev_vs.end() , -75.0 );
    
    srand( time( NULL));
    for(size_t c = 0; c<NCELL; c++) {
        ostringstream sstr;
        sstr << "io_spike_" << trial << "_" << c+1 << ".dat";
        
        ofs[c].open(sstr.str().c_str());
        
        for( size_t i = 0; i < CELLDIM; i++) {
            v[c*CELLDIM+i] = default_state[i] + ( 0.001*rand() ) / RAND_MAX;
        }
    }
    
    runge_kutta4< state_type2 , double , state_type2 , double , range_algebra > rk4;
    
    clock_t begin = clock();
    
    // perform 10000 transient steps
    integrate_n_steps( rk4 , io_model() , v , 0.0 , dt , 10000 );
    cout << "Perform 10000 transient steps ... done" << endl;
    //]
    
    double t = 0.0;
    size_t show_nstep = nstep/10;
    
    for (size_t n=0; n<nstep; ++n) {
        t = integrate_n_steps( rk4 , io_model() , v , t , dt , 10);
        
        for (size_t c=0; c<NCELL; c++) {
            if (vs(v,c) - prev_vs[c] >= threshold) ofs[c] << t << "\n";
            prev_vs[c] = vs(v,c);
        }
        
        if (n>0 && n%show_nstep==0) cout << "Generating IO spikes ... " << (n*100/nstep) << "%" << endl;
    }
    
    for (size_t c=0; c<NCELL; c++) { ofs[c].close(); }
    
    clock_t end = clock();
    double elapsed_sec = double (end - begin) / CLOCKS_PER_SEC;
    cout << "Total time [s]: " << elapsed_sec << endl;
    
    return 0;
}
