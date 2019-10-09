/*
 * io_raw_noise.cpp
 *
 * This program generates soma membrane potential of 3x3 IO neurons with synaptic noise
 * Using: ./io_raw_noise [$trial] [$gc] [$gi] [$synaptic_input]
 *
 * Input:
 *    + [$trial]: the trial index used to name the output files
 *    + [$nstep]: number of integrations, in 0.5 ms step size (e.g, nstep = 100000 ~ 50 s)
 *    + [$gc]: gap-junctional conductance
 *    + [$gi]: inhibitory conductance
 *    + [$freq]: frequency of synaptic noise
 * Output: <io_raw_$trial_$neuron.dat> ($neuron = 1..9)
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
#include "synaptic_input.hpp"

typedef boost::array< double , NCELL > ncell_type;
typedef boost::array< double , NCELL*CELLDIM > state_type2;

int main( int argc , char **argv ){
    if (argc < 6) {
        cout << "Error: do not enough required parameters" << endl;
        cout << "Usage: ./io_raw_noise [$trial] [$nstep] [$gc] [$gi] [$freq]" << endl;
        cout << "\t" << "+ [$trial]: the trial index used to name the output files" << endl;
        cout << "\t" << "+ [$nstep]: number of integrations, in 0.5 ms step size (e.g, nstep = 100000 ~ 50 s)" << endl;
        cout << "\t" << "+ [$gc]: gap-junctional conductanc" << endl;
        cout << "\t" << "+ [$gi]: inhibitory conductance" << endl;
        cout << "\t" << "+ [$freq]: frequency of synaptic noise" << endl;
        
        return -1;
    }
    
    // parse input parameters
    size_t trial = atoi(argv[1]);
    size_t nstep = atoi(argv[2]);
    gc = atof(argv[3]);
    gsyn_in = atof(argv[4]);
    double freq  = atof(argv[5]);
    
    state_type2 v;
    ofstream ofs[NCELL];
    
    cout << "Initialize the system and generate synaptic input ... ";
    srand( time( NULL));
    for(size_t c = 0; c<NCELL; c++) {
        ostringstream sstr;
        sstr << "io_raw_" << trial << "_" << c+1 << ".dat";
        
        ofs[c].open(sstr.str().c_str());
        
        for( size_t i = 0; i < CELLDIM; i++) {
            v[c*CELLDIM+i] = default_state[i] + ( 0.001*rand() ) / RAND_MAX;
        }
        
        for( size_t i=0; i<3; i++ ) {
            generate_synaptic_input(excitatory_noise[c][i], freq, nsynapse[i], nstep);
	    generate_synaptic_input(inhibitory_noise[c][i], freq, nsynapse[i], nstep);
        }
    }
    cout << "done" << endl;
    
    runge_kutta4< state_type2 , double , state_type2 , double , range_algebra > rk4;
    
    clock_t begin = clock();
    
    // perform 10000 transient steps
    integrate_n_steps( rk4 , io_model() , v , 0.0 , dt , 10000 );
    cout << "Perform 10000 transient steps ... done" << endl;
    //]
    
    double t = 0.0;
    size_t show_nstep = nstep/10;
    
    for (size_t n=0; n<nstep; ++n) {
        t = integrate_n_steps( rk4 , io_model_with_noise() , v , t , dt , 10);
        
        for (size_t c=0; c<NCELL; c++) {
            ofs[c] << vs(v,c) << "\n";
        }
        
        if (n>0 & n%show_nstep==0) cout << "Running IO model ... " << (n*100/nstep) << "%" << endl;
    }
    
    for (size_t c=0; c<NCELL; c++) { ofs[c].close(); }
    
    clock_t end = clock();
    double elapsed_sec = double (end - begin) / CLOCKS_PER_SEC;
    cout << "Total time [s]: " << elapsed_sec << endl;
    
    return 0;
}
