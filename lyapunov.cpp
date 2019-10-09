/*
 * lyapunov.cpp
 *
 * This program computes lyapunov exponents of the IO model
 * Using: ./lyapunov [$trial] [$gc] [$gi] [$synaptic_input]
 *
 * Input:
 *    + [$trial]: the trial index used to name the output files
 *    + [$nstep]: number of integrations, in 0.5 ms step size (e.g, nstep = 100000 ~ 50 s)
 *    + [$gc]: gap-junctional conductance
 *    + [$gi]: inhibitory conductance
 *    + [$synaptic_input]: synaptic input constant
 * Output: <lyapunov_$trial.dat>
 * Note: for other parameters see <io_model.hpp>
 * History
 *    + March 29, 2017 - Huu Hoang - version 1.0
 */

#include <iostream>
#include <fstream>
#include <ctime>

#include "io_model.hpp"

int main( int argc , char **argv )
{
    if (argc < 6) {
        cout << "Error: do not enough required parameters" << endl;
        cout << "Usage: ./io_spike [$trial] [$nstep] [$gc] [$gi] [$synaptic_input]" << endl;
        cout << "\t" << "+ [$trial]: the trial index used to name the output files" << endl;
        cout << "\t" << "+ [$nstep]: number of integrations, in 0.5 ms step size (e.g, nstep = 100000 ~ 50 s)" << endl;
        cout << "\t" << "+ [$gc]: gap-junctional conductanc" << endl;
        cout << "\t" << "+ [$gi]: inhibitory conductance" << endl;
        cout << "\t" << "+ [$synaptic_input]: synaptic input constant" << endl;
        
        return -1;
    }
    
    // parse input parameters
    size_t trial = atoi(argv[1]);
    size_t nstep = atoi(argv[2]);
    gc = atof(argv[3]);
    gsyn_in = atof(argv[4]);
    synaptic_input = atof(argv[5]);
    
    state_type x;
    lyap_type lyap;
    
    srand( time(NULL) );
    for(size_t c = 0; c<NCELL; c++) {
        for( size_t i = 0; i < CELLDIM; i++) {
            x[c*CELLDIM+i] = default_state[i] + ( 0.001*rand() ) / RAND_MAX;
        }
    }
    
    //[ integrate_transients_with_range
    runge_kutta4< state_type , double , state_type , double , range_algebra > rk4;
    
    // perform 10000 transient steps
    integrate_n_steps( rk4 , io_model() , std::make_pair( x.begin() , x.begin() + num_dim ) , 0.0 , dt , 10000 );
    cout << "Perform 10000 transient steps ... done" << endl;
    //]
    
    //[ lyapunov_full_code
    clock_t begin = clock();
    
    fill( x.begin()+num_dim , x.end() , 0.0 );
    for( size_t i=0 ; i<num_of_lyap ; ++i ) x[num_dim+num_dim*i+i] = 1.0;
    fill( lyap.begin() , lyap.end() , 0.0 );
    
    ostringstream sstr;
    ofstream ofs_lyap;
    
    sstr << "lyapunov_" << trial << ".dat";
    ofs_lyap.open(sstr.str().c_str());
    
    double t = 0.0;
    size_t show_nstep = nstep/10;
    for (size_t n=0; n<nstep; n++)
    {
        t = integrate_n_steps( rk4 , io_model_with_lyap , x , t , dt , 10);
        gram_schmidt< num_of_lyap >( x , lyap , num_dim );
        
        ofs_lyap << t;
        for (size_t i=0; i<num_of_lyap; ++i) ofs_lyap << "\t" << lyap[i]/t;
        ofs_lyap << "\n";
        
        if (n>0 & n%show_nstep==0) cout << "Calculating Lyapunov exponents ... " << (n*100/nstep) << "%" << endl;
    }
    cout << "Calculating Lyapunov exponents ... done >> [output]: " << sstr.str() << endl;
    //]
    
    ofs_lyap.close();
    
    clock_t end = clock();
    double elapsed_sec = double (end - begin) / CLOCKS_PER_SEC;
    cout << "Running time [s]: " << elapsed_sec << endl;
    
    return 0;
}
