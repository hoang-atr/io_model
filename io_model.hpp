//
//  io_net.h
//  io_simulation
//
//  Created by Huu Hoang on 4/20/13.
//  Copyright (c) 2013 __Ritsumeikan__. All rights reserved.
//
#ifndef _io_model_hpp
#define _io_model_hpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#include "gram_schmidt.hpp"

using namespace std;
using namespace boost::numeric::odeint;

#define     NCELL        9
#define     CELLDIM     14

/*******************************************************/
/*****************  somatic currents ******************/
/*******************************************************/
/*  aux functions soma Calcium low threshold (Manor)*/
#define	hcalinf(v)	(1./(1+exp((v+85.5)/8.5)))
#define	hcaltau(v) 	((20.*exp((v+160.)/30.)/(1.+exp((v+84.)/7.3)))+35.)
#define	scalinf(v)	(1./(1.+exp(-(v+61.)/4.2)))
#define	scaltau(v)	(1.0)  /*PNAS*/
/*#define scaltau(v)      (5.0)*/

#define hcalinfdot(v)   ((-1/8.5)*exp((v+85.5)/8.5)/((1+exp((v+85.5)/8.5))*(1+exp((v+85.5)/8.5))))
#define hcaltaudot(v)   ( ( (1.+exp((v+84.)/7.3)) * (20/30)*exp((v+160.)/30.) - (20/7.3)*exp((v+160.)/30.)*exp((v+84.)/7.3) ) / ( (1.+exp((v+84.)/7.3))*(1.+exp((v+84.)/7.3)) ) )
#define scalinfdot(v)   ( (1/4.2)*exp(-(v+61.)/4.2) / ( (1.+exp(-(v+61.)/4.2))*(1.+exp(-(v+61.)/4.2)) ) )

/* gna activation */
/*#define	alpham(v)	(0.1*(v+41.)/(1.-exp(-(v+41.)/10 ))) */
#define	alpham(v)	(0.1*(v+48.)/(1.-exp(-(v+48.)/3)))
#define	betam(v)   	(9*exp(-(v+66)/20.))
#define	minf(v)		(alpham(v)/(alpham(v)+betam(v)))

#define alphamdot(v) ( (1.-exp(-(v+48.)/3))*0.1 - (exp(-(v+48.)/3)/3)*0.1*(v+48.) )/((1.-exp(-(v+48.)/3))*(1.-exp(-(v+48.)/3)) )
#define	betamdot(v)  ( (-9/20)*exp(-(v+66)/20.) )
#define minfdot(v)  ( ((alpham(v)+betam(v))*alphamdot(v) - alpham(v)*(alphamdot(v)+betamdot(v))) / ((alpham(v)+betam(v))*(alpham(v)+betam(v))) )

/* gna inactivation */
#define	alphah(v) 	(5.*exp(-(v+60)/15))
#define	betah(v)	((v+50)/(1.-exp(-(v+50)/10)))
#define	hinf(v)		((alphah(v)/(alphah(v)+betah(v))))
#define	htau(v)		(300/(alphah(v)+betah(v)))

#define alphahdot(v) ( (-5/15)*exp(-(v+60)/15) )
#define betahdot(v) ( ((1.-exp(-(v+50)/10)) - (v+50)*exp(-(v+50)/10)/10) / ((1.-exp(-(v+50)/10))*(1.-exp(-(v+50)/10))) )
#define hinfdot(v) ( ((alphah(v)+betah(v))*alphahdot(v) - alphah(v)*(alphahdot(v)+betahdot(v))) / ((alphah(v)+betah(v))*(alphah(v)+betah(v))) )
#define htaudot(v) ( -300*(alphahdot(v)+betahdot(v)) / ((alphah(v)+betah(v))*(alphah(v)+betah(v))) )

/* gk activation (from rush and rinzel)*/
#define	alphan(v) 	(v+41)/(1.-exp(-(v+41)/10.))
#define	betan(v)	(12.5*exp(-(v+51)/80.0))
#define	ninf(v)		(alphan(v)/(alphan(v)+betan(v)))
#define	ntau(v)		(5/(alphan(v)+betan(v)))

#define alphandot(v)  (  ((1.-exp(-(v+41)/10)) - (v+41)*exp(-(v+41)/10)/10) / ((1.-exp(-(v+41)/10))*(1.-exp(-(v+41)/10))) )
#define betandot(v)     (-12.5/80.0)*exp(-(v+51)/80.0)
#define ninfdot(v)      ((alphan(v)+betan(v))*alphandot(v) - alphan(v)*(alphandot(v)+betandot(v)))/((alphan(v)+betan(v))*(alphan(v)+betan(v)))
#define ntaudot(v)      (-5*(alphandot(v)+betandot(v))/((alphan(v)+betan(v))*(alphan(v)+betan(v))))

/*  gkar anomalous rectifier activation (McCormick)*/
#define	arinf(v)	(1./(1.+exp((v+75)/5.5)))
#define	artau(v)	(1./(exp(-0.086*v-14.6)+exp(0.070*v-1.87)))

#define arinfdot(v) ((-1/5.5)*exp((v+75)/5.5)) / ((1.+exp((v+75)/5.5))*(1.+exp((v+75)/5.5)))
#define artaudot(v) ( -(-0.086*exp(-0.086*v-14.6) + 0.070*exp(0.070*v-1.87)) / ((exp(-0.086*v-14.6)+exp(0.070*v-1.87))*(exp(-0.086*v-14.6)+exp(0.070*v-1.87))) )
/*******************************************************/
/***************** dendritic currents ******************/
/*******************************************************/

/*  aux functions dendrite Calcium high threshold (Traub et al.)*/
#define	alphacah(v)	(1.6/(1.+exp(-(v-5.)/13.9)))
#define	betacah(v)	(0.02*(v+8.5)/(exp((v+8.5)/5.) -1.))
#define	infscah(v)	(alphacah(v)/(alphacah(v)+betacah(v)))
#define	tauscah(v)	(5./(alphacah(v)+betacah(v)))

#define alphacahdot(v)  (1.6/13.9)*exp(-(v-5.)/13.9) / ((1.+exp(-(v-5.)/13.9))*(1.+exp(-(v-5.)/13.9)))
#define betacahdot(v)   ( 0.02*(exp((v+8.5)/5.) -1.) - (0.02/5.)*(v+8.5)*exp((v+8.5)/5.) ) / ((exp((v+8.5)/5.) -1.)*(exp((v+8.5)/5.) -1.))
#define infscahdot(v)   (alphacahdot(v)*(alphacah(v)+betacah(v)) - alphacah(v)*(alphacahdot(v)+betacahdot(v))) / ((alphacah(v)+betacah(v))*(alphacah(v)+betacah(v)))
#define tauscahdot(v)   ( -5*(alphacahdot(v)+betacahdot(v)) / ((alphacah(v)+betacah(v))*(alphacah(v)+betacah(v))) )

/*  gk Ca dependant  (Traub et al. but different time constant)*/
#define min(a,b) ({ typeof (a) _a = (a); typeof (b) _b = (b); _a < _b ? _a : _b; })
#define	alpha_m(ca_d)	(min(0.00002*ca_d,0.01))
#define	beta_m			(0.015)
#define	mk_ca_inf(ca_d)	(alpha_m(ca_d)/(alpha_m(ca_d)+beta_m))
#define	mk_ca_tau(ca_d)	(1./(alpha_m(ca_d)+beta_m))

#define mk_ca_infdot(ca_d) (0.00002*ca_d < 0.01 ? (0.00002*beta_m)/((alpha_m(ca_d)+beta_m)*(alpha_m(ca_d)+beta_m)) : 0)
#define mk_ca_taudot(ca_d) (0.00002*ca_d < 0.01 ? -0.00002/((alpha_m(ca_d)+beta_m)*(alpha_m(ca_d)+beta_m)) : 0)

/****
 *	Cellular Dynamics
 ****/

/* state of soma */
#define	vs(v,c)		v[c*CELLDIM+0]
#define	hcal(v,c)	v[c*CELLDIM+1]
#define	scal(v,c)	v[c*CELLDIM+2]
#define	ns(v,c)		v[c*CELLDIM+3]
#define	h(v,c)		v[c*CELLDIM+4]
#define	ar(v,c)		v[c*CELLDIM+5]

/* state of dendrites */
#define	vd(v,c)		v[c*CELLDIM+6]
#define	scah(v,c)	v[c*CELLDIM+7]
#define	ca_d(v,c)   v[c*CELLDIM+8]
#define	mk_cad(v,c) v[c*CELLDIM+9]

/* state of spine*/
#define vsp1(v,c)   v[c*CELLDIM+10]
#define vsp2(v,c)   v[c*CELLDIM+11]
#define vsp3(v,c)   v[c*CELLDIM+12]
#define vsp4(v,c)   v[c*CELLDIM+13]

/* coupling and surface ratio soma / dendrite */
double const    ginter      =   0.13;
double const    p           =   0.14;
double const    ginter_sp   =   0.1;
double const    pp          =   0.05;

/* conductance leak */
double const    gl          =   0.015;
double const    vl          =   -10;

/* current Calcium */
double const    gcal        =   2.;
double const    gcah        =   4;
double const    vca         =   120;
double const    ca_tau      =   0.02;

/* conductance sodium  currents */
double const    gna         =   110;
double const    vna         =   55;

/* conductance potassium */
double const 	gk          =   18.;
double const    gk_cad      =   35.;
double const    gk_ar       =   0.15;
double const    vk          =   -75.;

/* synaptic parameters */
double const    vsynapex    =   -10.0;
double const    vsynapen_sp =   -10.0;
double const    vsynapin    =   -70.0;
double const    vsynapin_sp =   -70.0;

double const    dt          =   0.05;

double const default_state[CELLDIM] = {       // default states
    /* vs */		-72.55624,
    /* hcal */		0.1221349,
    /* scal */		0.0623159,
    /* ns */		0.0869847,
    /* h */			0.8096066,
    /* ar */		0.0737836,
    /* vd */		-72.52936,
    /* scah */		0.0046278,
    /* ca_d */		1.94452,
    /* mk_cad */    0.0037291,
    /* vsp1 */      -72.0,
    /* vsp2 */      -72.0,
    /* vsp3 */      -72.0,
    /* vsp4 */      -72.0 };
    
    double gc = 0.5;
    double gsyn_in = 0.01;
    double synaptic_input = 0.2;
    double Iapp = 0;
    
    vector<double> inhibitory_noise[NCELL][3];
    vector<double> excitatory_noise[NCELL][3];
    
    const int nsynapse[3] = {10, 80, 10};
    
    const double gsyn_en = .03;
    const double dend_ratio = 0.1;
    const double gc_table[18] = {0.3335, 0.2878, 0.2599, -0.5198, -0.0558,
    -0.4399, 0.6788, 0.7740, 0.7536, -0.5860, -0.3925, -0.2474, 0.0193,
    -0.1418, 0.1531, -0.7392, -0.8206, -0.1069 };
    const double gcal_table[9] = {-0.4071, -0.9254, 0.1864, 0.0502, -0.7294,
    0.3381, -0.4509, -0.3243, 0.8301};
    
    const size_t num_dim = NCELL*CELLDIM;
    const size_t num_of_lyap = NCELL*CELLDIM;
    const size_t M = num_dim + num_dim*num_dim;
    
    typedef boost::array< double , M > state_type;
    typedef boost::array< double , num_of_lyap > lyap_type;
    
    struct io_model
    {
        template< class State , class Deriv >
                void operator()( const State &x_ , Deriv &dxdt_ , double t ) const
        {
            typename boost::range_iterator< const State >::type v = boost::begin( x_ );
            typename boost::range_iterator< Deriv >::type v_dot = boost::begin( dxdt_ );
            
            double	Ical, Iap, Idend, Ils, Ikars, Isynapen_s, Isynapin_s;
            double	Icah, Ikcad, Isoma, Ild, Isynapen_d, Isynapin_d;
            double  Idend_sp1, Idend_sp2, Idend_sp3, Idend_sp4;
            double  Ispin1, Ispin2, Ispin3, Ispin4;
            double  Ilsp1, Ilsp2, Ilsp3, Ilsp4;
            double	Isynapen1, Isynapen2, Isynapen3, Isynapen4;
            double	Isynapin1, Isynapin2, Isynapin3, Isynapin4;
            double	Icoupl1, Icoupl2, Icoupl3, Icoupl4;
            
            for( size_t cell_nb = 0; cell_nb < NCELL; cell_nb++){
                switch(cell_nb){
                    
                    case 0:
                        Icoupl1 = (1+0.2*gc_table[0])*gc*(vsp1(v,0)-vsp3(v,1));
                        Icoupl2 = (1+0.2*gc_table[2])*gc*(vsp2(v,0)-vsp4(v,3));
                        Icoupl3 = (1+0.2*gc_table[12])*gc*(vsp3(v,0)-vsp1(v,2));
                        Icoupl4 = (1+0.2*gc_table[13])*gc*(vsp4(v,0)-vsp2(v,6));
                        break;
                        
                    case 1:
                        Icoupl1 = (1+0.2*gc_table[1])*gc*(vsp1(v,1)-vsp3(v,2));
                        Icoupl2 = (1+0.2*gc_table[3])*gc*(vsp2(v,1)-vsp4(v,4));
                        Icoupl3 = (1+0.2*gc_table[0])*gc*(vsp3(v,1)-vsp1(v,0));
                        Icoupl4 = (1+0.2*gc_table[17])*gc*(vsp4(v,1)-vsp2(v,7));
                        break;
                        
                    case 2:
                        Icoupl1 = (1+0.2*gc_table[12])*gc*(vsp1(v,2)-vsp3(v,0));
                        Icoupl2 = (1+0.2*gc_table[4])*gc*(vsp2(v,2)-vsp4(v,5));
                        Icoupl3 = (1+0.2*gc_table[1])*gc*(vsp3(v,2)-vsp1(v,1));
                        Icoupl4 = (1+0.2*gc_table[14])*gc*(vsp4(v,2)-vsp2(v,8));
                        break;
                        
                    case 3:
                        Icoupl1 = (1+0.2*gc_table[5])*gc*(vsp1(v,3)-vsp3(v,4));
                        Icoupl2 = (1+0.2*gc_table[7])*gc*(vsp2(v,3)-vsp4(v,6));
                        Icoupl3 = (1+0.2*gc_table[16])*gc*(vsp3(v,3)-vsp1(v,5));
                        Icoupl4 = (1+0.2*gc_table[2])*gc*(vsp4(v,3)-vsp2(v,0));
                        break;
                        
                    case 4:
                        Icoupl1 = (1+0.2*gc_table[6])*gc*(vsp1(v,4)-vsp3(v,5));
                        Icoupl2 = (1+0.2*gc_table[8])*gc*(vsp2(v,4)-vsp4(v,7));
                        Icoupl3 = (1+0.2*gc_table[5])*gc*(vsp3(v,4)-vsp1(v,3));
                        Icoupl4 = (1+0.2*gc_table[3])*gc*(vsp4(v,4)-vsp2(v,1));
                        break;
                        
                    case 5:
                        Icoupl1 = (1+0.2*gc_table[16])*gc*(vsp1(v,5)-vsp3(v,3));
                        Icoupl2 = (1+0.2*gc_table[9])*gc*(vsp2(v,5)-vsp4(v,8));
                        Icoupl3 = (1+0.2*gc_table[6])*gc*(vsp3(v,5)-vsp1(v,4));
                        Icoupl4 = (1+0.2*gc_table[4])*gc*(vsp4(v,5)-vsp2(v,2));
                        break;
                        
                    case 6:
                        Icoupl1 = (1+0.2*gc_table[10])*gc*(vsp1(v,6)-vsp3(v,7));
                        Icoupl2 = (1+0.2*gc_table[13])*gc*(vsp2(v,6)-vsp4(v,0));
                        Icoupl3 = (1+0.2*gc_table[15])*gc*(vsp3(v,6)-vsp1(v,8));
                        Icoupl4 = (1+0.2*gc_table[7])*gc*(vsp4(v,6)-vsp2(v,3));
                        break;
                        
                    case 7:
                        Icoupl1 = (1+0.2*gc_table[11])*gc*(vsp1(v,7)-vsp3(v,8));
                        Icoupl2 = (1+0.2*gc_table[17])*gc*(vsp2(v,7)-vsp4(v,1));
                        Icoupl3 = (1+0.2*gc_table[10])*gc*(vsp3(v,7)-vsp1(v,6));
                        Icoupl4 = (1+0.2*gc_table[8])*gc*(vsp4(v,7)-vsp2(v,4));
                        break;
                        
                    case 8:
                        Icoupl1 = (1+0.2*gc_table[15])*gc*(vsp1(v,8)-vsp3(v,6));
                        Icoupl2 = (1+0.2*gc_table[14])*gc*(vsp2(v,8)-vsp4(v,2));
                        Icoupl3 = (1+0.2*gc_table[11])*gc*(vsp3(v,8)-vsp1(v,7));
                        Icoupl4 = (1+0.2*gc_table[9])*gc*(vsp4(v,8)-vsp2(v,5));
                        break;
                }
                
                /* currents of somas */
                Ical = gcal*(1.02+0.05*gcal_table[cell_nb])*(scal(v,cell_nb)*scal(v,cell_nb)*scal(v,cell_nb))*hcal(v,cell_nb)*(vs(v,cell_nb)-vca);
                Iap = gna*(minf(vs(v,cell_nb))*minf(vs(v,cell_nb))*minf(vs(v,cell_nb)))*h(v,cell_nb)*(vs(v,cell_nb)-vna)+gk*(ns(v,cell_nb)*ns(v,cell_nb)*ns(v,cell_nb)*ns(v,cell_nb))*(vs(v,cell_nb)-vk);
                Idend = (ginter/p)*(vs(v,cell_nb)-vd(v,cell_nb));
                Ils = gl*(vs(v,cell_nb)-vl);
                Ikars = gk_ar*ar(v,cell_nb)*(vs(v,cell_nb)+43);
                
                Isynapen_s = (gsyn_en*synaptic_input)*(vs(v,cell_nb)-vsynapex);
                Isynapin_s = (gsyn_in*synaptic_input)*(vs(v,cell_nb)-vsynapin);
                
                /* currents of dendrites */
                Icah = gcah*scah(v,cell_nb)*scah(v,cell_nb)*(vd(v,cell_nb)-vca);
                Ikcad = gk_cad*mk_cad(v,cell_nb)*(vd(v,cell_nb)-vk);
                Isoma = (ginter/(1.-p-pp))*(vd(v,cell_nb)-vs(v,cell_nb));
                Ild = gl*(vd(v,cell_nb)-vl);
                
                Isynapen_d = (gsyn_en*synaptic_input)*(vd(v,cell_nb)-vsynapex);
                Isynapin_d = (gsyn_in*synaptic_input*dend_ratio)*(vd(v,cell_nb)-vsynapin);
                
                /* currents of spines */
                Ispin1 = (ginter_sp/(1.-p-pp))*(vd(v,cell_nb)-vsp1(v,cell_nb));
                Ispin2 = (ginter_sp/(1.-p-pp))*(vd(v,cell_nb)-vsp2(v,cell_nb));
                Ispin3 = (ginter_sp/(1.-p-pp))*(vd(v,cell_nb)-vsp3(v,cell_nb));
                Ispin4 = (ginter_sp/(1.-p-pp))*(vd(v,cell_nb)-vsp4(v,cell_nb));
                
                Isynapen1 = (gsyn_en*synaptic_input)*(vsp1(v,cell_nb)-vsynapen_sp);
                Isynapen2 = (gsyn_en*synaptic_input)*(vsp2(v,cell_nb)-vsynapen_sp);
                Isynapen3 = (gsyn_en*synaptic_input)*(vsp3(v,cell_nb)-vsynapen_sp);
                Isynapen4 = (gsyn_en*synaptic_input)*(vsp4(v,cell_nb)-vsynapen_sp);
                Isynapin1 = (gsyn_in*synaptic_input)*(vsp1(v,cell_nb)-vsynapin_sp);
                Isynapin2 = (gsyn_in*synaptic_input)*(vsp2(v,cell_nb)-vsynapin_sp);
                Isynapin3 = (gsyn_in*synaptic_input)*(vsp3(v,cell_nb)-vsynapin_sp);
                Isynapin4 = (gsyn_in*synaptic_input)*(vsp4(v,cell_nb)-vsynapin_sp);
                
                Idend_sp1 = (ginter_sp/(0.25*pp))*(vsp1(v,cell_nb)-vd(v,cell_nb));
                Idend_sp2 = (ginter_sp/(0.25*pp))*(vsp2(v,cell_nb)-vd(v,cell_nb));
                Idend_sp3 = (ginter_sp/(0.25*pp))*(vsp3(v,cell_nb)-vd(v,cell_nb));
                Idend_sp4 = (ginter_sp/(0.25*pp))*(vsp4(v,cell_nb)-vd(v,cell_nb));
                
                Ilsp1 = gl*(vsp1(v,cell_nb)-vl);
                Ilsp2 = gl*(vsp2(v,cell_nb)-vl);
                Ilsp3 = gl*(vsp3(v,cell_nb)-vl);
                Ilsp4 = gl*(vsp4(v,cell_nb)-vl);
                
                /* diff eq soma */
                vs(v_dot,cell_nb) =(-Ical-Iap-Ils-Ikars-Idend-Isynapin_s-Isynapen_s);
                hcal(v_dot,cell_nb) = (hcalinf(vs(v,cell_nb))-hcal(v,cell_nb))/(hcaltau(vs(v,cell_nb)));
                scal(v_dot,cell_nb) = (scalinf(vs(v,cell_nb))-scal(v,cell_nb))/(scaltau(vs(v,cell_nb)));
                ns(v_dot,cell_nb) = (ninf(vs(v,cell_nb))-ns(v,cell_nb))/(ntau(vs(v,cell_nb)));
                h(v_dot,cell_nb) = (hinf(vs(v,cell_nb))-h(v,cell_nb))/(htau(vs(v,cell_nb)));
                ar(v_dot,cell_nb) = (arinf(vs(v,cell_nb))-ar(v,cell_nb))/(artau(vs(v,cell_nb)));
                
                /* diff eq dendrite */
                vd(v_dot,cell_nb) = (-Icah-Ikcad-Ild-Isoma-Isynapen_d-Isynapin_d -Ispin1-Ispin2-Ispin3-Ispin4);
                scah(v_dot,cell_nb) = (infscah(vd(v,cell_nb))-scah(v,cell_nb))/(tauscah(vd(v,cell_nb)));
                ca_d(v_dot,cell_nb) = -1.01*Icah-0.02*ca_d(v,cell_nb);
                mk_cad(v_dot,cell_nb) = (mk_ca_inf(ca_d(v,cell_nb))-mk_cad(v,cell_nb))/(mk_ca_tau(ca_d(v,cell_nb)));
                
                /* diff eq spine */
                vsp1(v_dot,cell_nb) = (-Ilsp1-Idend_sp1-Isynapen1-Isynapin1-Icoupl1);
                vsp2(v_dot,cell_nb) = (-Ilsp2-Idend_sp2-Isynapen2-Isynapin2-Icoupl2);
                vsp3(v_dot,cell_nb) = (-Ilsp3-Idend_sp3-Isynapen3-Isynapin3-Icoupl3);
                vsp4(v_dot,cell_nb) = (-Ilsp4-Idend_sp4-Isynapen4-Isynapin4-Icoupl4);
            }
        }
    };
    
    struct io_model_cc
    {
        template< class State , class Deriv >
                void operator()( const State &x_ , Deriv &dxdt_ , double t ) const
        {
            typename boost::range_iterator< const State >::type v = boost::begin( x_ );
            typename boost::range_iterator< Deriv >::type v_dot = boost::begin( dxdt_ );
            
            double	Ical, Iap, Idend, Ils, Ikars, Isynapen_s, Isynapin_s;
            double	Icah, Ikcad, Isoma, Ild, Isynapen_d, Isynapin_d;
            double  Idend_sp1, Idend_sp2, Idend_sp3, Idend_sp4;
            double  Ispin1, Ispin2, Ispin3, Ispin4;
            double  Ilsp1, Ilsp2, Ilsp3, Ilsp4;
            double	Isynapen1, Isynapen2, Isynapen3, Isynapen4;
            double	Isynapin1, Isynapin2, Isynapin3, Isynapin4;
            double	Icoupl1, Icoupl2, Icoupl3, Icoupl4;
            const double Ihyp = -1;

	    int tidx = (int) t/0.5;
            
            for( size_t cell_nb = 0; cell_nb < NCELL; cell_nb++){
                switch(cell_nb){
                    
                    case 0:
                        Icoupl1 = (1+0.2*gc_table[0])*gc*(vsp1(v,0)-vsp3(v,1));
                        Icoupl2 = (1+0.2*gc_table[2])*gc*(vsp2(v,0)-vsp4(v,3));
                        Icoupl3 = (1+0.2*gc_table[12])*gc*(vsp3(v,0)-vsp1(v,2));
                        Icoupl4 = (1+0.2*gc_table[13])*gc*(vsp4(v,0)-vsp2(v,6));
                        break;
                        
                    case 1:
                        Icoupl1 = (1+0.2*gc_table[1])*gc*(vsp1(v,1)-vsp3(v,2));
                        Icoupl2 = (1+0.2*gc_table[3])*gc*(vsp2(v,1)-vsp4(v,4));
                        Icoupl3 = (1+0.2*gc_table[0])*gc*(vsp3(v,1)-vsp1(v,0));
                        Icoupl4 = (1+0.2*gc_table[17])*gc*(vsp4(v,1)-vsp2(v,7));
                        break;
                        
                    case 2:
                        Icoupl1 = (1+0.2*gc_table[12])*gc*(vsp1(v,2)-vsp3(v,0));
                        Icoupl2 = (1+0.2*gc_table[4])*gc*(vsp2(v,2)-vsp4(v,5));
                        Icoupl3 = (1+0.2*gc_table[1])*gc*(vsp3(v,2)-vsp1(v,1));
                        Icoupl4 = (1+0.2*gc_table[14])*gc*(vsp4(v,2)-vsp2(v,8));
                        break;
                        
                    case 3:
                        Icoupl1 = (1+0.2*gc_table[5])*gc*(vsp1(v,3)-vsp3(v,4));
                        Icoupl2 = (1+0.2*gc_table[7])*gc*(vsp2(v,3)-vsp4(v,6));
                        Icoupl3 = (1+0.2*gc_table[16])*gc*(vsp3(v,3)-vsp1(v,5));
                        Icoupl4 = (1+0.2*gc_table[2])*gc*(vsp4(v,3)-vsp2(v,0));
                        break;
                        
                    case 4:
                        Icoupl1 = (1+0.2*gc_table[6])*gc*(vsp1(v,4)-vsp3(v,5));
                        Icoupl2 = (1+0.2*gc_table[8])*gc*(vsp2(v,4)-vsp4(v,7));
                        Icoupl3 = (1+0.2*gc_table[5])*gc*(vsp3(v,4)-vsp1(v,3));
                        Icoupl4 = (1+0.2*gc_table[3])*gc*(vsp4(v,4)-vsp2(v,1));
                        break;
                        
                    case 5:
                        Icoupl1 = (1+0.2*gc_table[16])*gc*(vsp1(v,5)-vsp3(v,3));
                        Icoupl2 = (1+0.2*gc_table[9])*gc*(vsp2(v,5)-vsp4(v,8));
                        Icoupl3 = (1+0.2*gc_table[6])*gc*(vsp3(v,5)-vsp1(v,4));
                        Icoupl4 = (1+0.2*gc_table[4])*gc*(vsp4(v,5)-vsp2(v,2));
                        break;
                        
                    case 6:
                        Icoupl1 = (1+0.2*gc_table[10])*gc*(vsp1(v,6)-vsp3(v,7));
                        Icoupl2 = (1+0.2*gc_table[13])*gc*(vsp2(v,6)-vsp4(v,0));
                        Icoupl3 = (1+0.2*gc_table[15])*gc*(vsp3(v,6)-vsp1(v,8));
                        Icoupl4 = (1+0.2*gc_table[7])*gc*(vsp4(v,6)-vsp2(v,3));
                        break;
                        
                    case 7:
                        Icoupl1 = (1+0.2*gc_table[11])*gc*(vsp1(v,7)-vsp3(v,8));
                        Icoupl2 = (1+0.2*gc_table[17])*gc*(vsp2(v,7)-vsp4(v,1));
                        Icoupl3 = (1+0.2*gc_table[10])*gc*(vsp3(v,7)-vsp1(v,6));
                        Icoupl4 = (1+0.2*gc_table[8])*gc*(vsp4(v,7)-vsp2(v,4));
                        break;
                        
                    case 8:
                        Icoupl1 = (1+0.2*gc_table[15])*gc*(vsp1(v,8)-vsp3(v,6));
                        Icoupl2 = (1+0.2*gc_table[14])*gc*(vsp2(v,8)-vsp4(v,2));
                        Icoupl3 = (1+0.2*gc_table[11])*gc*(vsp3(v,8)-vsp1(v,7));
                        Icoupl4 = (1+0.2*gc_table[9])*gc*(vsp4(v,8)-vsp2(v,5));
                        break;
                }
                
                /* currents of somas */
                Ical = gcal*(1.02+0.05*gcal_table[cell_nb])*(scal(v,cell_nb)*scal(v,cell_nb)*scal(v,cell_nb))*hcal(v,cell_nb)*(vs(v,cell_nb)-vca);
                Iap = gna*(minf(vs(v,cell_nb))*minf(vs(v,cell_nb))*minf(vs(v,cell_nb)))*h(v,cell_nb)*(vs(v,cell_nb)-vna)+gk*(ns(v,cell_nb)*ns(v,cell_nb)*ns(v,cell_nb)*ns(v,cell_nb))*(vs(v,cell_nb)-vk);
                Idend = (ginter/p)*(vs(v,cell_nb)-vd(v,cell_nb));
                Ils = gl*(vs(v,cell_nb)-vl);
                Ikars = gk_ar*ar(v,cell_nb)*(vs(v,cell_nb)+43);
                
                Isynapen_s = (gsyn_en*excitatory_noise[cell_nb][0][tidx])*(vs(v,cell_nb)-vsynapex);
                Isynapin_s = (gsyn_in*inhibitory_noise[cell_nb][0][tidx])*(vs(v,cell_nb)-vsynapin);
                
                /* currents of dendrites */
                Icah = gcah*scah(v,cell_nb)*scah(v,cell_nb)*(vd(v,cell_nb)-vca);
                Ikcad = gk_cad*mk_cad(v,cell_nb)*(vd(v,cell_nb)-vk);
                Isoma = (ginter/(1.-p-pp))*(vd(v,cell_nb)-vs(v,cell_nb));
                Ild = gl*(vd(v,cell_nb)-vl);
                
                Isynapen_d = (gsyn_en*excitatory_noise[cell_nb][1][tidx])*(vd(v,cell_nb)-vsynapex);
                Isynapin_d = (gsyn_in*inhibitory_noise[cell_nb][1][tidx])*(vd(v,cell_nb)-vsynapin);
                
                /* currents of spines */
                Ispin1 = (ginter_sp/(1.-p-pp))*(vd(v,cell_nb)-vsp1(v,cell_nb));
                Ispin2 = (ginter_sp/(1.-p-pp))*(vd(v,cell_nb)-vsp2(v,cell_nb));
                Ispin3 = (ginter_sp/(1.-p-pp))*(vd(v,cell_nb)-vsp3(v,cell_nb));
                Ispin4 = (ginter_sp/(1.-p-pp))*(vd(v,cell_nb)-vsp4(v,cell_nb));
                
                Isynapen1 = (gsyn_en*excitatory_noise[cell_nb][2][tidx])*(vsp1(v,cell_nb)-vsynapen_sp);
                Isynapen2 = (gsyn_en*excitatory_noise[cell_nb][2][tidx])*(vsp2(v,cell_nb)-vsynapen_sp);
                Isynapen3 = (gsyn_en*excitatory_noise[cell_nb][2][tidx])*(vsp3(v,cell_nb)-vsynapen_sp);
                Isynapen4 = (gsyn_en*excitatory_noise[cell_nb][2][tidx])*(vsp4(v,cell_nb)-vsynapen_sp);
                Isynapin1 = (gsyn_in*inhibitory_noise[cell_nb][2][tidx])*(vsp1(v,cell_nb)-vsynapin_sp);
                Isynapin2 = (gsyn_in*inhibitory_noise[cell_nb][2][tidx])*(vsp2(v,cell_nb)-vsynapin_sp);
                Isynapin3 = (gsyn_in*inhibitory_noise[cell_nb][2][tidx])*(vsp3(v,cell_nb)-vsynapin_sp);
                Isynapin4 = (gsyn_in*inhibitory_noise[cell_nb][2][tidx])*(vsp4(v,cell_nb)-vsynapin_sp);
                
                Idend_sp1 = (ginter_sp/(0.25*pp))*(vsp1(v,cell_nb)-vd(v,cell_nb));
                Idend_sp2 = (ginter_sp/(0.25*pp))*(vsp2(v,cell_nb)-vd(v,cell_nb));
                Idend_sp3 = (ginter_sp/(0.25*pp))*(vsp3(v,cell_nb)-vd(v,cell_nb));
                Idend_sp4 = (ginter_sp/(0.25*pp))*(vsp4(v,cell_nb)-vd(v,cell_nb));
                
                Ilsp1 = gl*(vsp1(v,cell_nb)-vl);
                Ilsp2 = gl*(vsp2(v,cell_nb)-vl);
                Ilsp3 = gl*(vsp3(v,cell_nb)-vl);
                Ilsp4 = gl*(vsp4(v,cell_nb)-vl);
                
                /* diff eq soma */
                if (t>=1500 & t<=2000 & cell_nb==4) {
                    vs(v_dot,cell_nb) =(-Ical-Iap-Ils-Ikars-Idend-Isynapin_s-Isynapen_s)+Ihyp+Iapp;
                } else {
                    vs(v_dot,cell_nb) =(-Ical-Iap-Ils-Ikars-Idend-Isynapin_s-Isynapen_s)+Ihyp;
                }
                
                hcal(v_dot,cell_nb) = (hcalinf(vs(v,cell_nb))-hcal(v,cell_nb))/(hcaltau(vs(v,cell_nb)));
                scal(v_dot,cell_nb) = (scalinf(vs(v,cell_nb))-scal(v,cell_nb))/(scaltau(vs(v,cell_nb)));
                ns(v_dot,cell_nb) = (ninf(vs(v,cell_nb))-ns(v,cell_nb))/(ntau(vs(v,cell_nb)));
                h(v_dot,cell_nb) = (hinf(vs(v,cell_nb))-h(v,cell_nb))/(htau(vs(v,cell_nb)));
                ar(v_dot,cell_nb) = (arinf(vs(v,cell_nb))-ar(v,cell_nb))/(artau(vs(v,cell_nb)));
                
                /* diff eq dendrite */
                vd(v_dot,cell_nb) = (-Icah-Ikcad-Ild-Isoma-Isynapen_d-Isynapin_d -Ispin1-Ispin2-Ispin3-Ispin4)+Ihyp;
                scah(v_dot,cell_nb) = (infscah(vd(v,cell_nb))-scah(v,cell_nb))/(tauscah(vd(v,cell_nb)));
                ca_d(v_dot,cell_nb) = -1.01*Icah-0.02*ca_d(v,cell_nb);
                mk_cad(v_dot,cell_nb) = (mk_ca_inf(ca_d(v,cell_nb))-mk_cad(v,cell_nb))/(mk_ca_tau(ca_d(v,cell_nb)));
                
                /* diff eq spine */
                vsp1(v_dot,cell_nb) = (-Ilsp1-Idend_sp1-Isynapen1-Isynapin1-Icoupl1)+Ihyp;
                vsp2(v_dot,cell_nb) = (-Ilsp2-Idend_sp2-Isynapen2-Isynapin2-Icoupl2)+Ihyp;
                vsp3(v_dot,cell_nb) = (-Ilsp3-Idend_sp3-Isynapen3-Isynapin3-Icoupl3)+Ihyp;
                vsp4(v_dot,cell_nb) = (-Ilsp4-Idend_sp4-Isynapen4-Isynapin4-Icoupl4)+Ihyp;
            }
        }
    };
    
    struct io_model_with_noise
    {
        template< class State , class Deriv >
                void operator()( const State &x_ , Deriv &dxdt_ , double t ) const
        {
            typename boost::range_iterator< const State >::type v = boost::begin( x_ );
            typename boost::range_iterator< Deriv >::type v_dot = boost::begin( dxdt_ );
            
            double	Ical, Iap, Idend, Ils, Ikars, Isynapen_s, Isynapin_s;
            double	Icah, Ikcad, Isoma, Ild, Isynapen_d, Isynapin_d;
            double  Idend_sp1, Idend_sp2, Idend_sp3, Idend_sp4;
            double  Ispin1, Ispin2, Ispin3, Ispin4;
            double  Ilsp1, Ilsp2, Ilsp3, Ilsp4;
            double	Isynapen1, Isynapen2, Isynapen3, Isynapen4;
            double	Isynapin1, Isynapin2, Isynapin3, Isynapin4;
            double	Icoupl1, Icoupl2, Icoupl3, Icoupl4;
            
            int tidx = (int) t/0.5;
            
            for( size_t cell_nb = 0; cell_nb < NCELL; cell_nb++){
                switch(cell_nb){
                    
                    case 0:
                        Icoupl1 = (1+0.2*gc_table[0])*gc*(vsp1(v,0)-vsp3(v,1));
                        Icoupl2 = (1+0.2*gc_table[2])*gc*(vsp2(v,0)-vsp4(v,3));
                        Icoupl3 = (1+0.2*gc_table[12])*gc*(vsp3(v,0)-vsp1(v,2));
                        Icoupl4 = (1+0.2*gc_table[13])*gc*(vsp4(v,0)-vsp2(v,6));
                        break;
                        
                    case 1:
                        Icoupl1 = (1+0.2*gc_table[1])*gc*(vsp1(v,1)-vsp3(v,2));
                        Icoupl2 = (1+0.2*gc_table[3])*gc*(vsp2(v,1)-vsp4(v,4));
                        Icoupl3 = (1+0.2*gc_table[0])*gc*(vsp3(v,1)-vsp1(v,0));
                        Icoupl4 = (1+0.2*gc_table[17])*gc*(vsp4(v,1)-vsp2(v,7));
                        break;
                        
                    case 2:
                        Icoupl1 = (1+0.2*gc_table[12])*gc*(vsp1(v,2)-vsp3(v,0));
                        Icoupl2 = (1+0.2*gc_table[4])*gc*(vsp2(v,2)-vsp4(v,5));
                        Icoupl3 = (1+0.2*gc_table[1])*gc*(vsp3(v,2)-vsp1(v,1));
                        Icoupl4 = (1+0.2*gc_table[14])*gc*(vsp4(v,2)-vsp2(v,8));
                        break;
                        
                    case 3:
                        Icoupl1 = (1+0.2*gc_table[5])*gc*(vsp1(v,3)-vsp3(v,4));
                        Icoupl2 = (1+0.2*gc_table[7])*gc*(vsp2(v,3)-vsp4(v,6));
                        Icoupl3 = (1+0.2*gc_table[16])*gc*(vsp3(v,3)-vsp1(v,5));
                        Icoupl4 = (1+0.2*gc_table[2])*gc*(vsp4(v,3)-vsp2(v,0));
                        break;
                        
                    case 4:
                        Icoupl1 = (1+0.2*gc_table[6])*gc*(vsp1(v,4)-vsp3(v,5));
                        Icoupl2 = (1+0.2*gc_table[8])*gc*(vsp2(v,4)-vsp4(v,7));
                        Icoupl3 = (1+0.2*gc_table[5])*gc*(vsp3(v,4)-vsp1(v,3));
                        Icoupl4 = (1+0.2*gc_table[3])*gc*(vsp4(v,4)-vsp2(v,1));
                        break;
                        
                    case 5:
                        Icoupl1 = (1+0.2*gc_table[16])*gc*(vsp1(v,5)-vsp3(v,3));
                        Icoupl2 = (1+0.2*gc_table[9])*gc*(vsp2(v,5)-vsp4(v,8));
                        Icoupl3 = (1+0.2*gc_table[6])*gc*(vsp3(v,5)-vsp1(v,4));
                        Icoupl4 = (1+0.2*gc_table[4])*gc*(vsp4(v,5)-vsp2(v,2));
                        break;
                        
                    case 6:
                        Icoupl1 = (1+0.2*gc_table[10])*gc*(vsp1(v,6)-vsp3(v,7));
                        Icoupl2 = (1+0.2*gc_table[13])*gc*(vsp2(v,6)-vsp4(v,0));
                        Icoupl3 = (1+0.2*gc_table[15])*gc*(vsp3(v,6)-vsp1(v,8));
                        Icoupl4 = (1+0.2*gc_table[7])*gc*(vsp4(v,6)-vsp2(v,3));
                        break;
                        
                    case 7:
                        Icoupl1 = (1+0.2*gc_table[11])*gc*(vsp1(v,7)-vsp3(v,8));
                        Icoupl2 = (1+0.2*gc_table[17])*gc*(vsp2(v,7)-vsp4(v,1));
                        Icoupl3 = (1+0.2*gc_table[10])*gc*(vsp3(v,7)-vsp1(v,6));
                        Icoupl4 = (1+0.2*gc_table[8])*gc*(vsp4(v,7)-vsp2(v,4));
                        break;
                        
                    case 8:
                        Icoupl1 = (1+0.2*gc_table[15])*gc*(vsp1(v,8)-vsp3(v,6));
                        Icoupl2 = (1+0.2*gc_table[14])*gc*(vsp2(v,8)-vsp4(v,2));
                        Icoupl3 = (1+0.2*gc_table[11])*gc*(vsp3(v,8)-vsp1(v,7));
                        Icoupl4 = (1+0.2*gc_table[9])*gc*(vsp4(v,8)-vsp2(v,5));
                        break;
                }
                
                /* currents of somas */
                Ical = gcal*(1.02+0.05*gcal_table[cell_nb])*(scal(v,cell_nb)*scal(v,cell_nb)*scal(v,cell_nb))*hcal(v,cell_nb)*(vs(v,cell_nb)-vca);
                Iap = gna*(minf(vs(v,cell_nb))*minf(vs(v,cell_nb))*minf(vs(v,cell_nb)))*h(v,cell_nb)*(vs(v,cell_nb)-vna)+gk*(ns(v,cell_nb)*ns(v,cell_nb)*ns(v,cell_nb)*ns(v,cell_nb))*(vs(v,cell_nb)-vk);
                Idend = (ginter/p)*(vs(v,cell_nb)-vd(v,cell_nb));
                Ils = gl*(vs(v,cell_nb)-vl);
                Ikars = gk_ar*ar(v,cell_nb)*(vs(v,cell_nb)+43);
                
                Isynapen_s = (gsyn_en*excitatory_noise[cell_nb][0][tidx])*(vs(v,cell_nb)-vsynapex);
                Isynapin_s = (gsyn_in*inhibitory_noise[cell_nb][0][tidx])*(vs(v,cell_nb)-vsynapin);
                
                /* currents of dendrites */
                Icah = gcah*scah(v,cell_nb)*scah(v,cell_nb)*(vd(v,cell_nb)-vca);
                Ikcad = gk_cad*mk_cad(v,cell_nb)*(vd(v,cell_nb)-vk);
                Isoma = (ginter/(1.-p-pp))*(vd(v,cell_nb)-vs(v,cell_nb));
                Ild = gl*(vd(v,cell_nb)-vl);
                
                Isynapen_d = (gsyn_en*excitatory_noise[cell_nb][1][tidx])*(vd(v,cell_nb)-vsynapex);
                Isynapin_d = (gsyn_in*inhibitory_noise[cell_nb][1][tidx])*(vd(v,cell_nb)-vsynapin);
                
                /* currents of spines */
                Ispin1 = (ginter_sp/(1.-p-pp))*(vd(v,cell_nb)-vsp1(v,cell_nb));
                Ispin2 = (ginter_sp/(1.-p-pp))*(vd(v,cell_nb)-vsp2(v,cell_nb));
                Ispin3 = (ginter_sp/(1.-p-pp))*(vd(v,cell_nb)-vsp3(v,cell_nb));
                Ispin4 = (ginter_sp/(1.-p-pp))*(vd(v,cell_nb)-vsp4(v,cell_nb));
                
                Isynapen1 = (gsyn_en*excitatory_noise[cell_nb][2][tidx])*(vsp1(v,cell_nb)-vsynapen_sp);
                Isynapen2 = (gsyn_en*excitatory_noise[cell_nb][2][tidx])*(vsp2(v,cell_nb)-vsynapen_sp);
                Isynapen3 = (gsyn_en*excitatory_noise[cell_nb][2][tidx])*(vsp3(v,cell_nb)-vsynapen_sp);
                Isynapen4 = (gsyn_en*excitatory_noise[cell_nb][2][tidx])*(vsp4(v,cell_nb)-vsynapen_sp);
                Isynapin1 = (gsyn_in*inhibitory_noise[cell_nb][2][tidx])*(vsp1(v,cell_nb)-vsynapin_sp);
                Isynapin2 = (gsyn_in*inhibitory_noise[cell_nb][2][tidx])*(vsp2(v,cell_nb)-vsynapin_sp);
                Isynapin3 = (gsyn_in*inhibitory_noise[cell_nb][2][tidx])*(vsp3(v,cell_nb)-vsynapin_sp);
                Isynapin4 = (gsyn_in*inhibitory_noise[cell_nb][2][tidx])*(vsp4(v,cell_nb)-vsynapin_sp);
                
                Idend_sp1 = (ginter_sp/(0.25*pp))*(vsp1(v,cell_nb)-vd(v,cell_nb));
                Idend_sp2 = (ginter_sp/(0.25*pp))*(vsp2(v,cell_nb)-vd(v,cell_nb));
                Idend_sp3 = (ginter_sp/(0.25*pp))*(vsp3(v,cell_nb)-vd(v,cell_nb));
                Idend_sp4 = (ginter_sp/(0.25*pp))*(vsp4(v,cell_nb)-vd(v,cell_nb));
                
                Ilsp1 = gl*(vsp1(v,cell_nb)-vl);
                Ilsp2 = gl*(vsp2(v,cell_nb)-vl);
                Ilsp3 = gl*(vsp3(v,cell_nb)-vl);
                Ilsp4 = gl*(vsp4(v,cell_nb)-vl);
                
                /* diff eq soma */
                vs(v_dot,cell_nb) =(-Ical-Iap-Ils-Ikars-Idend-Isynapin_s-Isynapen_s);
                hcal(v_dot,cell_nb) = (hcalinf(vs(v,cell_nb))-hcal(v,cell_nb))/(hcaltau(vs(v,cell_nb)));
                scal(v_dot,cell_nb) = (scalinf(vs(v,cell_nb))-scal(v,cell_nb))/(scaltau(vs(v,cell_nb)));
                ns(v_dot,cell_nb) = (ninf(vs(v,cell_nb))-ns(v,cell_nb))/(ntau(vs(v,cell_nb)));
                h(v_dot,cell_nb) = (hinf(vs(v,cell_nb))-h(v,cell_nb))/(htau(vs(v,cell_nb)));
                ar(v_dot,cell_nb) = (arinf(vs(v,cell_nb))-ar(v,cell_nb))/(artau(vs(v,cell_nb)));
                
                /* diff eq dendrite */
                vd(v_dot,cell_nb) = (-Icah-Ikcad-Ild-Isoma-Isynapen_d-Isynapin_d -Ispin1-Ispin2-Ispin3-Ispin4);
                scah(v_dot,cell_nb) = (infscah(vd(v,cell_nb))-scah(v,cell_nb))/(tauscah(vd(v,cell_nb)));
                ca_d(v_dot,cell_nb) = -1.01*Icah-0.02*ca_d(v,cell_nb);
                mk_cad(v_dot,cell_nb) = (mk_ca_inf(ca_d(v,cell_nb))-mk_cad(v,cell_nb))/(mk_ca_tau(ca_d(v,cell_nb)));
                
                /* diff eq spine */
                vsp1(v_dot,cell_nb) = (-Ilsp1-Idend_sp1-Isynapen1-Isynapin1-Icoupl1);
                vsp2(v_dot,cell_nb) = (-Ilsp2-Idend_sp2-Isynapen2-Isynapin2-Icoupl2);
                vsp3(v_dot,cell_nb) = (-Ilsp3-Idend_sp3-Isynapen3-Isynapin3-Icoupl3);
                vsp4(v_dot,cell_nb) = (-Ilsp4-Idend_sp4-Isynapen4-Isynapin4-Icoupl4);
            }
        }
    };
    
    
    void io_model_with_lyap( const state_type &v , state_type &dvdt , double t )
    {
        double JIcal, JIap, JIdend, JIls, JIkars, JIsynapen_s, JIsynapin_s;
        double JIcah, JIkcad, JIsoma, JIld, JIsynapen_d, JIsynapin_d, JIspin;
        double JIlsp, JIdend_sp, JIsynapen, JIsynapin, JIcoupl1, JIcoupl2, JIcoupl3, JIcoupl4;
        
        io_model()( v , dvdt , t );
        
        for( size_t l=0 ; l<num_of_lyap ; ++l )
        {
            const double *pert = v.begin() + num_dim + l * num_dim;
            double *dpert = dvdt.begin() + num_dim + l * num_dim;
            
            for ( size_t cell_nb=0 ; cell_nb<NCELL; cell_nb++ ) {
                /* diff eq soma */
                JIcal = gcal*(1.02+0.05*gcal_table[cell_nb])*(scal(v,cell_nb)*scal(v,cell_nb)*scal(v,cell_nb))*hcal(v,cell_nb);
                JIap = gna*minf(vs(v,cell_nb))*minf(vs(v,cell_nb))*(3*vs(v,cell_nb)*minfdot(vs(v,cell_nb))+minf(vs(v,cell_nb)))*h(v,cell_nb) + gk*(ns(v,cell_nb)*ns(v,cell_nb)*ns(v,cell_nb)*ns(v,cell_nb));
                JIdend = (ginter/p);
                JIls = gl;
                JIkars = gk_ar*ar(v,cell_nb);
                JIsynapen_s = (gsyn_en*synaptic_input);
                JIsynapin_s = (gsyn_in*synaptic_input);
                
                vs(dpert,cell_nb) = (-JIcal-JIap-JIdend-JIls-JIkars-JIsynapen_s-JIsynapin_s)*vs(pert,cell_nb) +
                        (-gcal*(1.02+0.05*gcal_table[cell_nb])*(scal(v,cell_nb)*scal(v,cell_nb)*scal(v,cell_nb))*(vs(v,cell_nb)-vca))*hcal(pert,cell_nb) +
                        (-3*gcal*(1.02+0.05*gcal_table[cell_nb])*(scal(v,cell_nb)*scal(v,cell_nb))*hcal(v,cell_nb)*(vs(v,cell_nb)-vca))*scal(pert,cell_nb) +
                        (-4*gk*(ns(v,cell_nb)*ns(v,cell_nb)*ns(v,cell_nb))*(vs(v,cell_nb)-vk))*ns(pert,cell_nb) +
                        (-gna*(minf(vs(v,cell_nb))*minf(vs(v,cell_nb))*minf(vs(v,cell_nb)))*(vs(v,cell_nb)-vna))*h(pert,cell_nb) +
                        (-gk_ar*(vs(v,cell_nb)+43))*ar(pert,cell_nb) +
                        (ginter/p)*vd(pert,cell_nb);
                
                hcal(dpert,cell_nb) = ((hcaltau(vs(v,cell_nb))*hcalinfdot(vs(v,cell_nb)) - (hcalinf(vs(v,cell_nb))-hcal(v,cell_nb))*hcaltaudot(vs(v,cell_nb)))/(hcaltau(vs(v,cell_nb))*hcaltau(vs(v,cell_nb))))*vs(pert,cell_nb) +
                        (-1/hcaltau(vs(v,cell_nb)))*hcal(pert,cell_nb);
                
                scal(dpert,cell_nb) = (scalinfdot(vs(v,cell_nb))/scaltau(vs(v,cell_nb)))*vs(pert,cell_nb) +
                        (-1/scaltau(vs(v,cell_nb)))*scal(pert,cell_nb);
                
                ns(dpert,cell_nb) = ((ntau(vs(v,cell_nb))*ninfdot(vs(v,cell_nb))-(ninf(vs(v,cell_nb))-ns(v,cell_nb))*ntaudot(vs(v,cell_nb)))/(ntau(vs(v,cell_nb))*ntau(vs(v,cell_nb))))*vs(pert,cell_nb) +
                        (-1/ntau(vs(v,cell_nb)))*ns(pert,cell_nb);
                
                h(dpert,cell_nb) = ((htau(vs(v,cell_nb))*hinfdot(vs(v,cell_nb))-(hinf(vs(v,cell_nb))-h(v,cell_nb))*htaudot(vs(v,cell_nb)))/(htau(vs(v,cell_nb))*htau(vs(v,cell_nb))))*vs(pert,cell_nb) +
                        (-1/htau(vs(v,cell_nb)))*h(pert,cell_nb);
                
                ar(dpert,cell_nb) = ((artau(vs(v,cell_nb))*arinfdot(vs(v,cell_nb))-(arinf(vs(v,cell_nb))-ar(v,cell_nb))*artaudot(vs(v,cell_nb)))/(artau(vs(v,cell_nb))*artau(vs(v,cell_nb))))*vs(pert,cell_nb) +
                        (-1/artau(vs(v,cell_nb)))*ar(pert,cell_nb);
                
                /* diff eq dendrite */
                JIcah = gcah*scah(v,cell_nb)*scah(v,cell_nb);
                JIkcad = gk_cad*mk_cad(v,cell_nb);
                JIsoma = (ginter/(1.-p-pp));
                JIld = gl;
                JIsynapen_d = (gsyn_en*synaptic_input);
                JIsynapin_d = (gsyn_in*synaptic_input*dend_ratio);
                JIspin = (ginter_sp/(1.-p-pp));
                
                vd(dpert,cell_nb) = (ginter/(1.-p-pp))*vs(pert,cell_nb) +
                        (-JIcah-JIkcad-JIld-JIsoma-JIsynapen_d-JIsynapin_d-4*JIspin)*vd(pert,cell_nb) +
                        (-2*gcah*scah(v,cell_nb)*(vd(v,cell_nb)-vca))*scah(pert,cell_nb) +
                        (-gk_cad*(vd(v,cell_nb)-vk))*mk_cad(pert,cell_nb) +
                        (ginter_sp/(1.-p-pp))*vsp1(pert,cell_nb) +
                        (ginter_sp/(1.-p-pp))*vsp2(pert,cell_nb) +
                        (ginter_sp/(1.-p-pp))*vsp3(pert,cell_nb) +
                        (ginter_sp/(1.-p-pp))*vsp4(pert,cell_nb);
                
                scah(dpert,cell_nb) = ((tauscah(vd(v,cell_nb))*infscahdot(vd(v,cell_nb)) - (infscah(vd(v,cell_nb))-scah(v,cell_nb))*tauscahdot(vd(v,cell_nb)))/(tauscah(vd(v,cell_nb))*tauscah(vd(v,cell_nb))))*vd(pert,cell_nb) +
                        (-1/tauscah(vd(v,cell_nb)))*scah(pert,cell_nb);
                
                ca_d(dpert,cell_nb) = -JIcah*vd(pert,cell_nb) +
                        (-2*gcah*scah(v,cell_nb)*(vd(v,cell_nb)-vca))*scah(pert,cell_nb) +
                        -0.02*ca_d(pert,cell_nb);
                
                mk_cad(dpert,cell_nb) = ((mk_ca_tau(ca_d(v,cell_nb))*mk_ca_infdot(ca_d(v,cell_nb)) - (mk_ca_inf(ca_d(v,cell_nb))-mk_cad(v,cell_nb))*mk_ca_taudot(ca_d(v,cell_nb)))/(mk_ca_tau(ca_d(v,cell_nb))*mk_ca_tau(ca_d(v,cell_nb))))*ca_d(pert,cell_nb) +
                        (-1/mk_ca_tau(ca_d(v,cell_nb)))*mk_cad(pert,cell_nb);
                
                /* diff eq spine */
                JIlsp = gl;
                JIdend_sp = (ginter_sp/(0.25*pp));
                JIsynapen = (gsyn_en*synaptic_input);
                JIsynapin = (gsyn_in*synaptic_input);
                
                switch(cell_nb){
                    
                    case 0:
                        JIcoupl1 = (1+0.2*gc_table[0])*gc;
                        JIcoupl2 = (1+0.2*gc_table[2])*gc;
                        JIcoupl3 = (1+0.2*gc_table[12])*gc;
                        JIcoupl4 = (1+0.2*gc_table[13])*gc;
                        
                        vsp1(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl1)*vsp1(pert,cell_nb) +
                                JIcoupl1*vsp3(pert,1);
                        
                        vsp2(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl2)*vsp2(pert,cell_nb) +
                                JIcoupl2*vsp4(pert,3);
                        
                        vsp3(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl3)*vsp3(pert,cell_nb) +
                                JIcoupl3*vsp1(pert,2);
                        
                        vsp4(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl4)*vsp4(pert,cell_nb) +
                                JIcoupl4*vsp2(pert,6);
                        break;
                        
                    case 1:
                        JIcoupl1 = (1+0.2*gc_table[1])*gc;
                        JIcoupl2 = (1+0.2*gc_table[3])*gc;
                        JIcoupl3 = (1+0.2*gc_table[0])*gc;
                        JIcoupl4 = (1+0.2*gc_table[17])*gc;
                        
                        vsp1(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl1)*vsp1(pert,cell_nb) +
                                JIcoupl1*vsp3(pert,2);
                        
                        vsp2(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl2)*vsp2(pert,cell_nb) +
                                JIcoupl2*vsp4(pert,4);
                        
                        vsp3(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl3)*vsp3(pert,cell_nb) +
                                JIcoupl3*vsp1(pert,0);
                        
                        vsp4(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl4)*vsp4(pert,cell_nb) +
                                JIcoupl4*vsp2(pert,7);
                        break;
                        
                    case 2:
                        JIcoupl1 = (1+0.2*gc_table[12])*gc;
                        JIcoupl2 = (1+0.2*gc_table[4])*gc;
                        JIcoupl3 = (1+0.2*gc_table[1])*gc;
                        JIcoupl4 = (1+0.2*gc_table[14])*gc;
                        
                        vsp1(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl1)*vsp1(pert,cell_nb) +
                                JIcoupl1*vsp3(pert,0);
                        
                        vsp2(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl2)*vsp2(pert,cell_nb) +
                                JIcoupl2*vsp4(pert,5);
                        
                        vsp3(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl3)*vsp3(pert,cell_nb) +
                                JIcoupl3*vsp1(pert,1);
                        
                        vsp4(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl4)*vsp4(pert,cell_nb) +
                                JIcoupl4*vsp2(pert,8);
                        break;
                        
                    case 3:
                        JIcoupl1 = (1+0.2*gc_table[5])*gc;
                        JIcoupl2 = (1+0.2*gc_table[7])*gc;
                        JIcoupl3 = (1+0.2*gc_table[16])*gc;
                        JIcoupl4 = (1+0.2*gc_table[2])*gc;
                        
                        vsp1(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl1)*vsp1(pert,cell_nb) +
                                JIcoupl1*vsp3(pert,4);
                        
                        vsp2(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl2)*vsp2(pert,cell_nb) +
                                JIcoupl2*vsp4(pert,6);
                        
                        vsp3(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl3)*vsp3(pert,cell_nb) +
                                JIcoupl3*vsp1(pert,5);
                        
                        vsp4(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl4)*vsp4(pert,cell_nb) +
                                JIcoupl4*vsp2(pert,0);
                        break;
                        
                    case 4:
                        JIcoupl1 = (1+0.2*gc_table[6])*gc;
                        JIcoupl2 = (1+0.2*gc_table[8])*gc;
                        JIcoupl3 = (1+0.2*gc_table[5])*gc;
                        JIcoupl4 = (1+0.2*gc_table[3])*gc;
                        
                        vsp1(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl1)*vsp1(pert,cell_nb) +
                                JIcoupl1*vsp3(pert,5);
                        
                        vsp2(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl2)*vsp2(pert,cell_nb) +
                                JIcoupl2*vsp4(pert,7);
                        
                        vsp3(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl3)*vsp3(pert,cell_nb) +
                                JIcoupl3*vsp1(pert,3);
                        
                        vsp4(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl4)*vsp4(pert,cell_nb) +
                                JIcoupl4*vsp2(pert,1);
                        break;
                        
                    case 5:
                        JIcoupl1 = (1+0.2*gc_table[16])*gc;
                        JIcoupl2 = (1+0.2*gc_table[9])*gc;
                        JIcoupl3 = (1+0.2*gc_table[6])*gc;
                        JIcoupl4 = (1+0.2*gc_table[4])*gc;
                        
                        vsp1(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl1)*vsp1(pert,cell_nb) +
                                JIcoupl1*vsp3(pert,3);
                        
                        vsp2(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl2)*vsp2(pert,cell_nb) +
                                JIcoupl2*vsp4(pert,8);
                        
                        vsp3(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl3)*vsp3(pert,cell_nb) +
                                JIcoupl3*vsp1(pert,4);
                        
                        vsp4(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl4)*vsp4(pert,cell_nb) +
                                JIcoupl4*vsp2(pert,2);
                        break;
                        
                    case 6:
                        JIcoupl1 = (1+0.2*gc_table[10])*gc;
                        JIcoupl2 = (1+0.2*gc_table[13])*gc;
                        JIcoupl3 = (1+0.2*gc_table[15])*gc;
                        JIcoupl4 = (1+0.2*gc_table[7])*gc;
                        
                        vsp1(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl1)*vsp1(pert,cell_nb) +
                                JIcoupl1*vsp3(pert,7);
                        
                        vsp2(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl2)*vsp2(pert,cell_nb) +
                                JIcoupl2*vsp4(pert,0);
                        
                        vsp3(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl3)*vsp3(pert,cell_nb) +
                                JIcoupl3*vsp1(pert,8);
                        
                        vsp4(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl4)*vsp4(pert,cell_nb) +
                                JIcoupl4*vsp2(pert,3);
                        break;
                        
                    case 7:
                        JIcoupl1 = (1+0.2*gc_table[11])*gc;
                        JIcoupl2 = (1+0.2*gc_table[17])*gc;
                        JIcoupl3 = (1+0.2*gc_table[10])*gc;
                        JIcoupl4 = (1+0.2*gc_table[8])*gc;
                        
                        vsp1(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl1)*vsp1(pert,cell_nb) +
                                JIcoupl1*vsp3(pert,8);
                        
                        vsp2(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl2)*vsp2(pert,cell_nb) +
                                JIcoupl2*vsp4(pert,1);
                        
                        vsp3(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl3)*vsp3(pert,cell_nb) +
                                JIcoupl3*vsp1(pert,6);
                        
                        vsp4(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl4)*vsp4(pert,cell_nb) +
                                JIcoupl4*vsp2(pert,4);
                        break;
                        
                    case 8:
                        JIcoupl1 = (1+0.2*gc_table[15])*gc;
                        JIcoupl2 = (1+0.2*gc_table[14])*gc;
                        JIcoupl3 = (1+0.2*gc_table[11])*gc;
                        JIcoupl4 = (1+0.2*gc_table[9])*gc;
                        
                        vsp1(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl1)*vsp1(pert,cell_nb) +
                                JIcoupl1*vsp3(pert,6);
                        
                        vsp2(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl2)*vsp2(pert,cell_nb) +
                                JIcoupl2*vsp4(pert,2);
                        
                        vsp3(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl3)*vsp3(pert,cell_nb) +
                                JIcoupl3*vsp1(pert,7);
                        
                        vsp4(dpert,cell_nb) = (JIdend_sp)*vd(pert,cell_nb) +
                                (-JIlsp-JIdend_sp-JIsynapen-JIsynapin-JIcoupl4)*vsp4(pert,cell_nb) +
                                JIcoupl4*vsp2(pert,5);
                        break;
                }
            }
        }
    };
    
    
#endif
