//////////////////////////////////////////////////////////////////////////
////////////////              nanda.cxx         /////////////////////
//  Optimal control of treatment in a mathematical model of chronic myelogenous leukemia
//  Seema Nanda, Helen Moore, Suzanne Lenhart, Mathematical Biosciences 210 (2007) 143–156

#include "psopt.h"

/////////////////////////////  Strategy:        //////////////////////////
// start with bioreactor example and then use bryson's max range template (user)//
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   adouble x1tf = final_states[0];
   adouble x3tf = final_states[2];

   return 0.1*x3tf-100000*x1tf; // B3*C(tf)-B4*Tn(Tf)
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
     adouble u1 = controls[ 0 ];
     adouble u2 = controls[ 1 ];
     adouble x3 = states[ 2 ];

     return   (x3 + 500.0*u1*u1 + 250.0*u2*u2);

}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

    double kn = 0.066;
    double eta = 140.0;
    double Cmax = 160000.0;
    double dn = 0.35;
    double de = 0.40;
    double dc = 0.012;

    adouble u1 = controls[ 0 ];
    adouble u2 = controls[ 1 ];

    adouble Tn = states[ 0 ];
    adouble Te = states[ 1 ];
    adouble C  = states[ 2 ];

    derivatives[ 0 ] =  0.29 - dn*u2*Tn -kn*Tn*C/(C+eta);
    derivatives[ 1 ] =  0.39*kn*Tn*C/(C+eta) + 0.65*Te*C/(C+eta) - de*u2*Te - 0.079*C*Te;
    derivatives[ 2 ] =  (1-u1)*0.011*C*log(Cmax/C) - dc*u2*C - 0.058*C*Te;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{

   int i;

   for(i=0;i< 3; i++) {
          e[i] = initial_states[i];
   }



}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
// Single phase problem

}



////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Nanda et al 2007 CML";
    problem.outfilename                 = "nanda.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   					 = 1;
    problem.nlinkages            	 = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 3;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 3;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes         << 100;


    psopt_level2_setup(problem, algorithm);



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    int i;


    problem.phases(1).bounds.lower.states      << 0.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.states      << 2000.0, 100.0, 100000.0;


    problem.phases(1).bounds.lower.controls(0) = 0.0;
    problem.phases(1).bounds.lower.controls(1) = 1.0;



    problem.phases(1).bounds.upper.controls(0) = 0.9;
    problem.phases(1).bounds.upper.controls(1) = 2.5;


    double Tn0 = 1510.0;
    double Te0 = 10.0;
    double C0 = 10000.0;

    problem.phases(1).bounds.lower.events(0) = Tn0;
    problem.phases(1).bounds.lower.events(1) = Te0;
    problem.phases(1).bounds.lower.events(2) = C0;


    problem.phases(1).bounds.upper.events(0) = Tn0;
    problem.phases(1).bounds.upper.events(1) = Te0;
    problem.phases(1).bounds.upper.events(2) = C0;

//  chunk above equals this
//  problem.phases(1).bounds.lower.events     << 1510.0, 10.0, 10000.0;
//  problem.phases(1).bounds.upper.events     << 1510.0, 10.0, 10000.0;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 250.0;
    problem.phases(1).bounds.upper.EndTime      = 250.0;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 		= &integrand_cost;
    problem.endpoint_cost 	    	= &endpoint_cost;
    problem.dae             		= &dae;
    problem.events 					= &events;
    problem.linkages					= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			= problem.phases(1).nodes(0);
    int ncontrols          = problem.phases(1).ncontrols;
    int nstates            = problem.phases(1).nstates;

    MatrixXd x_guess    =  zeros(nstates,nnodes);  // follow bryson max range example here
    x_guess.row(0)  = Tn0*ones(1,nnodes);
    x_guess.row(1)  = Te0*ones(1,nnodes);
    x_guess.row(2)  = C0*ones(1,nnodes);

//    MatrixXd state_guess    =  zeros(nstates,nnodes);  // as this is getting cryptic
//    state_guess             <<   Tn0*ones(1,nnodes),Te0*ones(1,nnodes),C0*ones(1,nnodes);

//    MatrixXd control_guess  =  zeros(ncontrols,nnodes);
//    MatrixXd time_guess     =  ;

    problem.phases(1).guess.controls       = ones(ncontrols,nnodes);
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = linspace(0,250.0,nnodes);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////




    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-4;
//    algorithm.nlp_tolerance               = 1.e-5;  //bioreactor
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
//    algorithm.mesh_refinement             = "automatic";
//    algorithm.collocation_method          = "trapezoidal";
    algorithm.collocation_method          = "Hermite-Simpson";  //bioreactor
//    algorithm.defect_scaling              = "jacobian-based";  //commented in both
//    algorithm.diff_matrix                 = "central-differences"; //comment in bio, not in max
    algorithm.ode_tolerance               = 1.e-6;  // in max but not in bio



////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    MatrixXd x, u, t;
    MatrixXd lambda, H;


    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);
    lambda = solution.get_dual_costates_in_phase(1);
    H      = solution.get_dual_hamiltonian_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");
    Save(lambda,"lambda.dat");
    Save(H,"H.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////


    plot(t,u,problem.name+": control","time (s)", "control 1", "u1 u2");

    plot(t,x,problem.name+": state","time (s)", "state", "Tn Te C");


    plot(t,u,problem.name+": control","time (s)", "control 1", "u1 u2",
	      "pdf", "nanda_controls.pdf");

    plot(t,x,problem.name+": states","time (s)", "states", "Tn Te C",
	      "pdf", "nanda_states.pdf");

    plot(t,u,problem.name+": control","time (s)", "control 1", "u1 u2",
	      "png", "nanda_controls.png");


    plot(t,x,problem.name+": states","time (s)", "states", "Tn Te C",
	      "png", "nanda_states.png");



}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
