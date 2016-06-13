#include "psopt.h"
#include <iostream>

// Endpoint cost, returns nothing for this problem
adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
	return 0.0;
}

// Lagrange cost: return the integral of 1/2u^T*u
adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
	adouble ux = controls[ CINDEX(1) ];
	adouble uy = controls[ CINDEX(2) ];
	adouble uz = controls[ CINDEX(3) ];

	adouble l2Energy = 0.5*ux*ux + 0.5*uy*uy + 0.5*uz*uz;

	return l2Energy;
}

// Differential equation system
void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
	adouble xdot, ydot, zdot, xddot, yddot, zddot;
	
	// Define some constants
	double mu 	= 3.986e14;
	double a 	= 6678e3;
	double n 	= sqrt(mu/(a*a*a));
	double mass = 100.0;
	
	adouble x  = states[ CINDEX(1) ];
	adouble y  = states[ CINDEX(2) ];
	adouble z  = states[ CINDEX(3) ];
	adouble vx = states[ CINDEX(4) ];
	adouble vy = states[ CINDEX(5) ];
	adouble vz = states[ CINDEX(6) ];
	
	adouble ux = controls[ CINDEX(1) ];
	adouble uy = controls[ CINDEX(2) ];
	adouble uz = controls[ CINDEX(3) ];
	
	// Hill-clohessy-Wiltshire equations
	xdot  = vx;
	ydot  = vy;
	zdot  = vz;
	xddot = 3*n*n*x + 2*n*vy + ux/mass;
	yddot = -2*n*vx          + uy/mass;
	zddot = -n*n*z           + uz/mass;
	
	derivatives[ CINDEX(1) ] = xdot;
	derivatives[ CINDEX(2) ] = ydot;
	derivatives[ CINDEX(3) ] = zdot;
	derivatives[ CINDEX(4) ] = xddot;
	derivatives[ CINDEX(5) ] = yddot;
	derivatives[ CINDEX(6) ] = zddot;
}

// Boundary conditions
void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
	adouble x0  = initial_states[ CINDEX(1) ];
	adouble y0  = initial_states[ CINDEX(2) ];
	adouble z0  = initial_states[ CINDEX(3) ];
	adouble vx0 = initial_states[ CINDEX(4) ];
	adouble vy0 = initial_states[ CINDEX(5) ];
	adouble vz0 = initial_states[ CINDEX(6) ];
	
	adouble xf  = final_states[ CINDEX(1) ];
	adouble yf  = final_states[ CINDEX(2) ];
	adouble zf  = final_states[ CINDEX(3) ];
	adouble vxf = final_states[ CINDEX(4) ];
	adouble vyf = final_states[ CINDEX(5) ];
	adouble vzf = final_states[ CINDEX(6) ];
	
	e[ CINDEX(1) ] = x0;
	e[ CINDEX(2) ] = y0;
	e[ CINDEX(3) ] = z0;
	e[ CINDEX(4) ] =vx0;
	e[ CINDEX(5) ] =vy0;
	e[ CINDEX(6) ] =vz0;
	
	e[ CINDEX(7) ] = xf;
	e[ CINDEX(8) ] = yf;
	e[ CINDEX(9) ] = zf;
	e[ CINDEX(10)] =vxf;
	e[ CINDEX(11)] =vyf;
	e[ CINDEX(12)] =vzf;
}

// Define linkages, none for this problem
void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
	//
}

// Main function
int main(void)
{
	Alg algorithm;
	Sol solution;
	Prob problem;
	
	problem.name        = "Hill-Clohessy-Wiltshire Minimum Energy";
	problem.outfilename = "hcw_min_energy.txt";
	
	problem.nphases 	= 1;
	problem.nlinkages 	= 0;
	psopt_level1_setup(problem);
	
	problem.phases(1).nstates   = 6;
	problem.phases(1).ncontrols = 3;
	problem.phases(1).nevents   = 12;
	problem.phases(1).npath     = 0;
	problem.phases(1).nodes 	= 80;
	
	psopt_level2_setup(problem, algorithm);
	
	DMatrix x, u, t;
	DMatrix lambda, H;
	
	// Define constants
	double mu 	= 3.986e14; 		// m^3/s^2
	double a 	= 6678e3;   		// m
	double n 	= sqrt(mu/(a*a*a));	// rad/s
	
	// Initial and final times
	double t0 	= 0.0; 				// s
	double tf 	= 800.0;			// s
	
	// Initial conditions (m, m/s)
	double x0 	= 0.0;				// m
	double y0 	= -100.0;			// m
	double z0 	= 0.0;				// m
	double vx0 	= 0.0;				// m/s
	double vy0 	= 0.0;				// m/s
	double vz0 	= 0.0;				// m/s
	
	// Terminal conditions
	double xf 	= 0.0;				// m
	double yf 	= -2.0;				// m
	double zf 	= 0.0;				// m
	double vxf 	= 0.0;				// m/s
	double vyf 	= 0.0;				// m/s
	double vzf 	= 0.0;				// m/s
	
	// Upper and lower bounds on events (boundary conditions)
	problem.phases(1).bounds.lower.events(1) 	= x0;
	problem.phases(1).bounds.lower.events(2) 	= y0;
	problem.phases(1).bounds.lower.events(3) 	= z0;
	problem.phases(1).bounds.lower.events(4) 	= vx0;
	problem.phases(1).bounds.lower.events(5) 	= vy0;
	problem.phases(1).bounds.lower.events(6) 	= vz0;
	problem.phases(1).bounds.lower.events(7) 	= xf;
	problem.phases(1).bounds.lower.events(8) 	= yf;
	problem.phases(1).bounds.lower.events(9) 	= zf;
	problem.phases(1).bounds.lower.events(10) 	= vxf;
	problem.phases(1).bounds.lower.events(11) 	= vyf;
	problem.phases(1).bounds.lower.events(12) 	= vzf;
	
	problem.phases(1).bounds.upper.events(1) 	= x0;
	problem.phases(1).bounds.upper.events(2) 	= y0;
	problem.phases(1).bounds.upper.events(3) 	= z0;
	problem.phases(1).bounds.upper.events(4)	= vx0;
	problem.phases(1).bounds.upper.events(5) 	= vy0;
	problem.phases(1).bounds.upper.events(6) 	= vz0;
	problem.phases(1).bounds.upper.events(7) 	= xf;
	problem.phases(1).bounds.upper.events(8) 	= yf;
	problem.phases(1).bounds.upper.events(9) 	= zf;
	problem.phases(1).bounds.upper.events(10) 	= vxf;
	problem.phases(1).bounds.upper.events(11) 	= vyf;
	problem.phases(1).bounds.upper.events(12) 	= vzf;
	
	// Upper and lower bounds on initial and final times
	problem.phases(1).bounds.lower.StartTime = t0;
	problem.phases(1).bounds.upper.StartTime = t0;
	
	problem.phases(1).bounds.lower.EndTime = tf;
	problem.phases(1).bounds.upper.EndTime = tf;
	
	// Function registers
	problem.integrand_cost 	= &integrand_cost;
	problem.endpoint_cost 	= &endpoint_cost;
	problem.dae 			= &dae;
	problem.events 			= &events;
	problem.linkages 		= &linkages;
	
	// Initial guess
	int nnodes 		= problem.phases(1).nodes(1);
	int ncontrols 	= problem.phases(1).ncontrols;
	int nstates 	= problem.phases(1).nstates;
	
	DMatrix xguess = zeros(nstates,nnodes);
	
	xguess(1,colon()) = x0*ones(1,nnodes);
	xguess(2,colon()) = y0*ones(1,nnodes);
	xguess(3,colon()) = z0*ones(1,nnodes);
	xguess(4,colon()) = vx0*ones(1,nnodes);
	xguess(5,colon()) = vy0*ones(1,nnodes);
	xguess(6,colon()) = vz0*ones(1,nnodes);

	problem.phases(1).guess.controls 	= zeros(ncontrols,nnodes);
	problem.phases(1).guess.states 		= xguess;
	problem.phases(1).guess.time 		= linspace(t0,tf,nnodes);
	
	// Algorithm options
	algorithm.nlp_iter_max 	= 1000;
	algorithm.nlp_tolerance = 1.e-6;
	algorithm.nlp_method 	= "IPOPT";
	algorithm.scaling 		= "automatic";
	algorithm.derivatives 	= "automatic";
	algorithm.nlp_tolerance = 1.e-6;
	//algorithm.collocation_method = "Hermite-Simpson";
	
	// Solve problem
	psopt(solution,problem,algorithm);
	
	// Post-process solution
	DMatrix muE;
	
	x 		= solution.get_states_in_phase(1);
	u 		= solution.get_controls_in_phase(1);
	t 		= solution.get_time_in_phase(1);
	lambda 	= solution.get_dual_costates_in_phase(1);
	H 		= solution.get_dual_hamiltonian_in_phase(1);
	muE 	= solution.get_dual_events_in_phase(1);

	x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");
    lambda.Save("p.dat");
    H.Save("H.dat");
    
    plot(t,x,problem.name+": states", "time (s)", "states","x y z vx vy vz");
    plot(t,u,problem.name+": controls","time (s)", "controls","ux uy uz");
    plot(t,H,problem.name+": Hamiltonian","time (s)","H","H");
}





































