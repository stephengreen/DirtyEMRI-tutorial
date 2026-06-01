#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_ellint.h>
#include <algorithm>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <complex>
#include <cmath>

#include "Interpolant.h"
#include "global.h"
#include "Utility.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <iomanip>      // std::setprecision
#include <cstring>

#include "dIdt8H_5PNe10.h"
#include "ode.hh"
#include "KerrEquatorial.h"

#define pn5_Y
#define pn5_citation1 Pn5_citation
__deriv__
void pn5(double* pdot, double* edot, double* Ydot,
                  double* Omega_phi, double* Omega_theta, double* Omega_r,
                  double epsilon, double a, double p, double e, double Y, double* additional_args)
{
    // evaluate ODEs

    // the frequency variables are pointers!
    double x = Y_to_xI(a, p, e, Y);
    KerrGeoCoordinateFrequencies(Omega_phi, Omega_theta, Omega_r, a, p, e, x);

	int Nv = 10;
    int ne = 10;
    *pdot = epsilon * dpdt8H_5PNe10 (a, p, e, Y, Nv, ne);

    // needs adjustment for validity
    Nv = 10;
    ne = 8;
	*edot = epsilon * dedt8H_5PNe10 (a, p, e, Y, Nv, ne);

    Nv = 7;
    ne = 10;
    *Ydot = epsilon * dYdt8H_5PNe10 (a, p, e, Y, Nv, ne);

}


// Initialize flux data for inspiral calculations
void load_and_interpolate_flux_data(struct interp_params *interps, const std::string& few_dir){

	// Load and interpolate the flux data
    std::string fp = "few/files/FluxNewMinusPNScaled_fixed_y_order.dat";
    fp = few_dir + fp;
	ifstream Flux_file(fp);

    if (Flux_file.fail())
    {
        throw std::runtime_error("The file FluxNewMinusPNScaled_fixed_y_order.dat did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
    }

	// Load the flux data into arrays
	string Flux_string;
	vector<double> ys, es, Edots, Ldots;
	double y, e, Edot, Ldot;
	while(getline(Flux_file, Flux_string)){

		stringstream Flux_ss(Flux_string);

		Flux_ss >> y >> e >> Edot >> Ldot;

		ys.push_back(y);
		es.push_back(e);
		Edots.push_back(Edot);
		Ldots.push_back(Ldot);
	}

	// Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)
	sort( ys.begin(), ys.end() );
	ys.erase( unique( ys.begin(), ys.end() ), ys.end() );

	sort( es.begin(), es.end() );
	es.erase( unique( es.begin(), es.end() ), es.end() );

	Interpolant *Edot_interp = new Interpolant(ys, es, Edots);
	Interpolant *Ldot_interp = new Interpolant(ys, es, Ldots);

	interps->Edot = Edot_interp;
	interps->Ldot = Ldot_interp;

}


// Class to carry gsl interpolants for the inspiral data
// also executes inspiral calculations
SchwarzEccFlux::SchwarzEccFlux(std::string few_dir)
{
    interps = new interp_params;

    // prepare the data
    // python will download the data if
    // the user does not have it in the correct place
    load_and_interpolate_flux_data(interps, few_dir);
	//load_and_interpolate_amp_vec_norm_data(&amp_vec_norm_interp, few_dir);
}

#define SchwarzEccFlux_num_add_args 0
#define SchwarzEccFlux_spinless
#define SchwarzEccFlux_equatorial
#define SchwarzEccFlux_file1 FluxNewMinusPNScaled_fixed_y_order.dat
__deriv__
void SchwarzEccFlux::deriv_func(double* pdot, double* edot, double* xdot,
                  double* Omega_phi, double* Omega_theta, double* Omega_r,
                  double epsilon, double a, double p, double e, double x, double* additional_args)
{
    if ((6.0 + 2. * e) > p)
    {
        *pdot = 0.0;
        *edot = 0.0;
        *xdot = 0.0;
        return;
    }

    SchwarzschildGeoCoordinateFrequencies(Omega_phi, Omega_r, p, e);
    *Omega_theta = *Omega_phi;

    double y1 = log((p -2.*e - 2.1));

    // evaluate ODEs, starting with PN contribution, then interpolating over remaining flux contribution

	double yPN = pow((*Omega_phi),2./3.);

	double EdotPN = (96 + 292*Power(e,2) + 37*Power(e,4))/(15.*Power(1 - Power(e,2),3.5)) * pow(yPN, 5);
	double LdotPN = (4*(8 + 7*Power(e,2)))/(5.*Power(-1 + Power(e,2),2)) * pow(yPN, 7./2.);

	double Edot = -epsilon*(interps->Edot->eval(y1, e)*pow(yPN,6.) + EdotPN);

	double Ldot = -epsilon*(interps->Ldot->eval(y1, e)*pow(yPN,9./2.) + LdotPN);

	*pdot = (-2*(Edot*Sqrt((4*Power(e,2) - Power(-2 + p,2))/(3 + Power(e,2) - p))*(3 + Power(e,2) - p)*Power(p,1.5) + Ldot*Power(-4 + p,2)*Sqrt(-3 - Power(e,2) + p)))/(4*Power(e,2) - Power(-6 + p,2));

    // handle e = 0.0
	if (e > 0.)
    {
        *edot = -((Edot*Sqrt((4*Power(e,2) - Power(-2 + p,2))/(3 + Power(e,2) - p))*Power(p,1.5)*
            	  (18 + 2*Power(e,4) - 3*Power(e,2)*(-4 + p) - 9*p + Power(p,2)) +
            	 (-1 + Power(e,2))*Ldot*Sqrt(-3 - Power(e,2) + p)*(12 + 4*Power(e,2) - 8*p + Power(p,2)))/
            	(e*(4*Power(e,2) - Power(-6 + p,2))*p));
    }
    else
    {
        *edot = 0.0;
    }

    *xdot = 0.0;
}


// destructor
SchwarzEccFlux::~SchwarzEccFlux()
{

    delete interps->Edot;
    delete interps->Ldot;
    delete interps;


}


// -------------------------------------------------------------------------------------------
// Initialize flux data for inspiral calculations
// void load_and_interpolate_flux_data_Kerr(struct interp_params *interps, const std::string& few_dir){

// 	// Load and interpolate the flux data
//     std::string fp = "/few/files/FluxesEdot_tensor_1PNnorm.dat";//"/few/files/FluxesEdotResc_u_a.dat";
//     fp = few_dir + fp;
// 	ifstream Flux_file(fp);
//     cout << few_dir << endl;

//     if (Flux_file.fail())
//     {
//         throw std::runtime_error("The file FluxesEdot_tensor_1PNnorm did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
//     }

// 	// Load the flux data into arrays
// 	string Flux_string;
// 	vector<double> ys, as, Edots;
// 	double y, a, Edot;
// 	while(getline(Flux_file, Flux_string)){

// 		stringstream Flux_ss(Flux_string);

// 		Flux_ss >> Edot >> y >> a;

// 		ys.push_back(y);
// 		as.push_back(a);
// 		Edots.push_back(Edot);

// 	}

// 	// Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)
// 	sort( ys.begin(), ys.end() );
// 	ys.erase( unique( ys.begin(), ys.end() ), ys.end() );

// 	sort( as.begin(), as.end() );
// 	as.erase( unique( as.begin(), as.end() ), as.end() );

//     // notice that if you resort ys and a you have to change also Edots

// 	Interpolant *Edot_interp = new Interpolant(ys, as, Edots);

// 	interps->Edot = Edot_interp;

// }

void load_and_interpolate_flux_data_Kerr(struct interp_params *interps, const std::string& few_dir){

	// Load and interpolate the flux data
    std::string fp = "/few/files/Kerr_flux_minus_PN1_l30_Has.dat";
    fp = few_dir + fp;
	ifstream Flux_file(fp);
    cout << few_dir << endl;

    if (Flux_file.fail())
    {
        throw std::runtime_error("The file Kerr_flux_minus_PN1_l30_Has did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
    }

	// Load the flux data into arrays
	string Flux_string;
	vector<double> ys, as, Edots_Kerr;
	double y, a, Edot_Kerr;
	while(getline(Flux_file, Flux_string)){

		stringstream Flux_ss(Flux_string);

		Flux_ss  >> y >> a >> Edot_Kerr;

		ys.push_back(y);
		as.push_back(a);
		Edots_Kerr.push_back(Edot_Kerr);

	}

	// Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)
	sort( ys.begin(), ys.end() );
	ys.erase( unique( ys.begin(), ys.end() ), ys.end() );

	sort( as.begin(), as.end() );
	as.erase( unique( as.begin(), as.end() ), as.end() );

    // notice that if you resort ys and a you have to change also Edots
        Interpolant *Edot_interp = new Interpolant(ys, as, Edots_Kerr);

        interps->Edot = Edot_interp;

}



// Class to carry gsl interpolants for the inspiral data
// also executes inspiral calculations
KerrCircFlux::KerrCircFlux(std::string few_dir)
{
    interps = new interp_params;

    // prepare the data
    // python will download the data if
    // the user does not have it in the correct place
    load_and_interpolate_flux_data_Kerr(interps, few_dir);
	//load_and_interpolate_amp_vec_norm_data(&amp_vec_norm_interp, few_dir);
}


double KerrCircFlux::EdotPN(double r, double a)
{
    double y = pow(1./(sqrt(r*r*r) + a), 2./3.) ;
    double res = 6.4*pow(y,5); //- 23.752380952380953*pow(y,6) + 1.6*(50.26548245743669 - 11.*a)*pow(y,6.5) + (-31.54215167548501 + 13.2*pow(a,2))*pow(y,7) + 0.009523809523809525*(-25732.785425553997 - 2646.*a - 504.*pow(a,3))*pow(y,7.5) + 
        //(-649.6614141423464 + 260.32427983539094*a + 163.36281798666926*pow(a,2) - 32.13333333333333*pow(a,3))*pow(y,8.5) + pow(y,8)*(740.6829867239124 - 217.8170906488923*a + 7.758730158730159*pow(a,2) - 52.17523809523809*log(y));//
        //pow(y,9)*(-748.828100625135 - 515.5802343491364*a + 69.31499118165785*pow(a,2) + 5.2*pow(a,4) + 3.2*sqrt(1. - 1.*pow(a,2)) + 41.6*pow(a,2)*sqrt(1. - 1.*pow(a,2)) + 19.2*pow(a,4)*sqrt(1. - 1.*pow(a,2)) + 12.8*(a + 3.*pow(a,3)) + 168.77786848072563*log(y));
    return res;
}

double EdotPN_global(double r, double a)
{
    double y = pow(1./(sqrt(r*r*r) + a), 2./3.) ;
    double res = 6.4*pow(y,5);
    return res;
}

#define KerrCircFlux_num_add_args 3
#define KerrCircFlux_equatorial
#define KerrCircFlux_circular
#define KerrCircFlux_file1 FluxesEdot_tensor_1PNnorm.dat
__deriv__
void KerrCircFlux::deriv_func(double* pdot, double* edot, double* xdot,
                      double* Omega_phi, double* Omega_theta, double* Omega_r,
                      double epsilon, double a, double p, double e, double x, double* additional_args)
{
    double p_sep = get_separatrix(a, e, x);
    if (p_sep > p)
    {
        *pdot = 0.0;
        *edot = 0.0;
        *xdot = 0.0;
        return;
    }

    double u = log(p - p_sep + 3.9);
    // cout << additional_args[0] << "\t" << additional_args[1] << "\t" << additional_args[2] <<endl;

    // Accretion effects
    double F = pow(1-Sqrt(p_sep/p),1./4.);
    // Accretion corrections 
    // Amplitude, n, m, totM
    //              Ldot_gw  * (M/1e5)**(6./5.)
    // logA
    // double F1 = 1-Sqrt(p_sep/p);
    // double  quota = (additional_args[2]*log(1 - Sqrt(p_sep)/Sqrt(10)))/4.;
    // double slope = (additional_args[1]/10. - (additional_args[2]*Sqrt(p_sep))/(8.*Sqrt(10)*(-10 + Sqrt(10)*Sqrt(p_sep))));
    // double approx = (p-10.0)*slope + quota
    // double true_eff = additional_args[1] * log(p/10) + additional_args[2] * log(F1) /4.0;

    // cout << (approx - true_eff)/true_eff <<endl;
    
    // double TotA = exp(-3.5 * log(p) + log(additional_args[0]) + additional_args[1] * log(p/10) + additional_args[2] * log(F1) /4.0);
    double LdotAcc = 32./5. * pow(p, -7./2.) * additional_args[0] * pow(p/10., additional_args[1]) * pow(F,additional_args[2]);

    // cout << p_sep << F << LdotAcc <<endl;
    // evaluate ODEs, starting with PN contribution, then interpolating over remaining flux contribution
	// double Edot = interps->Edot->eval(u, a);
    double y = pow(1./(sqrt(p*p*p) + a), 2./3.) ;
    double Edot = interps->Edot->eval(u, a)*pow(y, 6) + EdotPN(p, a); //

    // cout << 1./(sqrt(p*p*p) + a) - Omega_phi <<endl;
    // double omega = *Omega_phi;
    KerrGeoCoordinateFrequencies(Omega_phi, Omega_theta, Omega_r, a, p, e, x);
    // cout << "om1 = \t" << *Omega_phi << endl;
    // KerrGeoEquatorialCoordinateFrequencies(Omega_phi, Omega_theta, Omega_r, a, p, e, x);
    // cout << "om3 = \t" << *Omega_phi << endl;
    // double Omega_phi2 = 1.0/ (a + pow(p, 1.5) );
    // cout << "om2 = \t" << Omega_phi2 << endl;
    // cout << "rel diff = \t" << 1 - *Omega_phi / Omega_phi2 << endl;
    // *Omega_r = 0.0;
    // *Omega_theta = 0.0;
    double Ldot = Edot/ *Omega_phi + LdotAcc;
    double dL_dp = (-3*Power(a,3) + Power(a,2)*(8 - 3*p)*Sqrt(p) + (-6 + p)*Power(p,2.5) + 3*a*p*(-2 + 3*p))/(2.*Power(2*a + (-3 + p)*Sqrt(p),1.5)*Power(p,1.75));

	*pdot = -epsilon*Ldot/dL_dp; //LdotAcc / (32./5. * pow(p, -7./2.))

    // cout << "p_sep = " << p_sep <<endl;
    // cout << "F = " << F <<endl;
    // cout << "Ldot aa= "<<LdotAcc <<endl;
    // cout << "dL_dp = "<< dL_dp <<endl;

    *edot = 0.0;

    *xdot = 0.0;
}


// destructor
KerrCircFlux::~KerrCircFlux()
{

    delete interps->Edot;
    delete interps;

}