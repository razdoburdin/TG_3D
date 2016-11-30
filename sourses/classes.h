//Copyright 2016 Dmitry N. Razdoburdin.
//
//This file is part of TG_3D. TG_3D is a program that calculats transient growth of linear perturbations in 3D shearing box by power iterations.
//
//TG_3D is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. TG_3D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with TG_3D. If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>

class parameters
{
private:
	
public:
	double n;
	double kx;
	double ky;
	double Topt;
	double dz;
	double q;
	double Lz;
	int Z;//Amount of points in z-direction
	double mu;
	double sigma;
	int nProcessors;
	double C;
	double cond1;//First condition of iterations ending
	double cond2;//Second condition of iterations ending
	int cond4;//Fourth condition of iterations ending

#if defined(LOGFILENAME) // Logfile
	FILE* LogFile;
#endif

	parameters (int, char**);
};

class background
{
public:
	std::vector <double> rho_0;
	std::vector <double> a2;
	std::vector <double> a2__rho_0__dz;
	std::vector <double> dp_0__dz__rho_0;

//with dz/2 ofset
	std::vector <double> rho_0_v;
	std::vector <double> a2_v;
	std::vector <double> drho_0__dz__rho_0_v;

	background (parameters); 

	void get_metric_type(char**);

	void get_boundary_type(char**);
};

class perturbation
{
private:
	double dt_calc (parameters, background,double, double);

	int forward (parameters, background, double, short);

	int backward (parameters, background, double, short);

	int boundary(parameters, background);

	double i_v (parameters, background, int);

	double i_v_x (parameters, background, int);

	double i_v_y (parameters, background, int);

	double i_v_z (parameters, background, int);

	double i_w (parameters, background, int);
	
public:
	std::vector <double> Im_v_x;
	std::vector <double> Im_v_y;
	std::vector <double> Re_v_z;
	std::vector <double> Re_w;

	double t;
	double t_kx;//tilde k_x
	double dt;

	perturbation (parameters, background, double);

	void initial_conditions (parameters);

	void write (parameters,background);

	int evolve (parameters, background, double, short);

	double norm (parameters, background);
		
	void norm_evolve (parameters, background, double, double);

	void gEz_evolve (parameters, background, double, double);

	void normalisation(parameters,background);

	void amplitude_unification(parameters,background);

	double kz_calc(parameters);

	void kz_max_calc(parameters, double*, double*);

	void norm_components_calc(parameters,  background, double*, double*, double*, double*, double*);

	void spectra_calc(parameters);

	void average_subtraction(parameters);

	double momentum_flux(parameters, background);

	double full_momentum_flux(parameters, background, double, double);
};

class optimal : public perturbation
{
private:
	double i_v (parameters, background, perturbation, int);

	double i_w (parameters, background, perturbation, int);

	double scalar_production (parameters, background, perturbation);

	int end_code;

	void save_vector(parameters, perturbation*);

	perturbation add(parameters, background, perturbation, double);

	void add_to_singular(parameters, background, perturbation, double);

	void singular_vectors_subtraction(parameters,background,optimal*,int);
	
public:

	double sigma;
		
	optimal(parameters, background, optimal*,int);

	void write_G(parameters, background, int);

	double Kz_0;

	double Kz_T;

	double kz_max_0;

	double F_kz_max_0;

	double kz_max_T;

	double F_kz_max_T;

//Weights of energy components
	double Ex_0;
	double Ex_T;
	double Ey_0;
	double Ey_T;
	double Ez_0;
	double Ez_T;
	double Ew_0;
	double Ew_T;

};

class G_T
{
private:
	std::vector <double> G;
	std::vector <double> kx;
	int max_i();

public:
	double Gmax;
	double kx_max;

	G_T(parameters, background);

	void write(parameters, background);
};
