//Copyright 2016 Dmitry N. Razdoburdin.
//
//This file is part of TG_3D. TG_3D is a program that calculats transient growth of linear perturbations in 3D shearing box by power iterations.
//
//TG_3D is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. TG_3D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with TG_3D. If not, see <http://www.gnu.org/licenses/>.

#include "classes.h"
#include "procedures.h"

int main(int N_data, char** char_data)
{
	parameters data(N_data,char_data);
	background bg(data);

//Chose one of code blocks:

//1. Evolution of initial condition forward in time and then backward in time with recording of perturbation profiles after each evolution.
//	int end_code;
//	perturbation vec(data,bg,0);
//	vec.initial_conditions(data);
//	end_code=vec.evolve(data,bg,data.Topt,0);
//	vec.write(data,bg);
//	end_code=vec.evolve(data,bg,0,0);
//	vec.write(data);

//2. Evolution of normalised initial condition forward in time with writing full norm and its vertical component to the "result" folder.
//	double Ex,Ey,Ez,Ew,Norm;
//	perturbation vec(data,bg,0);
//	vec.initial_conditions(data);
//	vec.normalisation(data,bg);
//	vec.gEz_evolve(data,bg,100,0.1);
//	vec.norm_components_calc(data,bg,&Ex,&Ey,&Ez,&Ew,&Norm);
//	fprintf(stderr,"Ex=%lf\tEy=%lf\tEz=%lf\tEw=%lf\tNorm=%lf\n",Ex,Ey,Ez,Ew,Norm);

//3. Calculation of integrated over time momentum flux through the box.
//	perturbation vec(data,bg,0);
//	vec.initial_conditions(data);
//	vec.normalisation(data,bg);
//	fprintf(stderr,"Flux=%lf\n",vec.full_momentum_flux(data,bg,200,0.1));
	
//4. Calculating of N singular values.
	optimal* singular_vectors;
	int j;
	int N=1;
	singular_vectors=(optimal*) malloc (N*sizeof(optimal));
	for(j=0;j<N;j++)
	{
		singular_vectors[j]=optimal(data,bg,singular_vectors,j);
		singular_vectors[j].write_G(data,bg,j+1);
	}

//5. Calculating of G as function of kx.
//	double kx_min=-25*data.ky;
//	double kx_max=10*data.ky;
//	G_kx(data,bg,kx_min,kx_max,0.5*data.ky);

//6. Calculating maximal G over all kx.
//	G_T G(data, bg);
//	G.write(data,bg);

	return 0;
}
