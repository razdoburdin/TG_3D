#include <omp.h>
#include <sys/time.h>

#include "functions.h"

double optimal::i_v (parameters data, background bg, perturbation vec, int j)
{
	return bg.rho_0_v[j]*(Im_v_x[j]*vec.Im_v_x[j]+Im_v_y[j]*vec.Im_v_y[j]+Re_v_z[j]*vec.Re_v_z[j]);
}

double optimal::i_w (parameters data, background bg, perturbation vec, int j)
{
	double I=0;
	if ((j==0) || (j==data.Z-1))
	{
		if (data.n>1)
		{
			I=0;
		}
		if (data.n==1)
		{
			I=0.5;
		}
		if (data.n<1)
		{
			fprintf(stderr,"Problem at the boundary!");
			I=0;
		}
	}
	else
	{
		I=bg.rho_0[j]/bg.a2[j]*(Re_w[j]*vec.Re_w[j]);
	}
	return I;
}

double perturbation::i_v (parameters data, background bg, int j)
{
	return bg.rho_0_v[j]*(Im_v_x[j]*Im_v_x[j]+Im_v_y[j]*Im_v_y[j]+Re_v_z[j]*Re_v_z[j]);
}

double perturbation::i_v_x (parameters data, background bg, int j)
{
	return bg.rho_0_v[j]*(Im_v_x[j]*Im_v_x[j]);
}

double perturbation::i_v_y (parameters data, background bg, int j)
{
	return bg.rho_0_v[j]*(Im_v_y[j]*Im_v_y[j]);
}

double perturbation::i_v_z (parameters data, background bg, int j)
{
	return bg.rho_0_v[j]*(Re_v_z[j]*Re_v_z[j]);
}

double perturbation::i_w (parameters data, background bg, int j)
{
	double I=0;
	if ((j==0) || (j==data.Z-1))
	{
		if (data.n>1)
		{
			I=0;
		}
		if (data.n==1)
		{
			I=0.5;
		}
		if (data.n<1)
		{
			fprintf(stderr,"Энергия расходится на границе!");
			I=0;
		}
	}
	else
	{
		I=bg.rho_0[j]/bg.a2[j]*(Re_w[j]*Re_w[j]);
	}
	return I;
}

int perturbation::backward (parameters data,background bg,double t_end,short silence)
{
	struct timeval time_start,time_end;
	char* time;
	if (silence==0)
	{
		gettimeofday(&time_start,NULL);
		fprintf(stderr,"Backward: ");
	}
	
	std::vector <double> Im_dot_v_x(data.Z);
	std::vector <double> Im_dot_v_y(data.Z);
	std::vector <double> Re_dot_v_z(data.Z);
	std::vector <double> Re_dot_w(data.Z);

	int i,j,error,i_max;
	double dzpPowMinusOne=1.0/data.dz;

	int persent=0;
	double t_start=t;
	if (silence==0)
	{
		fprintf(stderr,"\t\t\t");
		progress(persent);
	}

	omp_set_num_threads(data.nProcessors);
	dt=dt_calc(data,bg,t_start,t_end);
	i_max=(int) ((t_start-t_end)/dt);

	for(i=0;i<i_max;i++)
	{
		t=t_start-i*dt;
		t_kx=data.kx+data.q*data.ky*t;
#pragma omp parallel private(j) shared(i)
		{
//Second half-layer
#pragma omp for SIMD schedule(auto)
			for(j=0;j<data.Z-1;j++)
			{
				Im_dot_v_x[j]=0;
				Im_dot_v_x[j]+=(2-data.q)*Im_v_y[j];
				Im_dot_v_x[j]-=t_kx*(Re_w[j+1]+Re_w[j])*0.5;

				Im_dot_v_y[j]=0;
				Im_dot_v_y[j]-=2*Im_v_x[j];
				Im_dot_v_y[j]-=data.ky*(Re_w[j+1]+Re_w[j])*0.5;

				Re_dot_v_z[j]=0;
				Re_dot_v_z[j]-=(Re_w[j+1]-Re_w[j])*dzpPowMinusOne;

				Im_v_x[j]-=dt*Im_dot_v_x[j];
				Im_v_y[j]-=dt*Im_dot_v_y[j];
				Re_v_z[j]-=dt*Re_dot_v_z[j];
			}

//First half-layer
#pragma omp for SIMD schedule(auto) WAIT
			for(j=1;j<data.Z-1;j++)
			{
				Re_dot_w[j]=0;
				Re_dot_w[j]+=bg.a2[j]*t_kx*(Im_v_x[j]+Im_v_x[j-1])*0.5;
				Re_dot_w[j]+=bg.a2[j]*data.ky*(Im_v_y[j]+Im_v_y[j-1])*0.5;
				Re_dot_w[j]-=bg.a2[j]*(Re_v_z[j]-Re_v_z[j-1])*dzpPowMinusOne;
				Re_dot_w[j]-=bg.dp_0__dz__rho_0[j]*(Re_v_z[j]+Re_v_z[j-1])*0.5;
				Re_w[j]-=dt*Re_dot_w[j];
			}
#pragma omp single nowait
			{
				error=boundary(data,bg);
			}
		}
		if (error>0)
		{
			return error;
		}
		
//Changing percent value
		if (silence==0)
		{
			if ((t_start-t)/(t_start-t_end)>(persent+1)*0.01)
			{
				while ((t_start-t)/(t_start-t_end)>(persent+1)*0.01)
				{
					persent++;
				}
				progress(persent);
			}
		}
	}

	if (silence==0)
	{
		progress(100);
		gettimeofday(&time_end,NULL);
		timeformat(time_start,time_end,&time);
		fprintf(stderr," [%d,%d] done (%s)\n",data.Z,i,time);
		delete[] time;
	}
	return error;
}

void background::get_metric_type(char** s)
{
	sprintf((*s),"3D");
}
