#define STRINGIZE_NX(A) #A
#define STRINGIZE(A) STRINGIZE_NX(A)

#include <stdio.h>
#include <stdlib.h>
#include <argtable2.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <gsl/gsl_fft_real.h>

#include <sys/stat.h>
#include <dirent.h>
#include <sys/time.h>

#include "functions.h"

parameters::parameters(int N_data,char **char_data)
{
	struct arg_dbl *n_=arg_dbl0("n","n","<double>","Polytropic index. Required: no. Default: 1.5.");
	struct arg_dbl *kx_=arg_dbl0("k","kx","<double>","Wave-number in x direction. Required: no. Default: 0.");
	struct arg_dbl *ky_=arg_dbl1("k","ky","<double>","Wave-number in y direction. Required: yes.");
	struct arg_dbl *Topt_=arg_dbl0("T","Topt","<double>","Optimisation time. Required: no. Default: 1.");
	struct arg_dbl *dz_=arg_dbl1("z","dz","<double>","Step of discretization. Required: yes.");
	struct arg_dbl *C_=arg_dbl1("C","C","<double>","Constant for CFL condition. Required: yes.");
	struct arg_dbl *Lz_=arg_dbl0("L","Lz","<double>","Half-thickness of the isothermal flow. Required: no. Default: 1.0.");
	struct arg_dbl *q_=arg_dbl0("q","q","<double>","Shear rate. Required: no. Default: 1.5.");
	struct arg_dbl *mu_=arg_dbl0("q","mu","<double>","Position of initial condition. Required: no. Default: 0.2.");
	struct arg_dbl *sigma_=arg_dbl0("q","sigma","<double>","Size of initial condition. Required: no. Default: 0.1.");
	struct arg_int *cores_=arg_int0("c","cores","<int>","Number of openmp threads (0 --- all available). Required: no. Default: 0.");
	struct arg_dbl *cond1_=arg_dbl0("c","cond1","<double>","First conditions of iterations interruption. Required: no. Default: -5.");
	struct arg_dbl *cond2_=arg_dbl0("c","cond2","<double>","Second conditions of iterations interruption. Required: no. Default: -6.");
	struct arg_int *cond4_=arg_int0("c","cond4","<int>","Fourth conditions of iterations interruption. Required: no. Default: 500");
	struct arg_lit *help=arg_lit0("H","help","Help");
	struct arg_end *end=arg_end(16);

	n=1.5;
	kx=0;
	Topt=1;
	q=1.5;
	mu=0.2;
	sigma=0.1;
	Lz=1;
	cond1=-5;
	cond2=-6;
	cond4=500;
	nProcessors=omp_get_max_threads();

	void *argtable[]={n_,kx_,ky_,Topt_,dz_,C_,Lz_,q_,mu_, sigma_,cores_,cond1_,cond2_,cond4_,help,end};
	int nerrors=arg_parse(N_data,char_data,argtable);

	if ((*help).count>0)
	{
		arg_print_syntax (stdout, argtable, "\n\n");
		arg_print_glossary (stdout, argtable, "  %-10s %s\n");
		fprintf(stderr,"\t\t\t\t");
	}

	char* bg_type;
	bg_type=new char[15];
	get_bg_type(&bg_type);

	if (nerrors==0)
	{
		if ((*n_).count>0) n=(*n_).dval[0];
		if ((*q_).count>0) q=(*q_).dval[0];
		if ((*mu_).count>0) mu=(*mu_).dval[0];
		if ((*sigma_).count>0) sigma=(*sigma_).dval[0];
		if (((*Lz_).count>0) && (strcmp(bg_type,"isothermal")==0)) Lz=(*Lz_).dval[0];
		if (((*cores_).count>0) && ((*cores_).ival[0]>0)) nProcessors=(*cores_).ival[0];
		if ((*cond1_).count>0) cond1=(*cond1_).dval[0];
		if ((*cond2_).count>0) cond1=(*cond2_).dval[0];
		if ((*cond4_).count>0) cond4=(*cond4_).ival[0];
		kx=(*kx_).dval[0];
		ky=(*ky_).dval[0];
		Topt=(*Topt_).dval[0];
		dz=(*dz_).dval[0];
		C=(*C_).dval[0];
	}
	fprintf(stderr,"nerrors=%d\n",nerrors);

	delete [] bg_type;
	free(n_);
	free(kx_);
	free(ky_);
	free(Topt_);
	free(dz_);
	free(C_);
	free(mu_);
	free(sigma_);
	free(cores_);
	free(Lz_);
	free(q_);
	free(cond1_);
	free(cond2_);
	free(cond4_);
	free(help);
	free(end);

	Z=(int) (2*Lz/dz+1);
}

background::background (parameters data) : rho_0(data.Z), a2(data.Z), a2__rho_0__dz(data.Z), dp_0__dz__rho_0(data.Z), rho_0_v(data.Z), a2_v(data.Z), drho_0__dz__rho_0_v(data.Z)
{
	int j;
	double z;
	for(j=0;j<data.Z-1;j++)
	{
		z=coordinate(j,data.dz,data.Lz);
		rho_0[j]=rho_bg(z,data.n);
		a2[j]=a2_bg(z,data.n);
		a2__rho_0__dz[j]=a2_bg(z,data.n)/rho_bg(z,data.n)/data.dz;
		dp_0__dz__rho_0[j]=dp_0__dz__rho_0_bg(z);


		rho_0_v[j]=rho_bg(z+data.dz/2,data.n);
		a2_v[j]=a2_bg(z+data.dz/2,data.n);
		drho_0__dz__rho_0_v[j]=drho_0__dz__rho_0_bg(z+data.dz/2,data.n);
	}
	z=coordinate(data.Z-1,data.dz,data.Lz);
	rho_0[j]=rho_bg(z,data.n);
	a2[j]=a2_bg(z,data.n);
}

perturbation::perturbation (parameters data, background bg, double t_pert): Im_v_x(data.Z), Im_v_y(data.Z), Re_v_z(data.Z), Re_w(data.Z)
{
	t=t_pert;
	t_kx=data.kx+data.q*data.ky*t;
}

double perturbation::dt_calc (parameters data, background bg,double t_start, double t_end)
{
	double kz=2*M_PI/data.dz;
	double t_kx_start=data.kx+data.q*data.ky*t_start;
	double t_kx_end=data.kx+data.q*data.ky*t_end;

	double t_kx_max;
	if (t_kx_start*t_kx_start>t_kx_end*t_kx_end) t_kx_max=t_kx_start;
	else t_kx_max=t_kx_end;
	
	return data.C*2*M_PI*pow(bg.a2[(data.Z-1)/2]*(t_kx_max*t_kx_max+data.ky*data.ky+kz*kz),-0.5);
}

void perturbation::initial_conditions (parameters data)
{
	int j;
	double z;
	for(j=0;j<data.Z;j++)
	{
		z=coordinate(j,data.dz,data.Lz);
		Im_v_x[j]=0;
		Im_v_y[j]=0;
		Re_v_z[j]=0;
		Re_w[j]=exp(-pow(z-data.mu,2)*0.5*pow(data.sigma,-2));
	}

	t=0;
}

void perturbation::write (parameters data, background bg)
{
	chek_result_floder();
	FILE* output;
	char* name;
	name=new char[100];

	sprintf(name,"./result/q=%.3lf kx=%.2lf ky=%.2lf t=%.2lf",data.q,data.kx,data.ky,t);
	output=fopen(name,"w+t");

	int j;
	double z;
	for(j=0;j<data.Z-1;j++)
	{
		z=coordinate(j,data.dz,data.Lz);
		fprintf(output,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",z,Im_v_x[j],Im_v_y[j],Re_v_z[j],Re_w[j],i_v(data,bg,j)+i_w(data,bg,j));
	}
	z=coordinate(data.Z-1,data.dz,data.Lz);
	fprintf(output,"%lf\tgrid_ends\tgrid_ends\tgrid_ends\t%lf\n",z,Re_w[data.Z-1]);
	fclose(output);
	delete[] name;
}

int perturbation::forward (parameters data,background bg,double t_end,short silence)
{
	struct timeval time_start,time_end;
	char* time;
	if (silence==0)
	{
		gettimeofday(&time_start,NULL);
		fprintf(stderr,"Forward: ");
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
	i_max=(int) ((t_end-t_start)/dt);

	for (i=0;i<i_max;i++)
	{
		t=t_start+i*dt;
		t_kx=data.kx+data.q*data.ky*t;

#pragma omp parallel private(j)
		{
//First half-layer
#pragma omp for SIMD schedule(auto) WAIT
			for(j=1;j<data.Z-1;j++)
			{
				Re_dot_w[j]=0;
				Re_dot_w[j]+=bg.a2[j]*t_kx*(Im_v_x[j]+Im_v_x[j-1])*0.5;
				Re_dot_w[j]+=bg.a2[j]*data.ky*(Im_v_y[j]+Im_v_y[j-1])*0.5;
				Re_dot_w[j]-=bg.a2__rho_0__dz[j]*(bg.rho_0_v[j]*Re_v_z[j]-bg.rho_0_v[j-1]*Re_v_z[j-1]);
				Re_w[j]+=dt*Re_dot_w[j];
			}
#pragma omp single
			{
				error=boundary(data,bg);
			}

//Second half-layer
#pragma omp for SIMD schedule(auto) nowait
			for(j=0;j<data.Z-1;j++)
			{
				Im_dot_v_x[j]=0;
				Im_dot_v_x[j]+=2*Im_v_y[j];
				Im_dot_v_x[j]-=t_kx*(Re_w[j+1]+Re_w[j])*0.5;

				Im_dot_v_y[j]=0;
				Im_dot_v_y[j]-=(2-data.q)*Im_v_x[j];
				Im_dot_v_y[j]-=data.ky*(Re_w[j+1]+Re_w[j])*0.5;

				Re_dot_v_z[j]=0;
				Re_dot_v_z[j]-=(Re_w[j+1]-Re_w[j])*dzpPowMinusOne;

				Im_v_x[j]+=dt*Im_dot_v_x[j];
				Im_v_y[j]+=dt*Im_dot_v_y[j];
				Re_v_z[j]+=dt*Re_dot_v_z[j];
			}
		}
		if (error>0)
		{
			return error;
		}

//Changing percent value
		if (silence==0)
		{
			if ((t-t_start)/(t_end-t_start)>(persent+1)*0.01)
			{
				while ((t-t_start)/(t_end-t_start)>(persent+1)*0.01)
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
		fprintf(stderr," [%d,%d] done (%s)\n",data.Z,i_max,time);
		delete[] time;
	}

	return error;
}

int perturbation::evolve (parameters data, background bg, double t_end, short silence)
{
	int error=0;
	if (t<t_end) error=forward(data,bg,t_end,silence);
	else
	{
		if (t>t_end) error=backward(data,bg,t_end,silence);
	}
	return error;
}

double perturbation::norm (parameters data, background bg)
{
//Using of the trapezoidal rule to calculate integral. Because of the grid offset integral is divided into two parts.
	int j;
	double Iv,Iw;

	Iv=data.dz/2*(i_v(data,bg,0)+i_v(data,bg,data.Z-2));
	for(j=1;j<data.Z-2;j++)
	{
		Iv+=data.dz*i_v(data,bg,j);
	}

	Iw=data.dz/2*(i_w(data,bg,0)+i_w(data,bg,data.Z-1));
	for(j=1;j<data.Z-1;j++)
	{
		Iw+=data.dz*i_w(data,bg,j);
	}
	return Iv+Iw;
}

void perturbation::norm_evolve (parameters data, background bg,double t_end,double dT)
{
	struct timeval time_start,time_end;
	char* time;
	gettimeofday(&time_start,NULL);
	fprintf(stderr,"\nСчитаем g(t): ");
	
	chek_result_floder();
	FILE* output;
	char* name;
	name=new char[30];

	sprintf(name,"./result/g(t) t=(%.1lf,%.1lf)",t,t_end);
	output=fopen(name,"w+t");

	int persent=0;
	double t_start=t;
	fprintf(stderr,"\t\t\t");
	progress(persent);

	while (t<t_end)
	{
		fprintf(output,"%lf\t%lf\n",t,norm(data,bg));
		evolve(data,bg,t+dT,1);

//Changing percent value
		if ((t-t_start)/(t_end-t_start)>(persent+1)*0.01)
		{
			while ((t-t_start)/(t_end-t_start)>(persent+1)*0.01)
			{
				persent++;
			}
			progress(persent);
		}
	}
	fprintf(output,"%lf\t%lf\n",t,norm(data,bg));

	fclose(output);
	delete[] name;

	gettimeofday(&time_end,NULL);
	timeformat(time_start,time_end,&time);
	fprintf(stderr,"\tуспешно (%s)\n",time);
	delete[] time;
}

void perturbation::normalisation(parameters data,background bg)
{
	double I=sqrt(norm(data,bg));
	int j;
	for(j=0;j<data.Z;j++)
	{
		Im_v_x[j]=Im_v_x[j]/I;
		Im_v_y[j]=Im_v_y[j]/I;
		Re_v_z[j]=Re_v_z[j]/I;
		Re_w[j]=Re_w[j]/I;
	}
}

void perturbation::amplitude_unification(parameters data,background bg)
{
	double A=0;
	int j;
	for (j=0;j<data.Z;j++)
	{
		if (A<sqrt(Re_w[j]*Re_w[j])) A=sqrt(Re_w[j]*Re_w[j]);
	}
	for (j=0;j<data.Z;j++)
	{
		Im_v_x[j]=Im_v_x[j]/A;
		Im_v_y[j]=Im_v_y[j]/A;
		Re_v_z[j]=Re_v_z[j]/A;
		Re_w[j]=Re_w[j]/A;
	}
}

double perturbation::kz_calc(parameters data)
{
	double *Fv_x;
	Fv_x=(double*) malloc ((data.Z-1)*sizeof(double));
	double *Fv_y;
	Fv_y=(double*) malloc ((data.Z-1)*sizeof(double));
	double *Fv_z;
	Fv_z=(double*) malloc ((data.Z-1)*sizeof(double));
	double *Fw;
	Fw=(double*) malloc ((data.Z)*sizeof(double));

	int j;
	Fw[0]=Re_w[0];
	for (j=0;j<data.Z-1;j++)
	{
		Fv_x[j]=Im_v_x[j];
		Fv_y[j]=Im_v_y[j];
		Fv_z[j]=Re_v_z[j];
		Fw[j+1]=Re_w[j+1];
	}
	gsl_fft_real_wavetable * real;
	gsl_fft_real_workspace * work;
	work=gsl_fft_real_workspace_alloc(data.Z-1);
	real=gsl_fft_real_wavetable_alloc(data.Z-1);

	gsl_fft_real_transform(Fv_x,1,data.Z-1,real,work);
	gsl_fft_real_transform(Fv_y,1,data.Z-1,real,work);
	gsl_fft_real_transform(Fv_z,1,data.Z-1,real,work);
	gsl_fft_real_wavetable_free(real);
	gsl_fft_real_workspace_free (work);

	work=gsl_fft_real_workspace_alloc(data.Z);
	real=gsl_fft_real_wavetable_alloc(data.Z);
	gsl_fft_real_transform(Fw,1,data.Z,real,work);
	gsl_fft_real_wavetable_free(real);
	gsl_fft_real_workspace_free (work);

	double Fv_y_sum=0;
	for (j=0;j<data.Z-1;j++)
	{
		Fv_y_sum+=Fv_y[j]*Fv_y[j];
	}

	free(Fv_x);
	free(Fv_y);
	free(Fv_z);
	free(Fw);
	return Fv_y[0]*Fv_y[0]/Fv_y_sum;
}

void perturbation::spectra_calc(parameters data)
{
	double *Fv_x;
	Fv_x=(double*) malloc ((data.Z-1)*sizeof(double));
	double *Fv_y;
	Fv_y=(double*) malloc ((data.Z-1)*sizeof(double));
	double *Fv_z;
	Fv_z=(double*) malloc ((data.Z-1)*sizeof(double));
	double *Fw;
	Fw=(double*) malloc ((data.Z)*sizeof(double));

	int j;
	Fw[0]=Re_w[0];
	for (j=0;j<data.Z-1;j++)
	{
		Fv_x[j]=Im_v_x[j];
		Fv_y[j]=Im_v_y[j];
		Fv_z[j]=Re_v_z[j];
		Fw[j+1]=Re_w[j+1];
	}
	gsl_fft_real_wavetable * real;
	gsl_fft_real_workspace * work;
	work=gsl_fft_real_workspace_alloc(data.Z-1);
	real=gsl_fft_real_wavetable_alloc(data.Z-1);

	gsl_fft_real_transform(Fv_x,1,data.Z-1,real,work);
	gsl_fft_real_transform(Fv_y,1,data.Z-1,real,work);
	gsl_fft_real_transform(Fv_z,1,data.Z-1,real,work);
	gsl_fft_real_wavetable_free(real);
	gsl_fft_real_workspace_free (work);

	work=gsl_fft_real_workspace_alloc(data.Z);
	real=gsl_fft_real_wavetable_alloc(data.Z);
	gsl_fft_real_transform(Fw,1,data.Z,real,work);
	gsl_fft_real_wavetable_free(real);
	gsl_fft_real_workspace_free (work);

	FILE* Spectra;
	Spectra=fopen("./result/spectr_real","w+t");

	for (j=0;j<data.Z-1;j++)
	{
		fprintf(Spectra,"%lf\t%lf\t%lf\t%lf\t%lf\n",M_PI*j/data.Lz,Fv_x[j],Fv_y[j],Fv_z[j],Fw[j]);
	}

	free(Fv_x);
	free(Fv_y);
	free(Fv_z);
	free(Fw);
}

void perturbation::kz_max_calc(parameters data, double* kz_max, double* F_kz_max)
{
	double *Fv_y;
	Fv_y=(double*) malloc ((data.Z-1)*sizeof(double));

	int j;
	for (j=0;j<data.Z-1;j++)
	{
		Fv_y[j]=Im_v_y[j];
	}
	gsl_fft_real_wavetable * real;
	gsl_fft_real_workspace * work;
	work=gsl_fft_real_workspace_alloc(data.Z-1);
	real=gsl_fft_real_wavetable_alloc(data.Z-1);

	gsl_fft_real_transform(Fv_y,1,data.Z-1,real,work);
	gsl_fft_real_wavetable_free(real);
	gsl_fft_real_workspace_free (work);

	double Fv_y_sum=0;
	(*kz_max)=0;
	(*F_kz_max)=0;
	for (j=0;j<data.Z-1;j++)
	{
		Fv_y_sum+=Fv_y[j]*Fv_y[j];
		if ((*F_kz_max)<Fv_y[j]*Fv_y[j])
		{
			(*kz_max)=M_PI*j/data.Lz;
			(*F_kz_max)=Fv_y[j]*Fv_y[j];
		}
	}
	(*F_kz_max)=(*F_kz_max)/Fv_y_sum;

	for(j=0;j<10;j++)
	{
		fprintf(LOG,"F(kz=%.1lf pi)=%.8lf\n",j/data.Lz,Fv_y[j]*Fv_y[j]/Fv_y_sum);
	}

	free(Fv_y);
}

void perturbation::norm_components_calc(parameters data,  background bg, double* Ex, double* Ey, double* Ez, double* Ew)
{
//Using of the trapezoidal rule to calculate integral. Because of the grid offset integral is divided into two parts.
	int j;
	double Iv_x,Iv_y,Iv_z,Iw;

	Iv_x=data.dz/2*(i_v_x(data,bg,0)+i_v_x(data,bg,data.Z-2));
	Iv_y=data.dz/2*(i_v_y(data,bg,0)+i_v_y(data,bg,data.Z-2));
	Iv_z=data.dz/2*(i_v_z(data,bg,0)+i_v_z(data,bg,data.Z-2));
	for(j=1;j<data.Z-2;j++)
	{
		Iv_x+=data.dz*i_v_x(data,bg,j);
		Iv_y+=data.dz*i_v_y(data,bg,j);
		Iv_z+=data.dz*i_v_z(data,bg,j);
	}

	Iw=data.dz/2*(i_w(data,bg,0)+i_w(data,bg,data.Z-1));
	for(j=1;j<data.Z-1;j++)
	{
		Iw+=data.dz*i_w(data,bg,j);
	}
	(*Ex)=Iv_x/(Iv_x+Iv_y+Iv_z+Iw);
	(*Ey)=Iv_y/(Iv_x+Iv_y+Iv_z+Iw);
	(*Ez)=Iv_z/(Iv_x+Iv_y+Iv_z+Iw);
	(*Ew)=Iw/(Iv_x+Iv_y+Iv_z+Iw);
}

void perturbation::average_subtraction(parameters data)
{
	double Vx_av=0;
	double Vy_av=0;
	double W_av=0;
	int j;
	for (j=0;j<data.Z-2;j++)
	{
		Vx_av+=Im_v_x[j];
		Vy_av+=Im_v_y[j];
		W_av+=Re_w[j];
	}
	W_av+=Re_w[data.Z-1];
	Vx_av=Vx_av/(data.Z-1);
	Vy_av=Vy_av/(data.Z-1);
	W_av=W_av/data.Z;

	for (j=0;j<data.Z-2;j++)
	{
		Im_v_x[j]-=Vx_av;
		Im_v_y[j]-=Vy_av;
		Re_w[j]-=W_av;
	}
	Re_w[data.Z-1]-=W_av;
}

optimal::optimal(parameters data, background bg, optimal* singular_vectors, int N) : perturbation (data,bg,0.0)
{
	initial_conditions(data);
	int i;
	perturbation vec(data,bg,0);//q
	perturbation diff(data,bg,0);//A^+Aq-sigma^2q

	double scalar1=0;//scalar1=(Aq_i,Aq_i)=(q(t),q(t));
#if defined(TEST_OF_CONJUGATION)
	double scalar2=0;//scalar2=(A^+Aq_i,q_i)=(q_{i+1}(0),q_i(0));
#endif
	double scalar3=0;//scalar3=||A^+Aq-sigma^2q||^2
	double scalar4=0;//scalar3/scalar1 of the previous step
	double scalar5=0;//(scalar4-scalar3)/scalar4
#if defined(SIGNAL3)
	double scalar6=0;//scalar1 of the previous step
#endif
	end_code=0;
	i=0;

	while (end_code==0)
//	while (i<1)
	{
#if defined(LOGFILENAME)
		LOG=fopen(STRINGIZE(LOGFILENAME),"a+t");
#endif

		fprintf(LOG,"\nIteration %d\n",i+1);
		normalisation(data,bg);
		singular_vectors_subtraction(data,bg,singular_vectors,N);
		normalisation(data,bg);

		kz_max_calc(data,&kz_max_0,&F_kz_max_0);
		fprintf(LOG,"kz_max=%lf\tF(kz_max)=%lf\n",kz_max_0,F_kz_max_0);
//Forward
		save_vector(data,&vec);
		end_code=evolve(data,bg,data.Topt,0);

#if defined(SIGNAL3)
		scalar6=scalar1;
#endif
		scalar1=norm(data,bg);
		fprintf(LOG,"||Aq||^2=%lf\n",scalar1);

//Backward
		if (end_code==0)
		{
			end_code=evolve(data,bg,0,0);
		}

#if defined(TEST_OF_CONJUGATION)
		fprintf(LOG,"||A^{+}Aq||^2=%lf\n",norm(data,bg));
		scalar2=scalar_production(data,bg,vec);
		fprintf(LOG,"(A^{+}Aq,q)=%lf\n",scalar2);
		fprintf(LOG,"Conjugation=%lf\n",sqrt(pow((scalar1-scalar2)/scalar2,2)));
#endif

		diff=add(data,bg,vec,-scalar1);
		scalar4=scalar3/scalar1;
		scalar3=diff.norm(data,bg);
		fprintf(LOG,"lg ||A^{+}Aq-sigma^2 q||^2/sigma^2=%lf\n",log(scalar3/scalar1)/log(10));
		scalar5=(scalar4-scalar3/scalar1)/scalar4;
		if (scalar5<0) scalar5=-scalar5;
		fprintf(LOG,"lg change=%lf\n",log(scalar5)/log(10));
		
		if (scalar3/scalar1<pow(10,data.cond1) && i>5 && end_code==0)
		{
			end_code=1;
			sigma=scalar1;
		}
#if defined(SIGNAL2)
		if (scalar5<pow(10,data.cond2) && i>5 && end_code==0)
		{
			end_code=2;
			sigma=scalar1;
		}
#endif
#if defined(SIGNAL3)
		if (scalar6>scalar1 && i>5 && end_code==0)
		{
			end_code=3;
			sigma=scalar6;
		}
#endif
#if defined(SIGNAL4)
		if (i>data.cond4-2 && end_code==0)
		{
			end_code=4;
			sigma=scalar1;
		}
#endif

//		average_subtraction(data);
#if defined(LOGFILENAME)
		fclose(LOG);
#endif
		i++;
	}
#if defined(LOGFILENAME)
	LOG=fopen(STRINGIZE(LOGFILENAME),"a+t");
#endif
	fprintf(LOG,"Spectrum of optimal perturbation at t=0\n");
	Kz_0=kz_calc(data);
	kz_max_calc(data,&kz_max_0,&F_kz_max_0);
	norm_components_calc(data,bg,&Ex_0,&Ey_0,&Ez_0,&Ew_0);
	normalisation(data,bg);
	singular_vectors_subtraction(data,bg,singular_vectors,N);
	normalisation(data,bg);
	i=evolve(data,bg,data.Topt,0);
	fprintf(LOG,"Spectrum of optimal perturbation at t=Topt\n");
	Kz_T=kz_calc(data);
	kz_max_calc(data,&kz_max_T,&F_kz_max_T);
	norm_components_calc(data,bg,&Ex_T,&Ey_T,&Ez_T,&Ew_T);
	i=evolve(data,bg,0,0);
	normalisation(data,bg);
	singular_vectors_subtraction(data,bg,singular_vectors,N);
	normalisation(data,bg);
	fprintf(LOG,"Convergence. End code: %d\n",end_code);
#if defined(LOGFILENAME)
	fclose(LOG);
#endif
}

void optimal::singular_vectors_subtraction(parameters data,background bg,optimal* singular_vectors, int N)
{
		double fj;
		fprintf(LOG,"N=%d\n",N);
		int j;

		for (j=0;j<N;j++)
		{
			fj=scalar_production(data,bg,singular_vectors[j]);
			fprintf(LOG,"f[%d]=%lf\n",j,fj);
			add_to_singular(data,bg,singular_vectors[j],-fj);
		}
}

void optimal::write_G(parameters data,background bg,int j)
{
	FILE* G;
	char* bg_type;
	char* metric_type;
	char* boundary_type;
	char* name;
	bg_type=new char[15];
	metric_type=new char[10];
	boundary_type=new char[15];
	name=new char[50];

	get_bg_type(&bg_type);
	bg.get_metric_type(&metric_type);
	bg.get_boundary_type(&boundary_type);
	sprintf(name,"G_%s_%s_%s",bg_type,metric_type,boundary_type);
	G=fopen(name,"a+t");
	
	fprintf(G,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d",data.kx,data.ky,data.q,data.n,data.dz,data.C,data.Topt,sigma,data.cond1,end_code);
#if defined(G_OUTPUT_FULL)
	fprintf(G,"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",Kz_0,Kz_T,kz_max_0,kz_max_T,F_kz_max_0,F_kz_max_T,Ex_0,Ex_T,Ey_0,Ey_T,Ez_0,Ez_T,Ew_0,Ew_T);
	fprintf(G,"\t%d",j);	
#endif
	fprintf(G,"\n");
	
	fclose(G);
	delete[] bg_type;
	delete[] metric_type;
	delete[] boundary_type;
	delete[] name;
}

perturbation optimal::add(parameters data, background bg, perturbation vec, double a)
//Adding +a*vec
{
	int j;
	perturbation vec2(data,bg,0);
	for (j=0;j<data.Z-1;j++)
	{
		vec2.Im_v_x[j]=Im_v_x[j]+a*vec.Im_v_x[j];
		vec2.Im_v_y[j]=Im_v_y[j]+a*vec.Im_v_y[j];
		vec2.Re_v_z[j]=Re_v_z[j]+a*vec.Re_v_z[j];
		vec2.Re_w[j]=Re_w[j]+a*vec.Re_w[j];
	}
	vec2.Re_w[data.Z-1]=Re_w[data.Z-1]+a*vec.Re_w[data.Z-1];
	return vec2;
}

void optimal::add_to_singular(parameters data, background bg, perturbation vec, double a)
//Adding +a*vec 
{
	int j;
	for (j=0;j<data.Z-1;j++)
	{
		Im_v_x[j]=Im_v_x[j]+a*vec.Im_v_x[j];
		Im_v_y[j]=Im_v_y[j]+a*vec.Im_v_y[j];
		Re_v_z[j]=Re_v_z[j]+a*vec.Re_v_z[j];
		Re_w[j]=Re_w[j]+a*vec.Re_w[j];
	}
	Re_w[data.Z-1]=Re_w[data.Z-1]+a*vec.Re_w[data.Z-1];
}

double optimal::scalar_production(parameters data, background bg, perturbation vec)
{
	int j;
	double Iv,Iw;
	omp_set_num_threads(data.nProcessors);
	Iv=data.dz/2*(i_v(data,bg,vec,0)+i_v(data,bg,vec,data.Z-2));
#pragma omp parallel  for schedule(auto) private(j) shared(Iv) 
	for(j=1;j<data.Z-2;j++)
	{
		Iv+=data.dz*i_v(data,bg,vec,j);
	}

	Iw=data.dz/2*(i_w(data,bg,vec,0)+i_w(data,bg,vec,data.Z-1));
#pragma omp parallel  for schedule(auto) private(j) shared(Iw) 
	for(j=1;j<data.Z-1;j++)
	{
		Iw+=data.dz*i_w(data,bg,vec,j);
	}
	return Iv+Iw;	
}

void optimal::save_vector(parameters data, perturbation* vec)
{
	(*vec).t=t;
	(*vec).dt=dt;
	(*vec).t_kx=t_kx;

	int j;
	for (j=0;j<data.Z-1;j++)
	{
		(*vec).Im_v_x[j]=Im_v_x[j];
		(*vec).Im_v_y[j]=Im_v_y[j];
		(*vec).Re_v_z[j]=Re_v_z[j];
		(*vec).Re_w[j]=Re_w[j];
	}
	(*vec).Re_w[data.Z-1]=Re_w[data.Z-1];
}

G_T::G_T (parameters data, background bg) : G(6), kx(6)
{
	optimal* singular_vectors;
	singular_vectors=NULL;
	double kx_swing=-data.q*data.ky*data.Topt;

//Calculating of G(kx) in I points at the interval 2*kx_swing -- -kx_swing
	int i;
	int I=6;
	double dkx=-3*kx_swing/I;
	for(i=0;i<I+1;i++)
	{
		data.kx=2.0*kx_swing+dkx*i;
		fprintf(stderr,"kx=%lf\n",data.kx);
		optimal opt(data,bg,singular_vectors,0);
		opt.write_G(data,bg,1);
		G[i]=opt.sigma;
		kx[i]=data.kx;
	}

	double G_accur=0.01;//In units of G_max
	double kx_accur=0.01;//In units of kx_max

	size_t i_max;
	int end=0;

	while (end==0)
	{
		i_max=max_i();
		if (G[G.size()-1]>G[G.size()-2])
		{
			data.kx=kx.back()+5*dkx;
			fprintf(stderr,"kx=%lf\n",data.kx);
			optimal opt(data,bg,singular_vectors,0);
			opt.write_G(data,bg,1);
			G.resize(G.size()+1,opt.sigma);
			kx.resize(kx.size()+1,data.kx);		
		}
		else if (G[0]>G[1])
		{
			data.kx=kx.front()-5*dkx;
			fprintf(stderr,"kx=%lf\n",data.kx);
			optimal opt(data,bg,singular_vectors,0);
			opt.write_G(data,bg,1);
			G.insert(G.begin(),opt.sigma);
			kx.insert(kx.begin(),data.kx);
		}
		else if ((G[i_max]-G[i_max-1]>=G_accur*G[i_max]) && (kx[i_max]-kx[i_max-1]>=kx_accur*kx[i_max]))
		{
			data.kx=(kx[i_max]+kx[i_max-1])/2.0;
			fprintf(stderr,"kx=%lf\n",data.kx);
			optimal opt(data,bg,singular_vectors,0);
			opt.write_G(data,bg,1);
			G.insert(G.begin()+i_max,opt.sigma);
			kx.insert(kx.begin()+i_max,data.kx);
		}
		else if ((G[i_max]-G[i_max+1]>=G_accur*G[i_max]) && (kx[i_max+1]-kx[i_max]>=kx_accur*kx[i_max]))
		{
			data.kx=(kx[i_max]+kx[i_max+1])/2.0;
			fprintf(stderr,"kx=%lf\n",data.kx);
			optimal opt(data,bg,singular_vectors,0);
			opt.write_G(data,bg,1);
			G.insert(G.begin()+i_max+1,opt.sigma);
			kx.insert(kx.begin()+i_max+1,data.kx);
		}
		else
		{
			end=1;
			Gmax=G[i_max];
			kx_max=kx[i_max];
		}
	}
}

int G_T::max_i()
{
	double G_max=0;
	int i_max=0;

	size_t i;
	for (i=0;i<G.size();i++)
	{
		if (G[i]>G_max)
		{
			G_max=G[i];
			i_max=i;
		}
	}
	return i_max;
}

void G_T::write(parameters data, background bg)
{
	FILE* G;
	char* bg_type;
	char* metric_type;
	char* boundary_type;
	char* name;
	bg_type=new char[15];
	metric_type=new char[10];
	boundary_type=new char[15];
	name=new char[50];

//	optimal opt(data,bg,&end);
	get_bg_type(&bg_type);
	bg.get_metric_type(&metric_type);
	bg.get_boundary_type(&boundary_type);	
	sprintf(name,"G(T)_%s_%s_%s",bg_type,metric_type,boundary_type);
	G=fopen(name,"a+t");
	fprintf(G,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",data.ky,data.q,data.n,data.dz,data.C,data.Topt,Gmax,kx_max,data.cond1);

	delete[] bg_type;
	delete[] metric_type;
	delete[] boundary_type;
	delete[] name;	
	fclose(G);
}
