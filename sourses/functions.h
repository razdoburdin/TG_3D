#include "classes.h"
#include <stdio.h>
#include <math.h>

double coordinate(double j, double dz, double Lz);

double rho_bg(double z, double n);

double a2_bg(double z, double n);

double dp_0__dz__rho_0_bg(double z);

double drho_0__dz__rho_0_bg(double z,double n);

void chek_result_floder();

void timeformat(struct timeval time_start,struct timeval time_end,char **out);

void progress(int persent);

void get_bg_type(char** s);
