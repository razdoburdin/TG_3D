# Numerical code for calculations of transient growth of linear perturbations in 3D shearing box by power iterations.
See paper #PAPER for details

##Project structure:
+ configs/                                                                     *Configuration files:*
  +    compile_keys.conf                                             *Contain compilation parameters.*
  + link_keys.conf                                                         *Contain non-numerical parameters of calculation.*
  + params.conf                                                         *Contain numerical parameters of calculation.*

+ sourses/                                                                  *Sourses files:*
  + main.cpp                                                               *The main project file (uncoment one of five code blocks).*
  + classes.h                                                               *Header file with classes description.*
  + methods.cpp                                                        *Class methods which are common for all boundary conditions and metric type.*
  + methods_3D.cpp                                                 *Class methods for metric coincident with  acoustic energy of perturbations.*
  + methods_first.cpp                                               *Class methods for boundary conditions of first type.*
  + methods_second.cpp                                          *Class methods for boundary conditions of second type.*
  + methods_infinite.cpp                                             *Class methods for infinite flows.*
  + methods_periodic.cpp                                        *Class methods for periodic boundary conditions.*
  + functions.cpp                                                       *Functions which are common for all flow types.*
  + functions_homogeneous.cpp                            *Functions for homogeneous flow.*
  + functions_isothermal.cpp                                  *Functions for isothermal flow.*
  + functions_polytropic.cpp                                   *Functions for politropic flow.*
  + functions.h                                                           *Common header for all functions*.cpp files.*
  + procedures.cpp                                                    *Some procedures for main.cpp*
  + procedures.h                                                        *Header for procedures.cpp*

+ Makefile                                                                  *Makefile*
+ README                                                                  *This file*

##Compilation
To compile project run 'make' in terminal.

You need to have the following:
+ c++ compiler
+ make
+ gnu scientific library (gsl)
+ argtable

(Note that code was tested only for gcc)

##Runing the code
To start calculations run 'make task' in terminal.

##Configuration
###Non-numerical parameters (configs/link_keys.conf)
+ Choosing of the background flow:
BACKGROUND={homogeneous, polytropic, isothermal}
+ Choosing of the metric type (now only one type is avaliable)
METRIC=3D
+ Choosing of the boundary conditions:
BOUNDARY={first, second, periodic}
Boundary conditions of first type mean that W=0 at the boundaries, second --- dW/dz=0
+ Using vectorisation in parallel regions of the code (gcc >=4.9 is required)
VECTORIZE={yes, no}
+ Using additional conditions of iterations interruption (see paragraph "Conditions of iterations interruption" for details):
SIGNAL2={yes, no}
SIGNAL3={yes, no}
SIGNAL4={yes, no}
+ Check equality of (Aq, Aq)=(A^{\dag}Aq, q) after every iteration (may spend a lot of time)
TEST_OF_CONJUGATION={yes, no}
+ Output of short or full information at the file:
G_OUTPUT={short, full}
+ Write iterations log at the screen or save it into log file:
LOG_OUTPUT={stderr, TG_3D.log}

###Numerical parameters (configs/params.conf)
For setting numerical parametrs the follwing keys can be used:
+ Required parameters:
  + --ky, **double**, *Wave-number in y direction*
  + --dz, **double**, *Step of discretization*
  + --C, **double**, *Constant for CFL condition*
+ Not required parameters:
  + --n, **double**, *Polytropic index*, default: 1.5
  + --kx, **double**, *Wave-number in x direction*, default: 0
  + --Topt, **double**, *Optimisation time*, default: 1
  + --Lz, **double**, *Half-thickness of the isothermal flow*, , default: 1
  + --q, **double**, *Shear rate*, default: 1.5
  + --mu, **double**, *Position of initial condition*, default: 0.2
  + --sigma, **double**, *Size of initial condition*, default: 0.1
  + --cores, **int**, *Number of openmp threads (0 --- all avalible)*, default: 0
  + --cond1, **double**, *First conditions of iterations interruption*, default -5
  + --cond2, **double**, *Second conditions of iterations interruption*, default -6
  + --cond4, **int**, *Fourth conditions of iterations interruption*, default 500

##Output format
The transient amplification factor is recorded in file "G_$BACKGROUND_$METRIC_$BOUNDARY" in the following format (one lines for one calculation):
+ [1] Wave-number in x direction
+ [2] Wave-number in y direction
+ [3] Shear rate
+ [4] Polytropic index
+ [5] Step of discretization
+ [6] Constant for CFL condition
+ [7] Optimisation time
+ [8] Amplification factor
+ [9] First conditions of iterations interruption
+ [10] Condition that interrupted iterations during this calculations {1,2,3,4}.

By setting G_OUTPUT key to "full" additional information can be added to the output:
+ [11] F(kz=0, t=0)
+ [12] F(kz=0, t=Topt)
+ [13] kz_max(t=0)
+ [14] kz_max(t=Topt)
+ [15] F(kz_max(t=0), t=0)
+ [16] F(kz_max(t=Topt), t=Topt)
+ [17] Ex(t=0)
+ [18]    Ex(t=Topt)
+ [19] Ey(t=0)
+ [20] Ey(t=Topt)
+ [21] Ez(t=0)
+ [22] Ez(t=Topt)
+ [23] Ew(t=0)
+ [24] Ew(t=Topt)
+ [25] Number of curent singular value
Here F(kz, t) is a square of Fourier-amplitude of wave-number kz divided by sum of squares of all Fourier-amplitudes in the decomposition at time t,
kz_max(t) is wave-number corresponding to maximal F at time t,
E{x,y,z,w}(t) is energy component associated with {vx, vy, vz, w} in units of full energy of perturbation.

It is also possible to record perturbation profiles by using method "write" for classes "perturbation" and "optimal".
It records perturabtion in folder "result", file name format "q=%.3lf kx=%.2lf ky=%.2lf t=%.2lf".
Format of output:
+ [1] coordinate
+ [2] v_x
+ [3] v_y
+ [4] v_z
+ [5] w
+ [6] energy density associated with v_x, v_y and v_z
+ [7] energy density associated with w

##Conditions of iterations interruption:
To determine moment to interruption of iterations it is naturally to use determination of singular vector:
A^{\dag} A q = \sigma^2 q.
That leads to the first condition:
||A^{\dag} A q -\sigma^2 q||^2 / sigma^2 ** 10^{cond1}

Unfortunately sometimes first condition cannot be met with fixed dz and C due to numerical errors of integration.
Thats way it can be useful to use another criterions.
The second criterion interrupts iterations than changing of value X=||A^{\dag} A q -\sigma^2 q||^2 / sigma^2 become to slow:
(X_{i-1}-X{i})/X_{i-1}**10^{cond2}
The third one interrupts iterations when growth factor start decreasing.
And the fourth criterion interrupts iterations after cond4 iterations.
You can disable any of last three criterion in the link_keys.conf.
