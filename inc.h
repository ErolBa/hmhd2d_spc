/*Definitions for CPP preprocessor*/

/*Time a specific Kernel Only      */
/* 0: None                         */
/* 1: Timestepping                 */
/* 2: Restart output               */
/* 3: HDF output                   */
#define KERNEL_TIME 0


/*MPI-OpenMP Hybrid Mode*/
/* 0: MPI only          */
/* 1: MPI-OpenMP hybrid (default) */
#define HYBRID 1

/* Density Diffusion*/
/* 0: none             */
/* 1: Use 3-pt stencil */
/* 2: Use 5-pt stencil */
#define DIF_DEN  0


/*Hyperdiffusion for density*/
#define HDIF_DEN  1

/*Hyperdiffusion for Momentum*/
#define HDIF_MOM  1

/*Hyperdiffusion for Pressure*/
#define HDIF_PRE  1

/*Hyperdiffusion for Electron Pressure*/
#define HDIF_PE   1

/*Hyperdiffusion for B*/
#define HDIF_B    1

/*Velocity Dependent Hyperdiffusion                        */
/* 0: None                                                 */
/* 1: Use flow speed only                                  */
/* 2: Use flow speed + Sound Speed                         */
/* 3: Use flow speed + Alfven speed                        */
/* 4: Use flow speed + Fast magnetosonic wave speed        */
/*    this may work better for cases with strong shocks,   */
/*    such as Orszag-Tang                                  */ 

#define VDEPT_HDIF 1

/* 4-th order Hyperdiffusion                               */
/* 0: None                                                 */
/* 1: Uniform                                              */
/* 2: Allow Spatial variation according to pre-defined     */
/*    profile. This may help in dealing with tricky        */ 
/*    boundary conditions                                  */
#define HDIF4  0

/* Simple limiter                                            */
/* 0 : no limiter                                            */
/* 1 : If the hyperdiffusive flux amplifies local minimun    */
/*     or maximum, set to zero                               */
/* 2 : If the hyperdiffusive flux is not downgradient,       */
/*     set to 0 (stronger limiter)                           */
#define SIMPLE_LIMITER 0

/* Whether to use Limiter in each variable      */
/* 1: on, but only when SIMPLE_LIMITER>0        */ 
/* 0: off                                       */
/* It may be good to turn off limiter for B     */
/* to avoid excessive numerical reconnection    */
#define LIMITER_DEN 1
#define LIMITER_MOM 1
#define LIMITER_PRE 1
#define LIMITER_PE  1
#define LIMITER_B   0

/* Set Pressure to zero when negative */
#define POSITIVE_PRE 0
#define POSITIVE_PE  0

/* Set minimum density */
#define DEN_MIN 0

/* Set minimum temperature */
#define TEMP_MIN 0 



/*Gravity*/
#define GRAVITY 0

/* Fricition       */
/* 0: None                                                 */
/* 1: Uniform                                              */
/* 2: Allow Spatial variation according to pre-defined     */
/*    profile. This may help in dealing with tricky        */ 
/*    boundary conditions                                  */
#define FRICTION  0


/*Random Force*/
#define RANDOM_FORCE  0

/*Isothermal Equation of state    */
/*0: Step pressure                */
/*1: Isothermal, Assume Te=Ti     */
/*2: Isothermal, Te\=Ti           */ 
#define ISOTHERMAL  0


/*Jy Calculation  */
/*0: use -del2 psi   */
/*1: use d Bx/dz - d Bz/dx */
#define JY_CAL  0

/*Current Density Calculation      */
/* 0: Use 5-pt stencil             */
/* 1: Use 3-pt stencil             */
#define J_CAL 0


/*Magnetic Force Calculation                      */
/*0: Use div. stress tensor                       */
/*   (recommended, but may cause BC problems      */
/*1: Use JxB                                      */
#define JXB 0

/*Hall Term Calculation                           */
/*0: Use JxB (actually -ve x B)                   */
/*1: Use the same way as JXB calculation          */
/*   may be more stable when combining with       */
/*   calculating magnetic force with stree tensor */
#define HALL_CAL 0


/*Two and Half D, i.e. step By and Puy*/ 
#define TWO_AND_HALF  1

/*Step pe*/  
#define SPE  0

/*2-fluid effects. Two fluid implies Two and half D*/
#define TWOFLUID  0

/* Ambipolar Diffusion */
#define AMBIPOLAR 0

/*Include viscosity*/
#define VISCOSITY  1

/*Resistivity                           */
/*0: None                               */
/*1: Constant, Uniform                  */
/*2: Temperature-dependent, isotropic   */ 
/*3: Temperature-dependent, Anisotropic */ 
#define RESISTIVITY  1


/*Magnetic Diffusivity                              */
/*Implement resistivity in the form of eta*del2 B   */
/*Only valid for uniform resistivity                */
/*When turned on, the electric field in the code    */
/*does not include resistivity                      */ 
/*0: None                                           */
/*1: 5-pt stencil                                   */
/*2: 3-pt stencil                                   */
#define B_DIFF 0


/*Include hyper-resistivity*/
#define HYPER_RESISTIVITY  0

/*Include Viscous Heating*/
#define VISCOUS_HEATING 0


/*Include Ohmic heating*/
#define OHMIC_HEATING 0

/*Thermal exchange between ions and electrons*/
#define THERMAL_EXCHANGE 0

/*Thermal conductivity*/
#define THERMAL_CONDUCT 0

/*Add a potential background field*/
#define PSI0 0

/* Nonuniform in X*/
#define NONUNIFX  0

/* Nonuniform in Z */
#define NONUNIFZ  0

/* How first order derivative is treated.                                    */
/* 0 : 5 point interpolation                                                 */
/* 1 : Coordinate transform, d/dx=(dx/dx0)^(-1)*(d/dx0), where x0 is uniform.*/ 
#define FDTYPE 1

/* Time Stepping, 0: Trapezoidal leapfrog,                                   */ 
/* 1: SSP RK2, 2: SSP RK3 (Gottlieb and Shu 1998)                            */
/* 3: 4-stage SSP RK3 (Spiteri and Ruuth 2002)                               */ 
/* RK2 is generally unstable for waves, unless they are somehow damped.      */
/* Therefore RK2 is not very useful.                                         */
/* SSP RK3 in generally allows larger timesteps than trapezoidal leapfrog,   */
/* but is also more expansive per timestep (factor of 1.5).                  */
/* SSP RK3 is more efficient if the timestep is limited by dissipation,      */
/* as it allows a timestep about twice as large than trapezoidal leapfrog.   */
/* On the other hand, if timesteps are limited by waves, than trapzoidal     */
/* leapfrog has a slight edge, as the timestep in RK3 is only 1.2            */
/* times larger.                                                             */
/* 4-stage SSP RK3 is even more advantagous for dissipation limited timestep.*/ 
/* It is 33% more expensive per timestep compared to SSP RK3, but allows     */
/* a timestep about twice as large.                                          */
/* For wave limited problem, 4-stage SSP RK3 also allows slightly larger     */
/* timesteps, but probably not enough to offset the expense in each          */
/* timestep.                                                                 */
/* Overall RK3 should gives a better accuracy.                               */

#define TIMESTEPPING 0


/* If Two fluid, it has to be 2-1/2 D */
#if TWOFLUID==1
#define TWO_AND_HALF 1
#endif


/* 0 : use static arrays                                                 */
/* otherwise use ALLOCATABLE Arrays                                      */
/* Using Allocatable arrays allows for larger array sizes                */
/* Performance of allocatable arrays may be better or worse.             */
/* It is recommanded to do a few test runs to determine performance.     */

#define ALLOC 1

/*Problem Specific Parameters           */
/* 0: Only default qps                  */
/* Otherwise: Use probvars in prob.f    */

#define PROBVARS 0



