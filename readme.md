# hmhd2d_spc

This is a clone of original HMHD2D code (written by Yi-Min Huang), with some minor changes and addition of bootstrap-current effects. To compile, enter the build directory and run make, the executable will be in the same folder.

To run on SPC computers, load module "hdf5/1.10.4-intel19.1"


## Notices

If the code behaves strangely, lowering the optimization level may be one thing to try.

The user should provide a prob.f module, which defines the mesh, initial and boundary conditions. There are some existing modules for different physical systems in the package. One can copy the file (with filenames prof*.f) to prob.f

When use gfortran, the numbers in the input file should go like 1.D0, 1., 0.01, etc.  There must be a number following D. For example, 1.D  will confuse the code.

Parallel HDF5 can be slower than sequential one. Sequential output may 
  actually be prefered when memory is not a problem).  
  It is always a good practice to run a few small/short test runs for I/O 
  speed before throwing in big runs. 
  However, the problem of slow parallel output when binning data is largely
  solved in version 2.2

To enable hybrid MPI & OpenMP on NERSC, need to set 
  export MPICH_MAX_THREAD_SAFETY=multiple
  before executing the code  

Boundary condition calls of the same variable in different directions
  must be in the same OMP SECTION. Otherwise synchronization is NOT
  guaranteed.

To Do:
-- Add diagnostics for conservation laws checking

--------------------------------------------------------------------------
## Log

Version 3.21 - Fixed errors in bcz options 'ext4' and 'ext4_zero'.
               (Y.-M. Huang 11/19/2021)
             - add a boundary condition  option "ext2_sym'.
               (Y.-M. Huang 12/26/2021)
             - Slightly modified HISTCHEK routine.  
               (Y.-M. Huang 02/19/2022)
             - Fix a bug in "EXT3_ZERO" of BCZ 
               (Y.-M. Huang 02/28/2022)
             - Fix a typo nx->nz in step1d_z in stepping.
               (Y.-M. Huang 05/26/2022)

          

Version 3.20 - Make thermal conductivity depend on den, B, and T. 
               (Y.-M. Huang 11/20/2018)
             - Fixed an error in VDEPT_HDIF option in equation.f
               Options 2 & 3 were incorrectly placed.
               (Y.-M. Huang 07/11/2019)
               

Version 3.19 - Change the default boundary location from half-integer back 
               to full-integer grid point.
             - Clean up not not useful boundary conditon options.
               (Y.-M. Huang 08/09/2018)


Version 3.18 - Change the index order of FD coefficient arrays.
               (Y.-M. Huang 06/29/2018)
             - Add Ambipolar diffusion to the generalized Ohm's law
               (Y.-M. Huang 08/03/2018)
             - Add polynormial extrapolation options to BCs
             - Fix a typo in SMOM
               (Y.-M. Huang 08/07/2018)
             - Fix a typo in equation.f macro definition
               (Y.-M. Huang 08/08/2018)
             - Add more polynormial extrapolation options to BCs
               (Y.-M. Huang 08/09/2018)
               


Version 3.17 - Modify the friction part in SMOM so that the friction 
               profile is explicitly multiplied by the coefficient fric 
               there.
             - Allow problem specific parameters.
               (Y.-M. Huanng 02/22/2018) 
               
               (Y.-M. Huanng 02/22/2018) 
Version 3.16 - Revert back to the viscosity calculation of V.3.14
               (Y.-M. Huanng 01/16/2018
             - Move allocation of arrays to comm.f  
             - Add more measurements for HDF5 output time.  
               (Y.-M. Huanng 02/12/2018) 
               

Version 3.15 - Modify viscosity calculation.
               (Y.-M. Huang 06/11/2017)
             - Fixed 3-point 2nd order derivatives for uniform mesh
             - Add HDF WRITING TIME in output file
               (Y.-M. Huanng 01/04/2018)
             - Make the finite difference coefficients better sum
               to zero. 
               (Y.-M. Huanng 01/05/2018)               
             - Add an option to use allocatable arrays. This allows 
               using more memory in a given node, but may have positive 
               or negative impact on performace, depending on the machine.
               (Y.-M. Huanng 01/16/2018) 

Version 3.14 - Combine spux2, spuy2, spuz2 into smom2. 
             - Add viscous heating (Y.-M. Huang 06/10/2017)
             - Add vector version of hyperdiff and fix smom2 
             - Viscosity calculation is the same as HMHD3d V0.16
               (Y.-M. Huang 06/12/2017)

Version 3.13 - Do not use hdfwrite.f anymore. It is incompatible with
               hdf5-1.10.  (Y.-M. Huang 11/10/2016)
             - Add options to set minimum values of density and 
               temperature.
               (Y.-M. Huang 02/09/2017)

Version 3.12 - Add option to write Te and Ti in pdump. 
               (Yi-Min Huang 08/19/2016)


Version 3.11 - Clean up BCX and BCZ in SSUB.F.  Add 'asym1' option
               (Y.-M. Huang 04/25/2015)
             - Add auto-detection of compile machine and execute
               machines.
               (Y.-M. Huang 04/27/2015)
             - Fix a bug in calculating eta_perp  
               (Y.-M. Huang 04/28/2015)
             - Add anisotropic thermal conductivity to ion
               (Y.-M. Huang 04/29/2015)
             - Add temperature-dependent but isotropic resistivity  
               (Y.-M. Huang 04/30/2015)
             - Fix a bug that fin is not called after reading restart 
               (Y.-M. Huang 05/02/2015)
             - Fix a bug in bcx ext1 (Y.-M. Huang 05/03/2015)
             - Add 'asym2' to BC, and modify 'ext', and 'ext1'.
               (Y.-M. Huang 05/05/2015)
             - Add spatially dependent friction  
             - Merge full-integer boundary conditions into BCX and BCZ
               in 'symf', 'asymf', 'asymf2', 'extf1'
               (Y.-M. Huang 05/06/2015)
             - Move call restart to after setting up mesh and difference
               coefficients. 
               (Y.-M. Huang 05/12/2015)
             - Fix a bug that ihist was not initialized (hence was 0 but 
               should be 1) after setting initial conditions. This
               cause the first thist to be written at index-zero, which is 
               outside of the array and may pollute other variables if
               array bound checking is turned off.
               (Y.-M. Huang 06/19/2015)
             - Fix a typo in equation.f (Y.-M. Huang 06/25/2015)               

Version 3.10 - Add an option to time important kernels
             - Improve histchek & rsdump to calculate dump interval more 
               accurately. Make pdump sequence number continue 
               (instead of reset to 1) after restart. Change pdump 
               sequence number to 4 digits. Do not dump the 
               last frame unless pdump_last=.True.
               (Yi-Min Huang 09/05/2014)
             - Fix a bug in spre that Ohmic heating is included here even
               when pe is stepped separately.
               (Y.-M. Huang 09/06/2014)
             - No longer set idum=idum+iproc by default at bctparam.
               If the user wishes to use different random seeds 
               for different procs, that should be done when initialize
               in prob.f 
               (Y.-M. Huang 9/30/2014)
             - IMPORTANT !! Fix bugs in prob.f files that when 
               Alfven or magnetosonic wave speed are used for 
               hyper-diffusion coefficent, ghost cell of Bx and Bz 
               MUST be passed!
               (Y.-M. Huang 11/10/2014)
             - Change arguments of BCX, BCY,BCZ to strings for better
               readability. 
             - Modify 'ext1' in BCX and BCZ
               (Y.-M. Huang 04/24/2015)
               
               
               


Version 3.9  - Add an option to calculate Jy with (curl B)_y instead
               of -del2 psi. Calulation of Bx & Bz are moved to a new 
               subroutine CAL_B and the boundary condition
               for Bx and Bz are moved to a new subroutine FIN_B 
               (must be provided in PROB.F)
             - Add 3-pt stencil finite difference. Allow density diffusion
               to be calculated with 3-pt fD. Add option for calculating
               J  with 3-pt FD.
             - Add an option to implement resistivity as magnetic 
               diffusion rather than putting eta*J into the electric 
               field. 
             - Add an option to switch the calculation of magnetic force
               between using stress tensor or JxB. 
             - Add an option for Hall term calculation, either with
               JxB, or with the same way as the magnetic force.
             - Many improvements to hdftools.py.
               (Yi-Min Huang 02/26/2014)               
             - Fix a few typos in BC_Full/ssub.f
               (Yi-Min Huang 04/01/2014)  
             - Fix a typo in Equations.f
               (Yi-Min Huang 04/22/2014)
             - Write MPI configuration in output file.
               (Yi-Min Huang 09/05/2014)
               

Version 3.8  - Use OMP SECTIONS for halo exchange. Now each BC call is 
               handled by a single thread, but many BC calls can be 
               executed at the same time. To uniquely identify each
               MPI message, an extra "tag" is needed for each BC call,
               which needs to be an unique number within each fin_xxx
               function.
               (Yi-Min Huang 02/07/2014)
             - Change V-dependent coefficient abs(vm+vp) -> abs(vm)+abs(vp)
               This has effects only when vm and vp are of opposite sign.
               (Yi-Min Huang 02/18/2014)  


Version 3.7  - Add option to make dif4 spatially dependent. 
               (Yi-Min Huang 12/16/2013)
             - If RSDIR does not exist, use current directory instead.
               (Yi-Min Huang 12/17/2013)
             - Fix spatially-dependent dif4 in limited versions of 
               hyper-diffusion (Yi-Min Huang 01/21/2014)
             - Undo loop fusions in hyper-diffusion that will cause problem 
               in OpenMP (Yi-Min Huang 01/24/2014)
             - Some minor change in OpenMP implementation. Add OpenMP to 
               ssub.f of BC_Half (Yi-Min Huang 02/01/2014)  
             - Fix extrapolation BCs for BC_Half (Yi-Min Huang 02/03/2014)
Version 3.6  - Allow the simple limiter to switch on/off for each variable.
             - Add options to set pressure to zero when it becomes negative.
               (Yi-Min Huang 05/01/2012)  
             - Set a limit to the heat exchange rate, which is proportional
               to 1/di^2.  Excessive heat exchange rate casues numerical 
               instability but we only need something large enough to keep 
               Te~=Ti (Yi-Min Huang 9/13/2012)
             - Fix a typo in ssub.f (Yi-Min Huang 09/17/2012)

Version 3.5  - Add OpenMP to the code.  (Yi-Min Huang 08/31/2011)
             - Add gravity in Z (Yi-Min Huang 09/14/2011)
             - Modify io.f for some problems that seem to cause trouble in 
               Cray fortran compiler. [Each restart file needs a different id
               when reading and writing, also only proc0 write to the 
               output file.] (Yi-Min Huang 02/29/2012)
             - Add option to set RSDIR --- which is the directory the rstart 
               files reside.  This can be used to keep the home directory 
               clean. (Yi-Min Huang 03/05/2012)
             - Use MPI_Abort to ends the program when an error occurs.
               (Yi-Min Huang 03/14/2012) 
             - Do not write restart file should the code blow up.
               (Yi-Min Huang 03/15/2012)

Version 3.4  - Some tweak in performance.  Add option for Cray compiler, which
               gives significant performance gain on Hopper.
               (Yi-Min Huang 08/24/2011)

Version 3.3  - Remove DIAG.F90. The collection of diagnostics now go to 
               PROB.F (Yi-Min Huang 4/26/2011)
             - Add thermal exchange between ions and electrons. 
               Still primitive. The heat exchange rate is tied to resistivity.
               Since we have not distingished eta_perp and eta_para,
               there are some ambiguity here. (Yi-Min Huang 5/10/2011)
             - Add anisotropic resistivity and temperature-dependent 
               resistivity (Yi-Min Huang 5/23/2011)
             - change pt -> prei+pei, and add another variable ptb as 
               pt+byi^2/2 to avoid confusion (Yi-Min Huang 6/15/2011)
             - Add a simple limiter to ensure the hyperdiffusion is diffusive,
               or to ensure the hyperdiffusion won't enhance extrema.
               (Yi-Min Huang 6/20/2011)
             - Add hyperdiffusion to spuy2.  Previously this was missed.
               (Yi-Min Huang 7/28/2011)


Version 3.2  - Add Ohmic heating as an option. (Yi-Min Huang 04/12/2011)
             - Add an option to use coordinate transform instead of 
               interpolation for first order derivative 
               (Yi-Min Huang 04/12/2011)
             - Add INTEGRATE function in SSUB. (Yi-Min Huang 04/12/2011)             
             - Fix some bugs in velocity dependent hyper-diffusion 
               (mistakes in CPP preprocessing conditions)
             - Clean-up calculation of total pressure.
               (Yi-Min Huang 04/25/2011)
  
Version 3.1  - Move timestepping to a separate module.
               (Yi-Min Huang 03/09/2011)
             - Fix a bug in BCZ (typo nz->nx) (Yi-Min Huang 03/17/2011)

Version 3.0  - Move boundary to in between grid points. The canonical 
               boundary conditions are modified accordingly.
             - Add capability for a background potential psi.
             - Eliminate some MPI communicators that are only used by the 
               rotated Longcope problem.  Those should now go into the 
               prob.f file.  Two subroutine init_prob and finish_prob should 
               be defined in prob.f to initialize and finalize problem
               specific resources.  They can be left blank.
               (Yi-Min Huang 2/17/2011) 
             - Add another subroutine "restart" to prob.f, for problem 
               specific restart requirement.  (Yi-Min Huang 2/18/2011)
             - Change the default non-uniform grid mapping from cubic
               to fifth-order polynomial.  We enforce the second derivative
               of the mapping function to vanish at the edge of the simulation
               box.  This significantly reduces numerical noise when standard
               BCs (symmetric, antisymmetric, periodic) are employed.
               (Yi-Min Huang 3/01/2011)
             - Add output for div.v ('divv') and vorticity along y ('vorty')
               (Yi-Min Huang 3/09/2011)

Version 2.2  - Change integer kind to more portable way
               (Yi-Min Huang 11/24/2010)
             - Improve the efficiency of hdf5 binning output. The old
               problem of slow parallel output when binning data is
               largely solved.   
               (Yi-Min Huang 11/27/2010) 
             - Output the location and density value when divergence 
               happens. (Yi-Min Huang 11/29/2010) 
             - Add a script rstar to tar/gzip restart files
               (Yi-Min Huang 11/29/2010)          
             - Add a python script to read hdf5 data files. 
               (Yi-Min Huang 1/28/2010)
             - Fix a bug in z_exchange. (Yi-Min Huang 2/14/2011)
             - Change BYTE->IBYTE, WORD->IWORD, LWORD->ILWORD, otherwise
               gfortran will complain. (Yi-Min Huang 2/14/2011)


Version 2.1  - Move non-problem-specific code in prof.f to hmhd2d.f.
               Now prob.f only contains problem-specific code and therefore
               a lot cleaner. (Yi-Min Huang 11/20/2010)
             - Add another option of extrapolation based on Taylor expansion.
               (Yi-Min Huang 11/20/2010)
             - Add yet another option of extrapolation. 
               (Yi-Min Huang 11/22/2010)
             - Change integer kind to more portable way 
               (Yi-Min Huang 11/24/2010)

Version 2.0  - Cleanup version 1.x (Yi-Min Huang 10/30/2010)
             - Fix bugs in prob_coal_sym.f and prob_coal_quater.f
               (Yi-Min Huang 10/30/2010)
             - Minor bug fixing (Yi-Min Huang 11/20/2010)

Version 1.6  - Separate the physical problem to be solved (mesh, initial
               condition, boundary conditions) to a separate module PROB.
               (Yi-Min Huang 09/30/2010) 
             - Set intermediate times for multiple-step timestepping.
               This only has effects when the governing equation or 
               boundary conditions have time-dependence.  However, the 
               current implementation, although intuitive, may not be 
               optimal (see, Carpenter, Gottlieb, Abarbanel, and Don, 
               SIAM J. Scientific Compt. 1995, 1241-1252; Abarbanel,
               Gottlieb, and Carpenter, SIAM J. Sci. Comput. 1996, 777-782).  
               (Yi-Min Huang 10/04/2010).
             - Add 4-stage SSP RK3 as another timestepping method.
               (Yi-Min Huang 10/04/2010).                 
             - Get rid of the variables ifnonunifx & ifnonunifz.
               Now only use NONUNIFX and NONUNIFZ.  This to to 
               avoid ambiguity. (Yi-Min Huang 10/08/2010)
             - fix a bug in bczper (Yi-Min Huang 10/08/2010)
             - a dedicated variable rfamp for the random forcing amplitude.
               (Yi-Min Huang 10/08/2010)
             - set gamma=1 when ISOTHERMAL is defined.
               (Yi-Min Huang 10/12/2010)
             - Move the calculation of Velocity dependent hyperdiffusion
               coefficients fron FIN to a separate subroutine. Add the 
               option of adding wave speed to the flow speed.
               (Yi-Min Huang 10/13/2010)
             - Add a variable maxden, if density is above it, the code stops.
               (Yi-Min Huang 10/19/2010)
             - Fix a bug in bcx.  Add another option for extrapolation in
               bcx and bcz. (Yi-Min Huang 10/27/2010)
             - Add more options for velocity dependent hyperdiffusion.
               Now the wave speed can be choose from alfven speed, sound
               speed, and fast magnetosonic wave speed.
               (Yi-Min Huang 10/29/2010)

Version 1.5  - Add options for timestepping sheme. Trapezoidal leapfrog, 
               RK2, RK3.  We use SSP Runge-Kutta scheme (Gottlieb, Su, 
               and Tadmor SIAM Review 2001) (Yi-Min Huang 09/29/2010)

Version 1.4  - Get rid of a redundant assignment of byi -> temp.  
               (Yi-Min Huang 09/27/2010).

Version 1.3   - Add capbility of arranging procs into blocks.  This 
                potentially can improve data proximity.  For example,
                one may make the block size equal to the number of cores
                in one node, e.g. 4 cores -> block size = 2*2. 
                (Yi-Min Huang 06/24/2010) 
              - Change periodic BC back to the old convention. That is,
                grid point 3 is identical to grid point nxtot-2 instead 
                of nxtot-1.  This is necessary for potentials.  New periodic
                bc routines for potentials are added to ssub.f, 
                as bcxper_potential and bczper_potential.
                (Yi-Min Huang 07/01/2010)
              - Add a missing pei=prei when PE is not stepped.
                (Yi-Min Huang 07/11/2010)
              - Fix a bug in sequential output, caused by the newly 
                introduced block structure. (Yi-Min Huang 07/18/2010)   
              - Change extrapolation BC.  Add capability of having 
                different BCs on both sides. (Yi-Min Huang 08/27/2010).
              - Add hyper-resistivity (Yi-Min Huang 09/14/2010).
              - Add the option of communicate the interior only,
                leave the boundary for special treatment
                (Yi-Min Huang 09/24/2010).

Version 1.2   - When use sequential IO, allocate buffers in all PEs (not just
                PE0) even though those are not used.  This is to prevent
                problems on zaphod (Yi-Min Huang 04/23/2010)
              - Fix a bug in hyperdiffusion (Yi-Min Huang 05/19/2010)
              - Add iproc to idum, such that, each proc has a different 
                random seed. (Yi-Min Huang 05/28/2010)

Version 1.1   - Add di (ion skin depth), vex, vey, vez (electron velocity),
                pe (electron pressure), jxi, jzi
              - Add Hall effect and electron pressure gradient effects
              - Add my coalescence and random force (as an option)
              - Add options to dump v, ve, and j
                (Yi-Min Huang 4/04/2010)
              - Improve the way blow-up is checked (Yi-Min Huang 4/10/2010)    
              - Put data & coordinates/time into one file instead of two 
                separate files. (Yi-Min Huang 04/14/2010)
              - Add the capability of one-byte or two-byte binning. 
                (Yi-Min Huang 04/14/2010)

Version 1.0   - Improve HDF5 pdump.  This reduces some overhead, especially
                when using with many procs.
                (Yi-Min Huang 12/17/2009)
              - Comment out unnecessary communication for jyi. Also bxi and 
                bzi only when spuy is stepped.
                (Yi-Min Huang 12/18/2009) 
              - Fix a bug in coalbcs which is due to the modification of 
                periodic boundary conditions 
              - comment out vy=puy/den if BY is not stepped.
                (Yi-Min Huang 1/07/2010)
              - Change irsend --> isend.  This is for compatibility with 
                Bassi.  Also the performance seems slightly better in
                franklin and hopper. (Yi-Min Huang 1/09/2010)
              - Fix inc.h (Yi-Min Huang 1/13/2010)
              - Fix a bug in M4 macro for d/dz (Yi-Min Huang 1/13/2010)

Version 0.9   - Fix an inconsistency in periodic BC.  The 3rd grid point 
                is identical to nxtot-1 instead of nxtot-2 
              - Use mpi_irsend and mpi_irecv for preposting.  
                This allows better performance in Franklin.  
                Also now the number of blocks in each direction
                can be any number.  Previously it can only be either one or 
                an even numbers.   
                (Yi-Min Huang 12/16/2009)
Version 0.8   - Change the restart file name system to allow up to 99999 
                procs. The limit to 512 procs does not exist, at least in 
                franklin. (Yi-Min Huang 12/03/09)


Version 0.7   - Move CPP definitions to inc.h.  Now allows options to use
                simpler expression for derivatives when some directions are 
                known to be uniform. (Yi-Min Huang 10/15/2009) 
              - Slight change in the uniform grid finite difference.
                (Yi-Min Huang 11/11/09)
              - Slight change in hyper diffusion dx12 ->dx24, dz12 -> dz24
                vxm, fzm -> 2*vxm, 2*vzm and CPP
                (Yi-Min Huang 11/17/09)

Version 0.4   - Get rid of some redundancy in velocity dependent 
                hyperdiffusion.  vxmp,vxmm -> vxm, vzmp, vzmm -> vzm.
                This saves about 5% of the run time.
                 (Yi-Min Huang 9/25/09)
              - Further modify the hyperdif subroutine.  The 4th order 
                hyper diffusion coefficient now works differently from before.
                (Yi-Min Huang 9/28/09)
              - Now only proc0 check ifabort and ifend, then broadcast
                (Yi-Min Huang 10/03/09)

Version 0.3   - Fix a stupid typo in release 0.2  (Yi-Min Huang 9/25/09)


Version 0.2   - idum type -> INTEGER(4)
              - extrapolate jyi, bxi, bzi in fin.  Moves some terms 
                in sby2 and spuy2 to the "conservative terms"
                (Yi-Min Huang 09/05/2009)
              - fix a bug in bczasym (forget to set boundary values to zero)
                (Yi-Min Huang 09/11/2009)
              - add idum to nmlist (Yi-Min Huang 09/14/2009)
              - modify fin such that if velocity dependent hyperdif
                is not used, don't calculate vxmp, vxmm, vxmp, vzmm.
                (Yi-Min Huang 9/25/09)

Version 0.1   - First working MPI version. (Yi-Min Huang 09/01/2009)
              - Fix a bug of not broadcasting parameters. 
                (Yi-Min Huang 09/02/2009)
              - Add CPP preprocessor that make removing parts of the
                code easier (Yi-Min Huang 09/03/2009)  
              - This version should be identical to OMP version 0.4
