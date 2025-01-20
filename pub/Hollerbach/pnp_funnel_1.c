/* ============================ PNP1 library ================================ 

These functions let you set up a dynamic workspace adapted to your
current PNP problem (involving a number of ion species and a
geometrical grid) and solve the PNP equations in a given geometry and
using given boundary conditions and PNP parameters. The layout of a
pore geometry and the translation of far-bulk into domain-boundary
conditions are not topics of the PNP solver (they are facilitated by
other library modules).

PNP theory is solved here in one dimension, assuming homogeneity in
the others.  The metrics of the domain are expressed by an axial
coordinate and the area of the surface of homogeneity that belongs to
the axial coordinate.  Exact solutions are obtained for both
cylindrical and conical domains (with flat or spherically curved
surfaces, respectively). Approximate solutions are obtained for
domains that are piecewise cylindrical and conical (an example is a
domain that includes a cylindrical pore proper, tapered atria, and
hemispherical subdomains of the bulks). The PNP solver is given the
domain metrics as an array of axial grid coordinates and the
corresponding domain surface areas; it is otherwise unaware of the
shape of the domain, and thus capable of dealing with arbitrary
geometries. Furthermore, it works with the grid defined by the given
axial coordinates. This grid generally will be non-uniform (with
intervals made proportional to domain surface area).

The implemented theory makes provisions to account for two forms of
specific short-range interaction that ions may experience in a PNP
domain. (1) binding to saturable sites (which immobilizes the ion),
and (2) non-saturable partitioning (with potential effects on ion
mobility as expressed in the diffusion coefficient).  For instance,
form (1) may apply to protons that protonize structural carboxyl
groups, whereas form (2) may apply to alkali cations that interact
with an array of structural carbonyl groups. Form (2) is suitable also
for introducing Born energy variations associated with a non-uniform
dielectric permittivity (see below).

Several PNP parameters are allowed to vary along the axial coordinate
and need to be specified as axial profiles (on the grid of the axial
coordinate array):

  - area of domain surface (see above)
  - relative dielectric permittivity
  - diffusion coefficients of ions
  - 'partition energies' of ions (standard chemical potentials)
  - structural fixed-charge concentration

Saturable binding also is allowed to vary along the pore axis. This is
specified in a client-provided function that computes 'total' ion
concentrations along the pore from the given 'free' concentrations (a
simple equilibrium computation, since the 'free' concentrations in a
stationary flux are constant in time). A default function is provided
that merely sets the 'total' ion concentrations equal to the 'free'
ion concentrations, per species.

The PNP solver treats the bulk boundaries as providing
flux-independent ion concentrations and electric potentials. For PNP
domains that extend well into the bulks, the solver can use the
far-bulk values. Other boundary treatments (like Donnan or Guy-Chapman
theory) are possible; then the PNP solver is given the domain-side
boundary values produced by the respective boundary treatment.

The library comprises the front-end functions:

 build_PNP_ws:	dynamically allocate a workspace memory for the variables
		that transport information between the PNP solver and its
		clients. Several PNP workspaces can coexist and can be used
		randomly. The client is responsible for calling free_PNP_ws
		after use of the workspace memory is done.

 free_PNP_ws:	free a workspace memory after use

 init_PNP_var:	set up internal arrays for the PNP solver that need to be
		computed only once for many uses of the solver (but a change
		of geometry invalidates these internal arrays).

 solve_PNP:	solves the PNP equations for one parameter set

CAVEAT: This machinery defines its variables dynamically. You cannot expect
	your C compiler to protect you from indexing over bounds when you
	set up PNP input variables. Furthermore, the functions do not test
	your inputs: they execute them verbatim.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pnp_funnel.h"

/* Semaphore Constructor */

/*
   Unlike the previous version of this routine, the current version
   (a) uses malloc/free for dynamic memory management; and
   (b) does multiple allocations internally, to allocate memory of
   different types: this is so that any potential alignment restrictions
   are satisfied (malloc is guaranteed to satisfy all alignment restrictions)
*/

PNP_SMP *build_PNP_ws(int nions, int ngrid)
{
    int temp, kion, k;
    char *mem_free;
    PNP_SMP *smp;

    if ((smp = malloc(sizeof(PNP_SMP))) == NULL) {
	return(NULL);
    }
    smp->Nions = nions;
    smp->Ngrid = ngrid;

    smp->nbytes1 = 5*nions*sizeof(double *);
    if ((smp->mem1 = malloc(smp->nbytes1)) == NULL) {
	free(smp);
	return(NULL);
    }
    mem_free = smp->mem1;
    temp = nions*sizeof(double *);
    smp->Ds = (double **) mem_free;	mem_free += temp;
    smp->E0s = (double **) mem_free;	mem_free += temp;
    smp->Es = (double **) mem_free;	mem_free += temp;
    smp->Cs = (double **) mem_free;	mem_free += temp;
    smp->TCs = (double **) mem_free;	mem_free += temp;

    smp->nbytes2 = ((13 + 5*nions)*ngrid + 9*(ngrid - 2) + 4*nions)
		*sizeof(double);
    if ((smp->mem2 = malloc(smp->nbytes2)) == NULL) {
	free(smp->mem1);
	free(smp);
	return(NULL);
    }
    mem_free = smp->mem2;
    temp = ngrid*sizeof(double);
    smp->X = (double *) mem_free;	mem_free += temp;
    smp->Area = (double *) mem_free;	mem_free += temp;
    smp->Radius = (double *) mem_free;	mem_free += temp;
    smp->F = (double *) mem_free;	mem_free += temp;
    smp->EPS = (double *) mem_free;	mem_free += temp;
    smp->V = (double *) mem_free;	mem_free += temp;
    smp->WA = (double *) mem_free;	mem_free += temp;
    smp->WB = (double *) mem_free;	mem_free += temp;
    smp->WC = (double *) mem_free;	mem_free += temp;
    smp->WR = (double *) mem_free;	mem_free += temp;
    smp->WU = (double *) mem_free;	mem_free += temp;
    smp->WT = (double *) mem_free;	mem_free += temp;
    smp->H = (double *) mem_free;	mem_free += temp;
    for (k = 0; k < ngrid; ++k) {
	smp->X[k] = 0.0;
	smp->Area[k] = 0.0;
	smp->Radius[k] = 0.0;
	smp->F[k] = 0.0;
	smp->EPS[k] = 0.0;
	smp->V[k] = 0.0;
	smp->H[k] = 0.0;
    }
    for (kion = 0; kion < nions; ++kion) {
	smp->Ds[kion] = (double *) mem_free;	mem_free += temp;
	smp->E0s[kion] = (double *) mem_free;	mem_free += temp;
	smp->Es[kion] = (double *) mem_free;	mem_free += temp;
	smp->Cs[kion] = (double *) mem_free;	mem_free += temp;
	smp->TCs[kion] = (double *) mem_free;	mem_free += temp;
	for (k = 0; k < ngrid; ++k) {
	    smp->Ds[kion][k] = 0.0;
	    smp->E0s[kion][k] = 0.0;
	    smp->Es[kion][k] = 0.0;
	    smp->Cs[kion][k] = 0.0;
	    smp->TCs[kion][k] = 0.0;
	}
    }

    temp = nions*sizeof(double);
    smp->Z = (double *) mem_free;	mem_free += temp;
    smp->CL = (double *) mem_free;	mem_free += temp;
    smp->CR = (double *) mem_free;	mem_free += temp;
    smp->FLUX = (double *) mem_free;	mem_free += temp;
    for (kion = 0; kion < nions; ++kion) {
	smp->Z[kion] = 0.0;
	smp->CL[kion] = 0.0;
	smp->CR[kion] = 0.0;
	smp->FLUX[kion] = 0.0;
    }

    temp = (ngrid - 2)*sizeof(double);
    smp->al1 = (double *) mem_free;	mem_free += temp;
    smp->al2 = (double *) mem_free;	mem_free += temp;
    smp->al3 = (double *) mem_free;	mem_free += temp;
    smp->bt1 = (double *) mem_free;	mem_free += temp;
    smp->bt2 = (double *) mem_free;	mem_free += temp;
    smp->bt3 = (double *) mem_free;	mem_free += temp;
    smp->owa = (double *) mem_free;	mem_free += temp;
    smp->owb = (double *) mem_free;	mem_free += temp;
    smp->owc = (double *) mem_free;	mem_free += temp;

    return(smp);
}

/* Semaphore Destructor */

void free_PNP_var(PNP_SMP *smp)
{
    if (smp) {
	free(smp->mem2);
	free(smp->mem1);
	free(smp);
    }
}

/* ======================= Initializer of PNP variables ======================== 

   Initializes difference coefficients, voltage profile, and some PNP control
   variables. Use after building semaphore and before solving PNP.

   Note: this defaults the solver to the 'verbose' mode, in which it prints
   the maximal voltage change for each iteration as an index of convergence.
   The usefulness of this info hardly can be overstated.
*/

void init_PNP_var(PNP_SMP *smp)
{
    int k;
    double temp, *hp, *h1p, *al1p, *al2p, *al3p, *bt1p, *bt2p, *bt3p;

    for (k = 1; k < smp->Ngrid; ++k) {
	smp->H[k] = smp->X[k] - smp->X[k-1];
    }
    hp = smp->H + 1;
    h1p = hp + 1;
    al1p = smp->al1;
    al2p = smp->al2;
    al3p = smp->al3;
    bt1p = smp->bt1;
    bt2p = smp->bt2;
    bt3p = smp->bt3;
    for (k = 1; k < smp->Ngrid-1; ++k) {
	*al1p = 0.5/(*h1p);
	*al3p = -0.5/(*hp);
	*al2p++ = -(*al1p++ + *al3p++);
	temp = 2.0/(*h1p + *hp);
	*bt1p = temp/(*h1p++);
	*bt3p = temp/(*hp++);
	*bt2p++ = -(*bt1p++ + *bt3p++);
    }
    for (k = 0; k < smp->Ngrid; k++) {
	smp->V[k] = 0.0;
    }
    smp->MaxIter = 200;
    smp->Tolerance = 1e-8;
    smp->verbose = 1;
}

/* ========================== PNP solver ================================ */

/* ------------------------ solver library ------------------------------ */

/* find absolute maximum of an array */

static double maximum(double max, double *arr, int n)
{
    double temp;

    while (n-- > 0) {
	temp = *arr++;
	temp = (temp > 0.0) ? temp : -temp;
	max = (max > temp) ? max : temp;
    }
    return(max);
}

/*--- integrate array in place and return value of integral (results are *2) */

static double integrate(double *h, double *y, int n)
{
    double temp, sum;
    int k;

    sum = 0.0;
    temp = *y;
    *y++ = sum;
    h++;
    for (k = 1; k < n; k++) {
	sum += (temp + *y)* *h++;
	temp = *y;
	*y++ = sum;
    }
    return(sum);
}

/*--- compute energy and concentration profiles, and raw fluxes */

static void comp_EC_prof(PNP_SMP *smp)
{
    int kion, k;
    double *Ep, *Vp, *E0p, *WRp, *Dp, *Areap, *Cp, *cp, *wrp, *ep, z,
    temp, temp1, intL;

    for (kion = 0; kion < smp->Nions; kion++) {
	z = smp->Z[kion];
	Ep = smp->Es[kion];
	Vp = smp->V;
	E0p = smp->E0s[kion];
	WRp = smp->WR;
	Dp = smp->Ds[kion];
	Areap = smp->Area;
	temp = smp->kT_e;
	for (k = 0; k < smp->Ngrid; k++) {
	    *Ep = (*Vp++)*z + (*E0p++)/temp;
	    *WRp++ = exp(*Ep++)/((*Dp++)*(*Areap++));
	}
	intL = integrate(smp->H, smp->WR, smp->Ngrid);
	Cp = smp->Cs[kion];
	Ep = smp->Es[kion];
	temp = exp(Ep[0])*Cp[0];
	temp1 = exp(Ep[smp->Ngrid-1])*Cp[smp->Ngrid-1];
	cp = Cp + 1;
	wrp = smp->WR + 1;
	ep = Ep + 1;
	for (k = 1; k < smp->Ngrid-1; ++k) {
	    *cp++ = ((intL - *wrp)*temp + *wrp*temp1)/(exp(*ep++)*intL);
	    ++wrp;
	}
	smp->FLUX[kion] = (temp - temp1)*2.0/intL;
    }
    (*smp->comp_TC)(smp);
}

/*--- build discretized Poisson equations */

static void build_Poisson(PNP_SMP *smp)
{
    int kion, k;
    double *wap, *wbp, *wcp, *owap, *owbp, *owcp, *wrp, *Fp, temp, temp2,
    *Vp, z, z2, *vp, *tcp;

    Vp = smp->V;
    temp2 = 1.0/smp->Cscale;
    owap = smp->owa;
    owbp = smp->owb;
    owcp = smp->owc;
    wap = smp->WA + 1;
    wbp = smp->WB + 1;
    wcp = smp->WC + 1;
    wrp = smp->WR + 1;
    Fp = smp->F + 1;
    for (k = 1; k < smp->Ngrid - 1; k++) {
	*wap++ = *owap++;
	*wbp++ = *owbp++;
	*wcp++ = *owcp++;
	*wrp++ = -(*Fp++)*temp2;
    }
    for (kion = 0; kion < smp->Nions; kion++) {
	z = smp->Z[kion];
	z2 = z*z;
	wbp = smp->WB+1;
	wrp = smp->WR+1;
	vp = smp->V+1;
	tcp = smp->TCs[kion] + 1;
	for (k = 1; k < smp->Ngrid - 1; k++) {
	    temp = *tcp*z2;
	    *wbp++ -= temp;
	    *wrp++ -= temp*(*vp++) + z*(*tcp++);
	}
    }
    smp->WR[smp->Ngrid-2] -= smp->V[smp->Ngrid-1]*smp->WC[smp->Ngrid-2];
    smp->WR[1] -= smp->V[0]*smp->WA[1];
}

/*--- solve tridiagonal linear eqn system
   - wa holds the subdiagonal matrix elements (the first element is ignored)
   - wb holds the diagonal elements
   - wc holds the supradiagonal elements (the last element is ignored)
   - wr holds the right-hand sides
   - wu receives the solution
   - wg is an array for temporaries
   - returned int indicates whether system was solvable
*/

static int solvetridiag(PNP_SMP *smp)
{
    double *wap, *wbp, *wcp, *wrp, *wup, *wtp, bet;
    int j, n;

    wap = smp->WA + 1;
    wbp = smp->WB + 1;
    wcp = smp->WC + 1;
    wrp = smp->WR + 1;
    wup = smp->WU + 1;
    wtp = smp->WT + 1;
    n = smp->Ngrid - 2;

    if (wbp[0] == 0.0) {
	return(0);
    }
    bet = wbp[0];
    wup[0] = wrp[0]/bet;
    for (j = 1; j < n; j++) {
	wtp[1] = wcp[0]/bet;
	bet = wbp[1] - wap[1]*wtp[1];
	if (bet == 0.0) {
	    return(0);
	}
	wup[1] = (wrp[1] - wap[1]*wup[0])/bet;
	++wap;
	++wbp;
	++wcp;
	++wrp;
	++wup;
	++wtp;
    }
    j = n - 1;
    wup = smp->WU + j;
    wtp = smp->WT + n;
    while (j-- > 0) {
	wup[0] -= (*wtp--)*wup[1];
	--wup;
    }
    return(1);
}

/*--- compute error and new V */

static double new_V_and_err(PNP_SMP *smp)
{
    int k;
    double *wap, *wup, *vp;

    wap = smp->WA + 1;
    wup = smp->WU + 1;
    vp = smp->V + 1;
    for (k = 1; k < smp->Ngrid - 1; k++) {
	*wap++ = *wup - *vp;
	*vp++ = *wup++;
    }
    return(maximum(0.0, smp->WA + 1, smp->Ngrid - 2));
}


/* -------------------------- front end ------------------------------------

 - imports/exports data via semaphore structure (see PNP0.C)
 - the returned int signals convergence within the specified voltage
   tolerance
 - expects pore geometry to be defined in X, Area (you can use the geometry
   builder for setting these up)
 - expects boundary concentrations to be presented in CL, CR (it is your
   responsibility to compute these from bulk concentrations if, e.g.,
   surface charge effects are to be considered)
 - uses V profile from preceding instance as first approximation (adding
   a linear update for a change in VM); therefore, V should be nulled
   prior to the first instance (as is done automatically by the initializer
   function)
*/

int solve_PNP(PNP_SMP *smp)
{
    int k, kiter;
    double PoissonF, Fluxscale, temp, temp1;
    double *al1p, *al2p, *al3p, *bt1p, *bt2p, *bt3p;
    double *areap, *epsp, *owap, *owbp, *owcp, *Vp, *C;

/*-- design scaling factors */

    temp = maximum(0.0, smp->CL, smp->Nions);
    smp->Cscale = maximum(temp, smp->CR, smp->Nions);
    smp->kT = Boltzmann*(Celsius0 + smp->Celsius);
    smp->kT_e = smp->kT/e0;
    PoissonF = epsilon0*smp->kT_e*Liter /
	(e0*Avogadro*smp->Cscale*Angstrom*Angstrom);
    Fluxscale = cm2*Avogadro*smp->Cscale*Angstrom/Liter;

/*-- set up originals of left-hand arrays for Poisson */

    al1p = smp->al1;
    al2p = smp->al2;
    al3p = smp->al3;
    bt1p = smp->bt1;
    bt2p = smp->bt2;
    bt3p = smp->bt3;
    epsp = smp->EPS + 1;
    areap = smp->Area + 1;
    owap = smp->owa;
    owbp = smp->owb;
    owcp = smp->owc;
    for (k = 1; k < smp->Ngrid - 1; ++k) {
	temp = (epsp[0]/areap[0])
	    * ((*al1p)*areap[1] + (*al2p)*areap[0] + (*al3p)*areap[-1])
	    + ((*al1p)*epsp[1] + (*al2p)*epsp[0] + (*al3p)*epsp[-1]);
	*owcp++ = ((*epsp)*(*bt1p++) + temp*(*al1p++))*PoissonF;
	*owbp++ = ((*epsp)*(*bt2p++) + temp*(*al2p++))*PoissonF;
	*owap++ = ((*epsp)*(*bt3p++) + temp*(*al3p++))*PoissonF;
	++epsp;
	++areap;
    }

/*-- insert scaled boundary concentrations; update the V array of the previous
     instance for a change in VM (adding a linear correction to the profile)
*/

    for (k = 0; k < smp->Nions; k++) {
	C = smp->Cs[k];
	C[0] = smp->CL[k]/smp->Cscale;
	C[smp->Ngrid-1] = smp->CR[k]/smp->Cscale;
    }
    if ((temp = smp->VM/smp->kT_e - smp->V[0]) != 0.0) {
	Vp = smp->V;
	temp1 = smp->Ngrid - 1;
	for (k = 0; k < smp->Ngrid; k++) {
	    *Vp++ += temp*(temp1 - k)/temp1;
	}
    }

/*-- solve PNP by iteration */

    kiter = smp->MaxIter;
    while (kiter-- > 0) {
	comp_EC_prof(smp);
	build_Poisson(smp);
	if (!solvetridiag(smp)) {
	    return(0);
	}
	temp = new_V_and_err(smp);
	if (smp->verbose) {
	    fprintf(stderr, "dVmax = %.4e\n", temp);
	}
	if (temp < smp->Tolerance) {
	    break;
	}
    }

/* compute final concentration profiles, calibrated fluxes, and net current */

    comp_EC_prof(smp);
    smp->NetCurrent = 0.0;
    for (k = 0; k < smp->Nions; k++) {
	smp->FLUX[k] *= Fluxscale;
	smp->NetCurrent += smp->FLUX[k]*smp->Z[k]*e0/pA;
    }
    return(1);
}
