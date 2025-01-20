#ifndef PNP_FUNNEL_H
#define PNP_FUNNEL_H

/* ==================== PNP_in_C: header file ======================= */

/*------------------- PNP semaphore structure -------------------------

   - a simple pointer points to a profile array (dimension Ngrid)
     or an ion-related numerical array (dimension Nions)
   - a double pointer points to a pointer array (index: ion#)
     of ion-specific profiles
*/

struct pnp_smp_struct {

/* specification of memory used as PNP workspace */

    char *mem1;		/* pointer to base of chunk 1 of allocated */
    int nbytes1;	/* number of bytes in mem1 */
    char *mem2;		/* pointer to base of chunk 2 of allocated */
    int nbytes2;	/* number of bytes in mem2 */

/* descriptive arrays (not used in PNP computations) */

    double *Radius;	/* radius of pore               Angstrom */

/* basic dimensions */

    int Nions;		/* number of ions in system     -- */
    int Ngrid;		/* number of gridpoints         -- */

/* arrays submitted to PNP solver */

    double *X;		/* long-axis coordinate         Angstrom */
    double *Area;	/* area of diffusion front      Angstrom**2 */
    double *F;		/* fixed charge, signed         M */
    double *EPS;	/* rel dielectric permittivity  -- */
    double **Ds;	/* diffusion coefficients       cm**2/s */
    double **E0s;	/* basal free energies          eV */

    double *Z;		/* signed valencies             -- */
    double *CL;		/* left boundary conc's         M */
    double *CR;		/* right boundary conc's        M */

/* scalars submitted to PNP solver */

    double Celsius;	/* temperature                  Celsius */
    double VM;		/* transmembrane voltage        V */

    int MaxIter;	/* (Gummel) iteration limit     -- */
    double Tolerance;	/* tolerance                    kT/e */
    int verbose;	/* monitor iteration */

/* function for computing total concentrations */

    void (*comp_TC)(struct pnp_smp_struct *);

/* arrays computed by PNP solver */

    double **Cs;	/* free concentrations          Cscale */
    double **TCs;	/* total concentrations         Cscale */
    double *V;		/* voltage                      kT/e */
    double **Es;	/* el + basal free energies     kT */
    double *FLUX;	/* fluxes                       ions/s */

/* scalars computed by PNP solver */

    double Cscale;	/* concentration standard       M */
    double NetCurrent;	/* net pore current             pA */
    double kT;		/* kB*T                         AVs */
    double kT_e;	/* kB*T/e0                      V */

/* temporary arrays used by PNP solver */
    double *H;
    double *WA;
    double *WB;
    double *WC;
    double *WR;
    double *WU;
    double *WT;
    double *owa;
    double *owb;
    double *owc;
    double *al1;
    double *al2;
    double *al3;
    double *bt1;
    double *bt2;
    double *bt3;
};

typedef struct pnp_smp_struct PNP_SMP;

/* Natural Constants (from "A Physicist's Desk Reference", 2nd. ed.) */

#define Pi		3.141592653589793	/* -- */

#define epsilon0	8.85418781762e-12	/* As/(Vm) */
#define Boltzmann	1.380658e-23		/* AVs/Kelvin */
#define e0		1.60217733e-19		/* As */
#define Avogadro	6.0221367e23		/* 1/gmol */

#define Celsius0	273.15			/* Kelvin */

#define Angstrom	1.0e-10			/* m */
#define Liter		1.0E-3			/* m**3 */
#define cm2		1.0e-4			/* m**2 */
#define pA		1.0e-12			/* A */
#define mV		1.0e3			/* V */

/* prototypes of front end functions */

/* pnp_funnel_1.c */

PNP_SMP *build_PNP_ws(int nions, int ngrid);
void init_PNP_var(PNP_SMP *smp);
int solve_PNP(PNP_SMP *smp);

/* pnp_funnel_2.c */

PNP_SMP *build_cyl_geom(int nions, int ngrid, double poreradius, double porelength);

PNP_SMP *build_cab_geom(int nions, int npore, double poreradius, double *lengths, void comp_L_angle(double *, int), void comp_R_angle(double *, int), int *indices, int *dimensions);

#endif /* PNP_FUNNEL_H */
