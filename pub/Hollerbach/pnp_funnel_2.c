/*  PNP2 library:

The first part of this library sets up PNP domains of two basic geometries:

(1) PNP is set up for a cylindrical domain with flat boundaries to the
bulks. The grid is uniform.

(2) PNP is set up for a domain comprising a cylindrical pore proper,
two atria of variable taper, and two hemispherical subdomains of
bulk. The grid is non-uniform, with a mesh width that increases in
proportion to domain area. Domain area is computed as the area of the
spherical surface that intersects the domain axis at the respective
grid node and is normal to the walls of the domain at its edge.
Atrial tapers are defined by client-provided functions that compute
wall angle (with domain axis) for a given array of axis grid nodes.

NOTE: geometry builders implicitly invoke build_PNP_ws and init_PNP_var.
      A builder itself may allocate and de-allocate some memory in addition to
      that needed for the PNP workspace.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pnp_funnel.h"

static void def_comp_TC(PNP_SMP *smp)
{
    int kion, k;

    for (kion = 0; kion < smp->Nions; ++kion) {
	for (k = 0; k < smp->Ngrid; ++k) {
	    smp->TCs[kion][k] = smp->Cs[kion][k];
	}
    }
}

/* Builder of plain cylinder geometry */

PNP_SMP *build_cyl_geom(
   int nions,			/* number of ions */
   int ngrid,			/* number of gridpoints */
   double poreradius,		/* Angstrom */
   double porelength)		/* Angstrom */
{
    PNP_SMP *smp;
    double temp, temp1;
    int k;

    smp = build_PNP_ws(nions, ngrid);
    temp = porelength/((double) (ngrid-1));
    temp1 = Pi*poreradius*poreradius;
    for (k = 0; k < smp->Ngrid; k++) {
	smp->X[k] = temp*k;
	smp->Area[k] = temp1;
	smp->Radius[k] = poreradius;
    }
    init_PNP_var(smp);
    smp->comp_TC = def_comp_TC;
    return(smp);
}

/* Builder of bulk/atrium/cylinder geometries */

/* builder library */

/* compute surface of diffusion front and atrial radius:
  input: xs, phi, n, na, poreradius
  output: X, R, A, n, na

On input, na and n specify the number of equidistant grid nodes in the atrial
and the total domains (pertaining to xs, phi); on output, na and n contain
the numbers of non-uniform grid nodes (pertaining to X, R, A).
*/

static void comp_surface(double *xs, double *phi, double *rs, double *rp,
			 double *as, double *hs, double *X, double *R,
			 double *A, int *n, int *na, double poreradius)
{
    double temp, temp1, *rsp, xkv, *phip, *rpp;
    double a0, twoPi, *hsp, *asp, xstep, Xmax, y_1, r_1, a_1, a, r, y;
    int k, kx, kp, kv, NA, gotcha, done;

    xstep = xs[1] - xs[0];
    rs[0] = poreradius;
    temp = tan(phi[0]);
    phip = phi + 1;
    rsp = rs + 1;
    for (k = 1; k < *na; ++k) {
	temp1 = tan(*phip++);
	*rsp = rsp[-1] + 0.5*(temp + temp1)*xstep;
	temp = temp1;
	++rsp;
    }
    a0 = Pi*poreradius*poreradius;
    twoPi = 2.0*Pi;
    rsp = rs;
    phip = phi;
    hsp = hs;
    asp = as;
    rpp = rp;
    for (k = 0; k < *na; k++) {
	*rpp = (*rsp++)/sin(*phip);
	*hsp = (*rpp)*(1.0 - cos(*phip++));
	*asp++ = (*rpp++)*(*hsp++)*twoPi;
    }

    /* center and wall points in atrium: radius is defined */
    Xmax = xs[*n-1];
    kp = kx = kv = 0;
    xkv = 0.0;

  atrium:
    while (kx < *n) {
	kx++;
	if (xs[kx] >= xkv) {
	    break;
	}
    }
    if (kx >= *na) {
	goto cone;
    }
    while (kp < *na) {
	if (xs[kp] + hs[kp] >= xkv) {
	    break;
	} else {
	    kp++;
	}
    }
    if (!kp) {
	y_1 = 0.0;
	a_1 = a0;
    } else {
	y_1 = xs[kp-1] + hs[kp-1] - xkv;
	a_1 = as[kp-1];
    }
    a = as[kp];
    y = xs[kp] + hs[kp] - xkv;
    X[kv] = xkv;
    A[kv] = (a_1*y - a*y_1) / (y - y_1);

    if (!kx) {
	y_1 = 0.0;
	r_1 = poreradius;
    } else {
	y_1 = xs[kx-1];
	r_1 = rs[kx-1];
    }
    r = rs[kx];
    y = xs[kx];
    R[kv] = (r_1*y - r*y_1) / (y - y_1);
    xkv += xstep*A[kv] / a0;
    kv++;
    goto atrium;

    /* wall point in atrium, center point in bulk: no radius defined */
  cone:
    NA = kv;

  cone_1:
    gotcha = 1;
    while (xs[kp] + hs[kp] < xkv) {
	if ((kp + 1) >= *na) {
	    gotcha = 0;
	    break;
	}
	kp++;
    }
    if (!gotcha) {
	goto bulk;
    }
    if (!kp) {
	y_1 = 0.0;
	a_1 = a0;
    } else {
	y_1 = xs[kp-1] + hs[kp-1] - xkv;
	a_1 = as[kp-1];
    }
    a = as[kp];
    y = xs[kp] + hs[kp] - xkv;
    X[kv] = xkv;
    A[kv] = (a_1*y - a*y_1)/(y - y_1);
    xkv += xstep*A[kv]/a0;
    kv++;
    goto cone_1;

    /* bulk */
  bulk:
    done = (xkv > Xmax);
    if (done) {
	xkv = Xmax;
    }
    X[kv] = xkv;
    temp = xkv - xs[kp];
    A[kv++] = a = temp*temp*twoPi;
    if (done) {
	*n = kv;
	*na = NA;
	return;
    }
    xkv += xstep*a/a0;
    goto bulk;
}

/* ------------------------------- front end --------------------------------
   lays out the geometry for a PNP computation whose domain comprises a pore
   (a central cylindrical section flanked on both sides by widening atria) and
   extends hemispherically into both bulk domains.
*/

PNP_SMP *build_cab_geom(
    int nions,			/* number of ions */
    int npore,			/* number of gridpoints in central section */
    double poreradius,		/* in Angstrom, of central section */
    double *lengths,		/* in Angstrom, of left bulk, left atrium,
				   center, right atrium, right bulk */
    void comp_L_angle(double *, int),	/* fn computing left atrial angles */
    void comp_R_angle(double *, int),	/* fn computing right atrial angles */
    int *indices,		/* indices of pore subdomains in profiles */
    int *dimensions)		/* dimensions of subdomains in profiles */
/*  The angle functions compute 'angle' for n equal intervals spanning
    the atrium, with index 0 at the central end (where the atrium meets
    the cylinder).
    NOTE that 'angle' values span the range from 0 to pi/4 EXCLUDING these
    exact boundary values (which lead to floating point overflow).
*/
{
    double xstep, *phi, *rs, *as, *rp, *hs, *LA, *RA, *LR, *RR, *LX, *RX;
    double *xs, *xsp, *Xp, *Ap, *Rp, *xp, *ap, *rpp, porearea, xoffs;
    char *mem, *mem_free;
    int nla, nra, nlb, nrb, nl, nr, nb, na, temp, n, k, ngrid;
    PNP_SMP *smp;

    xstep = lengths[2] / (npore - 1);
    nlb = lengths[0] / xstep + 1.0e-8;
    nla = lengths[1] / xstep + 1.0e-8;
    nra = lengths[3] / xstep + 1.0e-8;
    nrb = lengths[4] / xstep + 1.0e-8;
    na = (nla > nra) ? nla : nra;
    nl = nla + nlb;
    nr = nra + nrb;
    n = (nl > nr) ? nl : nr;
    nb = (7*na + 5*n)*sizeof(double);
    if ((mem_free = malloc(nb)) == NULL) {
	return(NULL);
    }
    mem = mem_free;

    temp = na*sizeof(double);
    phi = (double *) mem_free;	mem_free += temp;
    rs = (double *) mem_free;	mem_free += temp;
    rp = (double *) mem_free;	mem_free += temp;
    as = (double *) mem_free;	mem_free += temp;
    hs = (double *) mem_free;	mem_free += temp;
    LR = (double *) mem_free;	mem_free += temp;
    RR = (double *) mem_free;	mem_free += temp;

    temp = n*sizeof(double);
    xs = (double *) mem_free;	mem_free += temp;
    LX = (double *) mem_free;	mem_free += temp;
    RX = (double *) mem_free;	mem_free += temp;
    LA = (double *) mem_free;	mem_free += temp;
    RA = (double *) mem_free;	mem_free += temp;

    for (k = 0, xsp = xs; k < nl; k++) {
	*(xsp++) = xstep*(k + 1);
    }
    comp_L_angle(phi, nla);
    comp_surface(xs, phi, rs, rp, as, hs, LX, LR, LA, &nl, &nla, poreradius);
    for (k = 0, xsp = xs; k < nr; k++) {
	*(xsp++) = xstep*(k + 1);
    }
    comp_R_angle(phi, nra);
    comp_surface(xs, phi, rs, rp, as, hs, RX, RR, RA, &nr, &nra, poreradius);

    ngrid = npore - 2 + nl + nr;
    indices[0] = 0;
    indices[1] = nl - nla;
    indices[2] = nl - 1;
    indices[3] = nl - 1 + npore;
    indices[4] = indices[3] + nra - 1;
    dimensions[0] = nl - nla;
    dimensions[1] = nla - 1;
    dimensions[2] = npore;
    dimensions[3] = nra - 1;
    dimensions[4] = nr - nra;

#if 0
    printf("ngrid = %d\n",ngrid);
    for (k = 0; k < 5; k++) {
	printf("%d, %d\n", indices[k], dimensions[k]);
    }
#endif /* 0 */

    if ((smp = build_PNP_ws(nions, ngrid)) == NULL) {
	return(NULL);
    }
    Xp = smp->X;
    xp = LX + nl;
    xoffs = xp[-1];
    for (k = 0; k < nl; k++) {
	*Xp++ = xoffs - *(--xp);
    }
    for (k = 0; k < (npore-2); k++) {
	*Xp++ = xoffs + (k + 1)*xstep;
    }
    xoffs = Xp[-1] + xstep;
    xp = RX;
    for (k = 0; k < nr; k++) {
	*Xp++ = xoffs + *xp++;
    }
    Ap = smp->Area;
    ap = LA + nl;
    for (k = 0; k < nl; k++) {
	*Ap++ = *(--ap);
    }
    porearea = Pi * poreradius * poreradius;
    for (k = 0; k < (npore-2); k++) {
	*Ap++ = porearea;
    }
    ap = RA;
    for (k = 0; k < nr; k++) {
	*Ap++ = *ap++;
    }
    Rp = smp->Radius;
    for (k = nla; k < nl; k++) {
	*Rp++ = 0.0;
    }
    rpp = LR + nla;
    for (k = 0; k < nla; k++) {
	*Rp++ = *(--rpp);
    }
    for (k = 0; k < (npore-2); k++) {
	*Rp++ = poreradius;
    }
    rpp = RR;
    for (k = 0; k < nra; k++) {
	*Rp++ = *rpp++;
    }
    for (k = nra; k < nr; k++) {
	*Rp++ = 0.0;
    }
#if 0
    printf("%e, %e, %e\n", smp->X[0], smp->Radius[0], smp->Area[0]);
    printf("%e, %e, %e\n", smp->X[smp->Ngrid/2], smp->Radius[smp->Ngrid/2],
	smp->Area[smp->Ngrid/2]);
    printf("%e, %e, %e\n", smp->X[smp->Ngrid-1], smp->Radius[smp->Ngrid-1],
	smp->Area[smp->Ngrid-1]);
#endif /* 0 */

    init_PNP_var(smp);
    free(mem);
    smp->comp_TC = def_comp_TC;
    return(smp);
}
