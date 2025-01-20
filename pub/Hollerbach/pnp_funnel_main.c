/* ======================== PNP tester =============================== */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pnp_funnel.h"

static void do_angle(double *phi, int n)
{
    int k;
    double pi_2 = Pi/2;

    for (k = 0; k < n; k++) {
	phi[k] = pi_2*((double) k + 1)/((double) n + 2);
    }
}

static void constant_angle(double *phi, int n)
{
    int k;
    double pi_4 = Pi/4.0;

    for (k = 0; k < n; k++) {
	phi[k] = pi_4;
    }
    phi[0] = 0.5*pi_4;
    phi[n-1] = 1.5*pi_4;
}

int main(int argc, char **argv)
{
    PNP_SMP *smp;
    char line[100];
    int sel, k, kf, kl, honor_break;
    double donnan_l, donnan_r, temp, v_lo, v_hi, v_step;

    int nions, ngrid, indices[5], dimensions[5];
    double poreradius, porelength, lengths[5];

    fprintf(stderr, "Possibilities are:\n    0\tquit this program\n\n    1\tsee allocation of cylindrical geometry (prints some stats)\n    2\tsee setup of cylindrical geometry (normally silent)\n    3\tsee running of cylindrical geometry (see Gummel iteration progress)\n    4\tsee normal order of cylindrical geometry: 1 then 2 then 3\n\n    11\tsee allocation of angled geometry (prints some stats)\n    12\tsee setup of angled geometry (normally silent)\n    13\tsee running of angled geometry (see Gummel iteration progress)\n    14\tsee normal order of angled geometry: 11 then 12 then 13\n\n    21\ta full IV curve for a cylindrical pore with donnan potentials,\n\twith user-specified parameters\n\n    31\ta full IV curve for a cylindrical pore with 45-degree atria\n\tand bulk bath, with user-specified parameters\n\nNormal progression of calls is 1, 2, 3, or 4; or 11, 12, 13, or 14; or 21; or 31\n\n");
    while (1) {
	fprintf(stderr, "enter choice: ");
	fgets(line, 100, stdin);
	if (sscanf(line, "%d", &sel) == 1) {
/* warning warning warning!
   The code in this switch messes with the flow of control; fall-through,
   where it occurs, is intentional, and the order of the cases must not
   be changed, nor should additional cases be added between existing ones,
   except where it's explicitly noted as being safe to do so.
   warning warning warning! */
	    honor_break = 1;
	    switch(sel) {
/* it is safe to insert other cases here */
	      case 0:
		goto done_loop;
/* it is safe to insert other cases here */
	      case 4:
		honor_break = 0;
	      case 1:
		nions = 2;
		ngrid = 1000;
		poreradius = 3.0;
		porelength = 10.0;
		smp = build_cyl_geom(nions, ngrid, poreradius, porelength);
		fprintf(stderr, "mem = %lx, (%lx, %d), (%lx, %d)\n",
		    (unsigned long) smp,
		    (unsigned long) smp->mem1, smp->nbytes1,
		    (unsigned long) smp->mem2, smp->nbytes2);
		if (honor_break) {
		    break;
		}
	      case 2:
		smp->Z[0] = -1.0;
		smp->Z[1] = 1.0;
		smp->CL[0] = 0.3;
		smp->CL[1] = 0.3;
		smp->CR[0] = 0.3;
		smp->CR[1] = 0.3;
		smp->VM = -0.1;
		smp->Celsius = 10.0;
		for (k = 0; k < smp->Ngrid; ++k) {
		    smp->F[k] = -5.0;
		    smp->EPS[k] = 80.0;
		    smp->Ds[0][k] = 2.0e-8;
		    smp->Ds[1][k] = 1.0e-8;
		    smp->E0s[0][k] = 0.0;
		    smp->E0s[1][k] = 0.0;
		}
		if (honor_break) {
		    break;
		}
	      case 3:
		smp->MaxIter = 50;
		smp->Tolerance = 1.0e-10;
		if (solve_PNP(smp)) {
		    printf("NetCurrent = %+.9e\n", smp->NetCurrent);
		    printf("Fluxes = %+.9e, %+.9e\n",
			smp->FLUX[0], smp->FLUX[1]);
		} else {
		    printf("No convergence!\n");
		}
		break;
/* it is safe to insert other cases here */
	      case 14:
		honor_break = 0;
	      case 11:
		nions = 2;
		ngrid = 1000;
		poreradius = 3.0;
		lengths[0] = 100;
		lengths[1] = 20;
		lengths[2] = 10;
		lengths[3] = 20;
		lengths[4] = 100;
		smp = build_cab_geom(nions, ngrid, poreradius, lengths,
		    do_angle, do_angle, indices, dimensions);
		fprintf(stderr, "mem = %lx, (%lx, %d), (%lx, %d)\n",
		    (unsigned long) smp,
		    (unsigned long) smp->mem1, smp->nbytes1,
		    (unsigned long) smp->mem2, smp->nbytes2);
		if (honor_break) {
		    break;
		}
	      case 12:
		smp->Z[0] = -1.0;
		smp->Z[1] = 1.0;
		smp->CL[0] = 0.3;
		smp->CL[1] = 0.3;
		smp->CR[0] = 0.3;
		smp->CR[1] = 0.3;
		smp->VM = -0.1;
		smp->Celsius = 10.0;
		for (k = 0; k < smp->Ngrid; ++k) {
		    smp->F[k] = 0.0;
		    smp->EPS[k] = 80.0;
		    smp->Ds[0][k] = 2.0e-6;
		    smp->Ds[1][k] = 1.0e-6;
		    smp->E0s[0][k] = 0.0;
		    smp->E0s[1][k] = 0.0;
		}
		kf = indices[2];
		kl = kf + dimensions[2];
		for (k = kf; k < kl; ++k) {
		    smp->F[k] = -5.0;
		    smp->Ds[0][k] = 2.0e-8;
		    smp->Ds[1][k] = 1.0e-8;
		}
		if (honor_break) {
		    break;
		}
	      case 13:
		smp->MaxIter = 50;
		smp->Tolerance = 1.0e-10;
		if (solve_PNP(smp)) {
		    printf("NetCurrent = %+.9e\n", smp->NetCurrent);
		    printf("Fluxes = %+.9e, %+.9e\n",
			smp->FLUX[0], smp->FLUX[1]);
		} else {
		    printf("No convergence!\n");
		}
		break;
/* it is safe to insert other cases here */
	      case 21:
/* fixing ions: one type of charge +1, one type of charge -1 */
		nions = 2;
		ngrid = 1000;
		fprintf(stderr, "enter pore radius in angstroms: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &poreradius) != 1) {
		    fprintf(stderr, "error: bad pore radius -- aborting\n");
		    break;
		}
		fprintf(stderr, "enter pore length in angstroms: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &porelength) != 1) {
		    fprintf(stderr, "error: bad pore length -- aborting\n");
		    break;
		}
		smp = build_cyl_geom(nions, ngrid, poreradius, porelength);
		smp->Z[0] = +1.0;
		smp->Z[1] = -1.0;
		fprintf(stderr, "enter left concentration of ions, in M/l: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad concentration -- aborting\n");
		    break;
		}
		smp->CL[0] = temp;
		smp->CL[1] = temp;
		fprintf(stderr, "enter right concentration of ions, in M/l: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad concentration -- aborting\n");
		    break;
		}
		smp->CR[0] = temp;
		smp->CR[1] = temp;
		fprintf(stderr, "enter diffusion coefficient of positive ion, in cm^2/sec: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad coefficient -- aborting\n");
		    break;
		}
		for (k = 0; k < smp->Ngrid; ++k) {
		    smp->Ds[0][k] = temp;
		}
		fprintf(stderr, "enter diffusion coefficient of negative ion, in cm^2/sec: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad coefficient -- aborting\n");
		    break;
		}
		for (k = 0; k < smp->Ngrid; ++k) {
		    smp->Ds[1][k] = temp;
		}
		fprintf(stderr, "enter ambient temperature, in degrees Celsius: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad temperature -- aborting\n");
		    break;
		}
		smp->Celsius = temp;
		fprintf(stderr, "enter relative dielectric permittivity inside pore: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad permittivity -- aborting\n");
		    break;
		}
		for (k = 0; k < smp->Ngrid; ++k) {
		    smp->EPS[k] = temp;
		}
		fprintf(stderr, "\tdriver routine handles linear variations in fixed charge density\nenter fixed charge at left end of pore, in M/l: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &donnan_l) != 1) {
		    fprintf(stderr, "error: bad charge -- aborting\n");
		    break;
		}
		fprintf(stderr, "enter fixed charge at right end of pore, in M/l: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &donnan_r) != 1) {
		    fprintf(stderr, "error: bad charge -- aborting\n");
		    break;
		}
		kl = smp->Ngrid - 1;
		for (k = 0; k <= kl; ++k) {
		    smp->F[k] = (donnan_l*(kl - k) + donnan_r*k)/kl;
		}
		donnan_l = log((donnan_l + sqrt(donnan_l*donnan_l +
		    4.0*smp->CL[0]*smp->CL[1]))/(2.0*smp->CL[1]));
		donnan_r = log((donnan_r + sqrt(donnan_r*donnan_r +
		    4.0*smp->CR[0]*smp->CR[1]))/(2.0*smp->CR[1]));
		for (k = 0; k < nions; ++k) {
		    smp->CL[k] *= exp(-smp->Z[k]*donnan_l);
		    smp->CR[k] *= exp(-smp->Z[k]*donnan_r);
		}
/*the following two lines need to be put in so that we get the correct scale*/
		smp->kT = Boltzmann*(Celsius0 + smp->Celsius);
		smp->kT_e = smp->kT/e0;
		donnan_l *= smp->kT_e;
		donnan_r *= smp->kT_e;
		fprintf(stderr, "\tDonnan potentials in volts:\nleft\t%.6e\nright\t%.6e\n",
		    donnan_l, donnan_r);
		fprintf(stderr, "\tsetting basal free energies to zero\n");
		for (k = 0; k < smp->Ngrid; ++k) {
		    smp->E0s[0][k] = 0.0;
		    smp->E0s[1][k] = 0.0;
		}
		fprintf(stderr, "enter initial transmembrane voltage, in V: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &v_lo) != 1) {
		    fprintf(stderr, "error: bad voltage -- aborting\n");
		    break;
		}
		fprintf(stderr, "enter final transmembrane voltage, in V: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &v_hi) != 1) {
		    fprintf(stderr, "error: bad voltage -- aborting\n");
		    break;
		}
		fprintf(stderr, "enter transmembrane voltage step, in V: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &v_step) != 1) {
		    fprintf(stderr, "error: bad voltage -- aborting\n");
		    break;
		}
		if (v_lo > v_hi) {
		    temp = v_lo;
		    v_lo = v_hi;
		    v_hi = temp;
		}
		if (v_step < 0.0) {
		    v_step = -v_step;
		}
		if (v_step == 0.0 && v_lo != v_hi) {
		    fprintf(stderr, "error: zero voltage step, but nonzero voltage range -- resetting\n");
		    v_step = 0.01*(v_hi - v_lo);
		}
		temp = v_lo;
		smp->MaxIter = 50;
		smp->Tolerance = 1.0e-10;
		fprintf(stderr, "\tsolver verboseness turned off temporarily\n");
		smp->verbose = 0;
		fprintf(stderr, "outputs: transmembrane voltage in mV, pore current in pA, fluxes in ions/sec\n");
		while (1) {
		    smp->VM = temp + donnan_l - donnan_r;
		    printf("%+.5e", mV*temp);
		    if (solve_PNP(smp)) {
			printf(" %+.9e %+.9e %+.9e\n", smp->NetCurrent,
			    smp->FLUX[0], smp->FLUX[1]);
		    } else {
			printf(" -> No convergence!\n");
			break;
		    }
		    if (temp >= v_hi) {
			break;
		    }
		    temp += v_step;
		}
		smp->verbose = 1;
		break;
	      case 31:
/* fixing ions: one type of charge +1, one type of charge -1 */
		nions = 2;
		ngrid = 2500;
		fprintf(stderr, "enter pore radius in angstroms: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &poreradius) != 1) {
		    fprintf(stderr, "error: bad pore radius -- aborting\n");
		    break;
		}
		fprintf(stderr, "enter pore length in angstroms: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &porelength) != 1) {
		    fprintf(stderr, "error: bad pore length -- aborting\n");
		    break;
		}
		fprintf(stderr, "enter atrium length in angstroms: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad atrium length -- aborting\n");
		    break;
		}
		lengths[0] = 100.0;
		lengths[1] = temp;
		lengths[2] = porelength;
		lengths[3] = temp;
		lengths[4] = 100.0;
		smp = build_cab_geom(nions, ngrid, poreradius, lengths,
		    constant_angle, constant_angle, indices, dimensions);
		kf = indices[2];
		kl = kf + dimensions[2];
		smp->Z[0] = +1.0;
		smp->Z[1] = -1.0;
		fprintf(stderr, "enter left concentration of ions, in M/l: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad concentration -- aborting\n");
		    break;
		}
		smp->CL[0] = temp;
		smp->CL[1] = temp;
		fprintf(stderr, "enter right concentration of ions, in M/l: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad concentration -- aborting\n");
		    break;
		}
		smp->CR[0] = temp;
		smp->CR[1] = temp;
		fprintf(stderr, "enter bath diffusion coefficient of positive ion, in cm^2/sec: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad coefficient -- aborting\n");
		    break;
		}
		for (k = 0; k < smp->Ngrid; ++k) {
		    smp->Ds[0][k] = temp;
		}
		fprintf(stderr, "enter bath diffusion coefficient of negative ion, in cm^2/sec: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad coefficient -- aborting\n");
		    break;
		}
		for (k = 0; k < smp->Ngrid; ++k) {
		    smp->Ds[1][k] = temp;
		}

		fprintf(stderr, "enter pore diffusion coefficient of positive ion, in cm^2/sec: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad coefficient -- aborting\n");
		    break;
		}
		for (k = kf; k < kl; ++k) {
		    smp->Ds[0][k] = temp;
		}
		fprintf(stderr, "enter pore diffusion coefficient of negative ion, in cm^2/sec: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad coefficient -- aborting\n");
		    break;
		}
		for (k = kf; k < kl; ++k) {
		    smp->Ds[1][k] = temp;
		}
/* ok to here */
		fprintf(stderr, "enter ambient temperature, in degrees Celsius: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &temp) != 1) {
		    fprintf(stderr, "error: bad temperature -- aborting\n");
		    break;
		}
		smp->Celsius = temp;
		for (k = 0; k < smp->Ngrid; ++k) {
		    smp->EPS[k] = 80.0;
		}
		fprintf(stderr, "\tdriver routine handles linear variations in fixed charge density\nenter fixed charge at left end of pore, in M/l: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &donnan_l) != 1) {
		    fprintf(stderr, "error: bad charge -- aborting\n");
		    break;
		}
		fprintf(stderr, "enter fixed charge at right end of pore, in M/l: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &donnan_r) != 1) {
		    fprintf(stderr, "error: bad charge -- aborting\n");
		    break;
		}
		kl = dimensions[2] - 1;
		for (k = 0; k <= kl; ++k) {
		    smp->F[k+kf] = (donnan_l*(kl - k) + donnan_r*k)/kl;
		}
		fprintf(stderr, "\tsetting basal free energies to zero\n");
		for (k = 0; k < smp->Ngrid; ++k) {
		    smp->E0s[0][k] = 0.0;
		    smp->E0s[1][k] = 0.0;
		}
		fprintf(stderr, "enter initial transmembrane voltage, in V: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &v_lo) != 1) {
		    fprintf(stderr, "error: bad voltage -- aborting\n");
		    break;
		}
		fprintf(stderr, "enter final transmembrane voltage, in V: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &v_hi) != 1) {
		    fprintf(stderr, "error: bad voltage -- aborting\n");
		    break;
		}
		fprintf(stderr, "enter transmembrane voltage step, in V: ");
		fgets(line, 100, stdin);
		if (sscanf(line, "%lf", &v_step) != 1) {
		    fprintf(stderr, "error: bad voltage -- aborting\n");
		    break;
		}
		if (v_lo > v_hi) {
		    temp = v_lo;
		    v_lo = v_hi;
		    v_hi = temp;
		}
		if (v_step < 0.0) {
		    v_step = -v_step;
		}
		if (v_step == 0.0 && v_lo != v_hi) {
		    fprintf(stderr, "error: zero voltage step, but nonzero voltage range -- resetting\n");
		    v_step = 0.01*(v_hi - v_lo);
		}
		temp = v_lo;
		smp->MaxIter = 50;
		smp->Tolerance = 1.0e-10;
		fprintf(stderr, "\tsolver verboseness turned off temporarily\n");
		smp->verbose = 0;
		fprintf(stderr, "outputs: transmembrane voltage in mV, pore current in pA, fluxes in ions/sec\n");
		while (1) {
		    smp->VM = temp;
		    printf("%+.5e", mV*temp);
		    if (solve_PNP(smp)) {
			printf(" %+.9e %+.9e %+.9e\n", smp->NetCurrent,
			    smp->FLUX[0], smp->FLUX[1]);
		    } else {
			printf(" -> No convergence!\n");
			break;
		    }
		    if (temp >= v_hi) {
			break;
		    }
		    temp += v_step;
		}
		smp->verbose = 1;
		break;
/* it is safe to insert other cases here */
	      default:
		fprintf(stderr, "sorry, bad choice\n");
		break;
/* it is safe to insert other cases here */
	    }
	} else {
	    fprintf(stderr, "error: bad input\n");
	}
    }
  done_loop:
    return(0);
}
