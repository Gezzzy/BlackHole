#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Odes.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double G = 1.0;
double M = 1.0;
double E;   // conserved energy

double* rightSide(int N, double lambda, double *y)
{
    double *yt = malloc(sizeof(double)*N);

    double r      = y[0];
    double theta  = y[1];
    double phi    = y[2];
    double vr     = y[3];
    double vtheta = y[4];
    double vphi   = y[5];

    double f = 1.0 - (2.0*G*M)/r;

    // Reconstruct ṫ from conserved energy
    double tdot = E / f;

    // Equations of motion

    // dr/dλ
    yt[0] = vr;

    // dθ/dλ
    yt[1] = vtheta;

    // dφ/dλ
    yt[2] = vphi;

    // dvr/dλ (radial acceleration)
    yt[3] =
        - (G*M)/(r*r*r) * (r - 2*G*M) * tdot*tdot
        + (G*M)/(r*(r - 2*G*M)) * vr*vr
        + (r - 2*G*M)*(vtheta*vtheta + sin(theta)*sin(theta)*vphi*vphi);

    // d(vθ)/dλ
    yt[4] =
        - (2.0/r) * vtheta * vr
        + sin(theta)*cos(theta)*vphi*vphi;

    // d(vφ)/dλ
    yt[5] =
        - (2.0/r) * vphi * vr
        - 2.0 * (cos(theta)/sin(theta)) * vtheta * vphi;

    return yt;
}

int main()
{
    int n = 20000;
    int N = 6;

    double x0 = 0.0;
    double xf = 1600.0;

    double r0 = 10.0;
    double theta0 = M_PI / 3.0;
    double phi0 = 0.0;

    double vr0 = 0.0;
    double vtheta0 = 0.0;
    double vphi0 = 0.05;

    double f0 = 1.0 - (2.0*G*M)/r0;

    E = sqrt(
        f0 * (1.0
        + (vr0*vr0)/f0
        + r0*r0*(vtheta0*vtheta0 + sin(theta0)*sin(theta0)*vphi0*vphi0))
    );

    double y0V[6] = { r0, theta0, phi0, vr0, vtheta0, vphi0 };

    double *solution = SolveOdeSystemRK4(n, N, x0, y0V, xf, rightSide);

    for(int i = 0; i < n; i++)
    {
        double r      = solution[0*n + i];
        double theta  = solution[1*n + i];
        double phi    = solution[2*n + i];
        double vr     = solution[3*n + i];
        double vtheta = solution[4*n + i];
        double vphi   = solution[5*n + i];

        double lambda = solution[N*n + i];

        printf("%f:%f:%f:%f:%f:%f:%f\n", lambda, r, theta, phi, vr, vtheta, vphi);
    }

    return 0;
}
