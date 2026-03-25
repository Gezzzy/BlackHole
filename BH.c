#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Odes.h"

int main(int argc, char* argv[])
{
    int n = 16000;

    double M = 1;
    double G = 1;

    double td0 = 0;
    double t0 = 0;

    double thetad0 = 0;
    double theta0 = 3.14/2;

    double rd0 = 0;
    double r0 = 2*G*M*3;

    double phid0 = 1;
    double phi0 = 0;

    int N = 8;
    
    double x0 = 0;
    double xf = 25;
    double y0V[8] = { td0, thetad0, rd0, phid0, t0, theta0, r0, phi0 };

    double f0 = 1.0 - (2.0*G*M)/r0;

    double E0 = f0 * td0;
    double L0 = r0*r0 * sin(theta0)*sin(theta0) * phid0;

    double *rightSide(int N, double x, double *y)
    {
        double *yt = malloc(sizeof(double)*N);

        //Define variables:
        double td = y[0];
        double thetad = y[1];
        double rd = y[2];
        double phid = y[3];
        double t = y[4];
        double theta = y[5];
        double r = y[6];
        double phi = y[7];

        //Define EoM:
        double vt = td;
        double vtheta = thetad;
        double vr = rd;
        double vphi = phid;
        double at = -(2*G*M)/( r * (r - 2*G*M) ) * rd * td;
        double atheta = - (2/r) * thetad * rd + sin(theta) * cos(theta) * phid * phid;
        double ar = - ((G*M)/(r*r*r))*(r - 2*G*M) * td*td + (G*M)/(r*(r - 2*G*M)) * rd* rd + (r - 2*G*M)*(thetad*thetad + sin(theta)*sin(theta)*phid*phid);
        double aphi = - ((2)/(r)) * phid * rd - 2 * (cos(theta)/sin(theta)) * thetad * phid;

        //Redefine yt
        yt[0] = at;
        yt[1] = atheta;
        yt[2] = ar;
        yt[3] = aphi;
        yt[4] = vt;
        yt[5] = vtheta;
        yt[6] = vr;
        yt[7] = vphi;

        return yt;
    }

    double *solution = SolveOdeSystemRK4(n, N, x0, y0V, xf, rightSide);
    for(int i = 0; i < n; i++)
    {
        double vt = solution[0*n + i];
        double vtheta = solution[1*n + i];
        double vr = solution[2*n + i];
        double vphi = solution[3*n + i];
        double t = solution[4*n + i];
        double theta = solution[5*n + i];
        double r = solution[6*n + i];
        double phi = solution[7*n + i];

        double x = r*sin(theta)*cos(phi);
        double y = r*sin(theta)*sin(phi);
        double z = r*cos(theta);

        double lambda = solution[N*n + i];
        printf("%f:%f:%f:%f\n", lambda, x, y, z);
    }
}
