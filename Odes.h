#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

double complex* SolveOdeSystemStepMethodComplex(int n, int N, double x0, double complex y0[n], double xf, double complex *(*Phi) (int, double, double complex*, double, double complex * (*ff) (int, double, double complex*)), double complex * (*f) (int, double, double complex*))
{
    /*
     * Implementa un generico Step method utilizzando la funzione
     * \Phi(x, y, h, f) passata come parametro
     *
     */

    // n numero di step, N dimensione vettore
    double h = (xf - x0)/n;

    double complex *solution = malloc(sizeof(double complex)*n*(N+1));
    // y_1 per n valori, y_2 per n valori e cosi via

    double complex yLast[N];    //contiene i componenti dello step i - 1
    for(int i = 0; i < N; i++)
    {
        yLast[i] = y0[i];
        solution[n*i] = yLast[i];
    }
    solution[N*n] = x0;

    for(int i = 1; i < n; i++)
    {
        double complex * phi = Phi(N, x0 + i*h, yLast, h, f);

        //printf("%d/%d\n", i, n);
        // Calcola step successivo
        for(int j = 0; j < N; j++)
            solution[j*n + i] = yLast[j] + h*phi[j];
        //printf("Calcolato\n");

        // Salva lo step attuale
        for(int j = 0; j < N; j++)
            yLast[j] = solution[j*n + i];

        // Salva i vari xi
        solution[N*n + i] = x0 + i*h;
    }

    return solution;
}

double* SolveOdeSystemStepMethod(int n, int N, double x0, double y0[n], double xf, double *(*Phi) (int, double, double*, double, double * (*ff) (int, double, double*)), double * (*f) (int, double, double*))
{
    /*
     * Implementa un generico Step method utilizzando la funzione
     * \Phi(x, y, h, f) passata come parametro
     *
     */

    // n numero di step, N dimensione vettore
    double h = (xf - x0)/n;

    double *solution = malloc(sizeof(double)*n*(N+1));
    // y_1 per n valori, y_2 per n valori e cosi via

    double yLast[N];    //contiene i componenti dello step i - 1
    for(int i = 0; i < N; i++)
    {
        yLast[i] = y0[i];
        solution[n*i] = yLast[i];
    }
    solution[N*n] = x0;

    for(int i = 1; i < n; i++)
    {
        double * phi = Phi(N, x0 + i*h, yLast, h, f);

        // Calcola step successivo
        for(int j = 0; j < N; j++)
            solution[j*n + i] = yLast[j] + h*phi[j];

        // Salva lo step attuale
        for(int j = 0; j < N; j++)
            yLast[j] = solution[j*n + i];

        // Salva i vari xi
        solution[N*n + i] = x0 + i*h;
    }

    return solution;
}


double* SolveOdeSystemEuler(int n, int N, double x0, double y0[n], double xf, double *(*f) (int, double, double*))
{
    double *Phi(int N, double x, double *y, double h, double* (*f)(int, double, double*))
    {
        return (*f)(N, x, y);
    }

    double *solution = SolveOdeSystemStepMethod(n, N, x0, y0, xf, Phi, f);
    return solution;
}

double* SolveOdeSystemRK2(int n, int N, double x0, double y0[n], double xf, double *(*f) (int, double, double*))
{
    double *Phi(int N, double x, double *y, double h, double* (*f)(int, double, double*))
    {
        // Devo calcolare le nuove y => y + 0.5*h*(*f)(h, x, y)
        double *funcValue = (*f)(N, x, y);

        double *newY = malloc(sizeof(double) * N);
        for(int i = 0; i < N; i++)
            newY[i] = y[i] + 0.5*h*funcValue[i];

        return (*f)(N, x + 0.5*h, newY);
    }

    double *solution = SolveOdeSystemStepMethod(n, N, x0, y0, xf, Phi, f);
    return solution;
}

double complex* SolveOdeSystemRK2Complex(int n, int N, double x0, double complex y0[n], double xf, double complex *(*f) (int, double, double complex*))
{
    double complex *Phi(int N, double x, double complex *y, double h, double complex* (*f)(int, double, double complex*))
    {
        // Devo calcolare le nuove y => y + 0.5*h*(*f)(h, x, y)
        double complex *funcValue = (*f)(N, x, y);
        double complex *newY = malloc(sizeof(double complex) * N);
        for(int i = 0; i < N; i++)
            newY[i] = y[i] + 0.5*h*funcValue[i];

        double complex *ret = (*f)(N, x + 0.5*h, newY);

        free(newY);
        free(funcValue);

        return ret;
    }

    double complex *solution = SolveOdeSystemStepMethodComplex(n, N, x0, y0, xf, Phi, f);
    return solution;
}

double* SolveOdeSystemRK4(int n, int N, double x0, double y0[n], double xf, double *(*f) (int, double, double*))
{
    double *Phi(int N, double x, double *y, double h, double* (*f)(int, double, double*))
    {
        double *k1 = (*f)(N, x, y);

        double *y1 = malloc(sizeof(double)*N);

        for(int i = 0; i < N; i++)
            y1[i] = y[i] + 0.5*h*k1[i];

        double *k2 = (*f)(N, x + h/2, y1);

        for(int i = 0; i < N; i++)
            y1[i] = y[i] + 0.5*h*k2[i];

        double *k3 = (*f)(N, x + 0.5*h, y1);

        for(int i = 0; i < N; i++)
            y1[i] = y[i] + h*k3[i];

        double *k4 = (*f)(N, x + 0.5*h, y1);

        double *ret = malloc(sizeof(double)*N);
        for(int i = 0; i < N; i++)
            ret[i] = (k1[i]  + 2*k2[i]  + 2*k3[i] + k4[i])/6;

        return ret;
    }

    double *solution = SolveOdeSystemStepMethod(n, N, x0, y0, xf, Phi, f);
    return solution;
}


double complex* SolveOdeSystemRK4Complex(int n, int N, double x0, double complex y0[n], double xf, double complex *(*f) (int, double, double complex*))
{
    double complex *Phi(int N, double x, double complex *y, double h, double complex* (*f)(int, double, complex double*))
    {
        double complex *k1 = (*f)(N, x, y);

        double complex *y1 = malloc(sizeof(double complex)*N);

        for(int i = 0; i < N; i++)
            y1[i] = y[i] + 0.5*h*k1[i];

        double complex *k2 = (*f)(N, x + h/2, y1);

        for(int i = 0; i < N; i++)
            y1[i] = y[i] + 0.5*h*k2[i];

        double complex *k3 = (*f)(N, x + 0.5*h, y1);

        for(int i = 0; i < N; i++)
            y1[i] = y[i] + h*k3[i];

        double complex *k4 = (*f)(N, x + 0.5*h, y1);

        double complex *ret = malloc(sizeof(double complex)*N);
        for(int i = 0; i < N; i++)
            ret[i] = (k1[i]  + 2*k2[i]  + 2*k3[i] + k4[i])/6;

        free(k1);
        free(k2);
        free(k3);
        free(k4);
        free(y1);

        return ret;
    }

    double complex *solution = SolveOdeSystemStepMethodComplex(n, N, x0, y0, xf, Phi, f);
    return solution;
}

double* SolveOdeSystemRK4SBAGLIATO(int n, int N, double x0, double y0[n], double xf, double *(*f) (int, double, double*))
{
    double *Phi(int N, double x, double *y, double h, double* (*f)(int, double, double*))
    {
        double *k1 = (*f)(N, x, y);

        double *y1 = malloc(sizeof(double)*N);

        for(int i = 0; i < N; i++)
            y1[i] = y[i] + 0.5*h*k1[i];

        double *k2 = (*f)(N, x + h/2, y1);

        for(int i = 0; i < N; i++)
            y1[i] = y[i] + 0.5*h*k2[i];

        double *k3 = (*f)(N, x + 0.5*h, y1);

        for(int i = 0; i < N; i++)
            y1[i] = y[i] + h*k2[i]; //errore qui

        double *k4 = (*f)(N, x + 0.5*h, y1);

        double *ret = malloc(sizeof(double)*N);
        for(int i = 0; i < N; i++)
            ret[i] = (k1[i]  + 2*k2[i]  + 2*k3[i] + k4[i])/6;

        return ret;
    }

    double *solution = SolveOdeSystemStepMethod(n, N, x0, y0, xf, Phi, f);
    return solution;
}


double *SolveSecondOrderOdeRK2(double x0, double xf, double y0, double yd0, int n, double (*f)(double, double, double))
{
    int N = 2;
    double *func(int N, double x, double y[N])
    {
        double *ret = malloc(sizeof(double)*N);
        ret[0] = y[1];
        ret[1] = f(x, y[0], y[1]);
        return ret;
    }
    double y0V[2] = { y0, yd0 };

    double *solution = SolveOdeSystemRK2(n, N, x0, y0V, xf, func);
    return solution;
}


double *SolveSecondOrderOdeRK4(double x0, double xf, double y0, double yd0, int n, double (*f)(double, double, double))
{
    int N = 2;
    double *func(int N, double x, double y[N])
    {
        double *ret = malloc(sizeof(double)*N);
        ret[0] = y[1];
        ret[1] = f(x, y[0], y[1]);
        return ret;
    }
    double y0V[2] = { y0, yd0 };

    double *solution = SolveOdeSystemRK4(n, N, x0, y0V, xf, func);
    return solution;
}

double *SolveSecondOrderOdeRK4SBAGLIATO(double x0, double xf, double y0, double yd0, int n, double (*f)(double, double, double))
{
    int N = 2;
    double *func(int N, double x, double y[N])
    {
        double *ret = malloc(sizeof(double)*N);
        ret[0] = y[1];
        ret[1] = f(x, y[0], y[1]);
        return ret;
    }
    double y0V[2] = { y0, yd0 };

    double *solution = SolveOdeSystemRK4SBAGLIATO(n, N, x0, y0V, xf, func);
    return solution;
}

double *SolveSecondOrderOdeEuler(double x0, double xf, double y0, double yd0, int n, double (*f)(double, double, double))
{
    int N = 2;
    double *func(int N, double x, double y[N])
    {
        double *ret = malloc(sizeof(double)*N);
        ret[0] = y[1];
        ret[1] = f(x, y[0], y[1]);
        return ret;
    }
    double y0V[2] = { y0, yd0 };

    double *solution = SolveOdeSystemEuler(n, N, x0, y0V, xf, func);
    return solution;
}
