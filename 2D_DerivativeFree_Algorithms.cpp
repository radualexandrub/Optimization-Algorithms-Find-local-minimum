#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

using namespace std;

FILE *logfile1;
FILE *logfile2;
FILE *logfile3;

/** Functions declarations **/
void Cautare_directa(double f(double, double), double alfa, double tol, double x_now[]);
void Cautare_coord(double f(double,double), double pas, double alfa, double eps, double x[]);

/** Test functions **/
double f1_paraboloid(double x, double y) {
    return x*x + y*y - x*y - x - y; }

double f2_paraboloid(double x, double y) {
    return x*x + y*y - x*y - 2*x - 2*y; }

double f3_rosenbrock(double x, double y) {
    return pow((1-x), 2) + 100*pow((y-x*x), 2); }

/************************
** Algorithms / Methods **
************************/
void Cautare_directa(double f(double, double),
                     double alfa,
                     double tol,
                     double x_now[]) {
    int i, j, k;
    int n_dim = 2, n_dir = 3;
    double step[n_dim];
    double dir[n_dir][n_dim] = {{1, 0}, {-0.5, sqrt(3)/2}, {-0.5, -sqrt(3)/2}};     // ne declaram o matrice pt directii
    double x_next[n_dim];

    printf("Coord start: %.3f, %.3f\n", x_now[0], x_now[1]);
    printf("\tf(x1,x2)=%.3f\n", f(x_now[0], x_now[1]));
    fprintf(logfile1, "x1,x2,f(x1 x2)\n");
    fprintf(logfile1, "%f, %f, %f\n", x_now[0],x_now[1], f(x_now[0], x_now[1]));

    int nr_iter = 10;
    for (j = 0; j < nr_iter; j++) {
        for (i = 0; i < n_dir; i++) { // Parcurgem cele 3 directii opuse la 120 grade

            /* a_scale: Inmultire array cu scalar */
            for (k = 0; k < n_dim; k++)
                step[k] = alfa*dir[i][k];

            /* a_add: Adunare array cu array */
            for (k = 0; k < n_dim; k++)
                x_next[k] = x_now[k] + step[k];

            if (i==1) // Afisare pe o singura directie
            {
            printf("i:%d Directia %d:\n\t%.3f, %.3f\n\tf(x1,x2)=%.3f\n\tPAS: %.9f %.9f\n",
                   j, i,
                   x_next[0], x_next[1],
                   f(x_next[0], x_next[1]),
                   step[0], step[1]);
            }
            fprintf(logfile1, "%f, %f, %f\n", x_next[0], x_next[1], f(x_next[0], x_next[1]));

            if (f(x_next[0], x_next[1]) < f(x_now[0], x_now[1])){
                /* Atribuire array */
                for (int k = 0; k < n_dim; k++)
                    x_now[k] = x_next[k];
            } else {
                 alfa = alfa/1.5;        // Paraboloid
                 //alfa = alfa/1.01;   // Rosenbrock
            }
        }
        if (alfa < tol) {
            break;
        }
    }
}

void Cautare_coord(double f(double,double),
                    double pas,
                    double alfa,
                    double eps,
                    double x[]) {
    printf("Coord start: %.3f, %.3f\n", x[0], x[1]);
    fprintf(logfile2, "x1,x2,f(x1 x2)\n");
    fprintf(logfile2, "%f, %f, %f\n", x[0],x[1], f(x[0], x[1]));

    int directie, nr_iter = 10;
    double f0;
    for (int j = 0; j < nr_iter; j++) {
        // Pentru fiecare coordonata folosim un algoritm de optimizare 1D: Line search
        for(int i = 0; i < 2; i ++) {
            printf("j = %d\n\tx = (%.3f, %.3f)\n\tf(x1,x2)=%.3f\tPAS: %.6f\n", j,
                   x[0], x[1],
                   f(x[0], x[1]),
                   pas);
            fprintf(logfile2, "%f, %f, %f\n", x[0], x[1], f(x[0], x[1]));

            f0 = f(x[0], x[1]);

            directie = 1;
            x[i]= x[i] + pas*directie;

            if (f(x[0], x[1]) > f0) {
                directie = -1;
                x[i] = x[i] + 2*pas*directie;
            }

            while(f0 >= f(x[0], x[1])) {
                f0 = f(x[0], x[1]);
                x[i] = x[i] + pas*directie;
            }
        }
        pas=pas*alfa;
        if (pas < eps) {
            break;
        }
    }
}


/**************
***** MAIN ****
**************/
int main() {
    logfile1 = fopen("log1_direct.csv","w");
    logfile2 = fopen("log2_coord.csv","w");
    if(logfile1==NULL || logfile2==NULL) {
        cout << "Could not open one of files sorry.\n";
        return 1;
    }

    printf("*** CAUTARE DIRECTA ***\n");
    double alfa = 2;
    double tol = 1e-4;
    double x_start[2] = {0, 2}; //punct de pornire
    // Cautare_directa(f1_paraboloid, alfa, tol, x_start);

    printf("\n*** CAUTARE PE COORDONATE ***\n");
    double pas = 0.1;
    tol = 1e-8;
    x_start[0] = 0;
    x_start[1] = 0;
    alfa = 0.5;
    Cautare_coord(f2_paraboloid, pas, alfa, tol, x_start);

    return 0;
}
