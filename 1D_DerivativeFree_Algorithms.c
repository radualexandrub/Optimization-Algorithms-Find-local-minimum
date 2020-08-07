#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define TOL 1e-3
#define PI 3.141592

FILE *logfile_grid;
FILE *logfile_dichotomy;
FILE *logfile_golden_ratio;
FILE *logfile_interpolare_parab;

double func1_sin(double x);
double func2_pol(double x);
double func3_exp(double x);
double func4_exp(double x);

double func1_sin(double x) { return sin(x); }
double func2_pol(double x) { return 3*pow(x,2)+2*x-3; }
double func3_exp(double x) { return -exp(pow((x-2),2)/2); }
double func4_exp(double x) { return -x*exp(pow(x,2)/2); }

double grid_search(double f(double), double start, double end, double tol);
double dichotomy(double f(double), double a, double b, double tol);
double golden_ratio(double f(double), double a, double b, double tol);
double interpolare_parab(double f(double), double a, double b, double tol);

int main() {
    logfile_grid = fopen("log_grid.csv","w");
    logfile_dichotomy = fopen("log_dichotomy_sin.csv","w");
    logfile_golden_ratio = fopen("log_raportaur_sin.csv", "w");
    logfile_interpolare_parab = fopen("log_interpolare_exp3.csv", "w");

    if(logfile_grid == NULL ) return 1;

    double op = grid_search(func1_sin, 0, 2*PI, 0.05);
    printf("Optimum_grid: x=%f, f(x)=%f \n", op, func1_sin(op));

    double op_dichotomy = dichotomy(func1_sin, 0, 2*PI, TOL);
    printf("Optimum_dihot: %f, %f \n", op_dichotomy, func1_sin(op_dichotomy));

    double op_raportaur = golden_ratio(func1_sin, 0, 2*PI, TOL);
    printf("Optimum_golden_ratio: %f, %f \n", op_raportaur, func1_sin(op_raportaur));

    double op_interpolare = interpolare_parab(func1_sin, 0, 2*PI, TOL);
    printf("Opt_Interpolare_parab: %f, %f \n", op_interpolare, func1_sin(op_interpolare));

    return 0;
}


/* ***********************
****** ALGORITMII ********
*********************** */
double grid_search(double f(double), double start, double end, double tol){
    double optim = f(start);
    double x, arg_min;
    
    int pasi = 0;

    for(x=start; x < end; x+=tol) {
        if(f(x) < optim) {
            optim = f(x);
            arg_min = x;
        }
        fprintf(logfile_grid,"%f, %f, %f, %f \n", x, f(x), arg_min, f(arg_min));
        // optim = f(arg_min);

        pasi += 1;
    }
    printf("%d iteratii\n", pasi);
    return arg_min;
}

double dichotomy(double f(double), double a, double b, double tol) {
    double c, d, Eps;
    Eps = tol/4;
    while(fabs(a-b) > tol) {
        c = (a+b)/2 - Eps;
        d = (a+b)/2 + Eps;

        if (f(c) > f(d)) {
            a = c;
        } else {
            b = d;
        }
        fprintf(logfile_dichotomy,"%f, %f, %f, %f \n", a, f(a),b,f(b));
     }
    return (a+b)/2;
}


double golden_ratio(double f(double), double a, double b, double tol) {
    double R = (1+sqrt(5))/2;
    double sect_c = a + (b-a)/R;
    double sect_d = b - (b-a)/R;

    while(fabs(a-b)/2 > tol) {

        if (f(sect_c) < f(sect_d)) {
            a = sect_d;
            sect_d = sect_c;
            sect_c = a + (b-a)/R;
        } else {
            b = sect_c;
            sect_c = sect_d;
            sect_d = b - (b-a)/R;
        }
        fprintf(logfile_golden_ratio,"%f, %f, %f, %f \n", a, f(a),b,f(b));
     }
    return (a+b)/2;
}


double sectionare(double f(double), double a, double c, double b) {
    double x1 = a;
    double x2 = c;
    double x3 = b;

    double denom = (x1 - x2)*(x1 - x3)*(x2 - x3);
    double d21 = f(x2) - f(x1);
    double d12 = f(x1) - f(x2);
    double d13 = f(x1) - f(x3);
    double d31 = f(x3) - f(x1);
    double d32 = f(x3) - f(x2);
    double d23 = f(x2) - f(x3);

    double A = (x3*d21 + x2*d13 + x1*d32) / denom;
    double B = (pow(x3,2)*d12 + pow(x2,2)*d31 + pow(x1,2)*d23)/denom;

    double vertex = -B/(2*A);
    return vertex;
}

double interpolare_parab(double f(double), double a, double b, double tol) {
    double R = (1+sqrt(5))/2;
    double c, d;

    while(fabs(a-b)/2 > tol) {
        c = a + (b-a)/R;
        d = sectionare(f,a,c,b);

        if (d>c) {
            float tmp = c;
            c = d;
            d = tmp;
        }
        //printf("x");
        if (f(c) < f(d)) {
            a = d;
        } else {
            b = c;
        }
        fprintf(logfile_interpolare_parab,"%f, %f, %f, %f \n", a, f(a),b,f(b));
     }
    return (a+b)/2;
}
