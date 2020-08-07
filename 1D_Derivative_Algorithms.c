#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define TOL 0.0005

FILE *logfile;
FILE *logfile2;
FILE *logfile3;
FILE *logfile4;

double func1(double x);
double func1D(double x);
double func1DD(double x);

double func2_sin(double x);
double func2_sinD(double x);
double func2_sinDD(double x);

double func1(double x) {
    return -x*exp(-x*x); }
double func1D(double x) {
    return exp(-x*x)*(2*x*x-1); }
double func1DD(double x) {
    return exp(-x*x)*(6*x-4*x*x*x);}

double func2_sin(double x) { return sin(x); }
double func2_sinD(double x) { return cos(x); }
double func2_sinDD(double x) { return -sin(x); }

double gradientDescent(double f(double), double fd(double), double start, double alfa, int nr_iteratii);
double newton(double f(double), double fD(double),double fDD(double), double start, int nr_iteratii);
double secant(double f(double), double fd(double), double fdd(double), double start, double alfa, double nr_iteratii);
double falsi(double f(double), double fd(double), double a, double b, double tol);

double f(double x) {
	return (x - log(x));
}
double fd(double x) {
	return (x-1)/x;
}
double fdd(double x) {
	return 1/(x*x);
}

/**************
***** MAIN ****
**************/
int main() {
    logfile = fopen("log_grad.csv","w");
    logfile2 = fopen("log_newton.csv","w");
    logfile3 = fopen("log_secant.csv","w");
    logfile4 = fopen("log_falsi.csv","w");

    if(logfile==NULL || logfile2==NULL || logfile3==NULL || logfile4==NULL) {
        printf("Could not open one of files\n");
        return 1;
    }

    // Parametri ale functillor
    float start = 0; //0.1... 20...
    float alfa = 0.3; // 0.3... 1...
    int nr_iteratii_max = 20;
    ///////////////////////////////

    printf("*** GRADIENT DESCENDENT ***\n");
    // double minim = gradientDescent(func1, func1D, start, alfa, nr_iteratii_max);
    printf("\n*** METODA NEWTON ***\n");
    // double minim2 = newton(f, fd, fdd, alfa, nr_iteratii_max);

    printf("\n*** METODA SECANTEI ***\n");
    double minim3 = secant(func2_sin, func2_sinD, func2_sinDD, start, alfa, nr_iteratii_max);

    printf("\n*** REGULA FALSI ***\n");
    float a = -2;
    float b = 2;
    // double minim4 = falsi(func1, func1D, a, b, TOL);
}

/*****************
**** FUNCTIONS ***
******************/
double gradientDescent(double f(double),
                          double fd(double),
                          double start, double alfa, int nr_iteratii) {
    double x = start;
    double y;
    for(int i = 0; i < nr_iteratii; i++) {
        x = x - alfa * fd(x);
        if (fabs(alfa * fd(x)) < 0.00001)
            break;

        printf("x=%f\n", x);
        printf("alphaf=%f\n", alfa * fd(x));

        fprintf(logfile, "%f, %f\n", x, f(x));
    }
    return x;
}

double newton(double f(double),
              double fD(double),
              double fDD(double),
              double start, int nr_iteratii){
    double x = start;
    double Eps = 1e-7;
    double pas = 1/fabs(fDD(x));

    fprintf(logfile2,"x, f(x)\n");
    for(int i = 0; i < nr_iteratii; i++) {
        x = x - fD(x)/fabs(fDD(x));
        if(pas < Eps)
            break;

        printf("x=%f, f(x)=%f \n",x, f(x));
        fprintf(logfile2,"%f, %f\n", x, f(x));
    }
    return x;
}

double secant(double f(double),
                double fd(double),
                double fdd(double),
                double start, double alfa, double nr_iteratii){
    double aux;
    double xn_1 = start;
    double x = xn_1-alfa*fd(x);
    double Eps = 1e-7;

    fprintf(logfile3,"x, f(x)\n");

    for(int i=0; i<nr_iteratii; i++){
        if (fabs(fd(x)) < Eps)
            break;

        aux = x;
        //x = x-fd(x)*( (x-xn_1)/(fd(x)-fd(xn_1)) );
        x = xn_1 - (fd(xn_1)/fdd(x));
        xn_1 = aux;

        printf("x=%f, f(x)=%f \n",x, f(x));
        fprintf(logfile3,"%f, %f\n", x, f(x));
    }
    return x;
}

double falsi(double f(double), double fd(double), double a, double b, double tol) {
    fprintf(logfile4,"a, f(a), b, f(b)\n");

    while (fabs(a-b) > tol) { // adica intervalul este peste o toleranta
        double c = a - fd(a)*( (a-b)/(fd(a)-fd(b)) );

        if (c <= a + tol || c > b - tol)
            c = (a + b) / 2;
        if (fd(c) < 0)
            a = c;
        else b = c;
        printf("a=%f, f(a)=%f, b=%f, f(b)=%f\nx_opt=%f\n", a, f(a), b, f(b), (a + b) / 2);
        fprintf(logfile4,"%f, %f, %f, %f\n", a, f(a), b, f(b));
    }
    return ((a+b)/2);
}
