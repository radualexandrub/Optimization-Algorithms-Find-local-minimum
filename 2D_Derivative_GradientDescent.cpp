#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

using namespace std;

FILE *logfile3;

void Gradient_Descent(double f(double,double),double fDx(double, double), double fDy(double, double), double alfa, double tol, double x_now[]);

/** Test Functions **/
double f1_paraboloid(double x, double y) {
    return x*x + y*y - x*y - x - y; }
double f1_paraboloid_Dx(double x, double y) {
    return 2*x - y - 1; }
double f1_paraboloid_Dy(double x, double y) {
    return -x + 2*y - 1; }

double f2_rosenbrock(double x, double y) {
    return pow((1-x), 2) + 100*pow((y-x*x), 2); }
double f2_rosenbrock_Dx(double x, double y) {
    return 2*( 200*x*x*x - 200*x*y + x - 1 ); }
double f2_rosenbrock_Dy(double x, double y) {
    return 200*( y - x*x ); }

void Gradient_Descent(double f(double, double),
    double fDx(double, double),
    double fDy(double, double),
    double alfa,
    double tol,
    double x_now[]) {

    printf("Coord start: %.3f, %.3f\n", x_now[0], x_now[1]);
    fprintf(logfile3, "x1,x2,f(x1 x2)\n");
    fprintf(logfile3, "%f, %f, %f\n", x_now[0],x_now[1], f(x_now[0], x_now[1]));

    double x_new[2];
	double distanta;
	int k, nr_iter = 10;

	for (k = 0; k < nr_iter; k++) {
		x_new[0] = x_now[0] - fDx(x_now[0], x_now[1])*alfa;
		x_new[1] = x_now[1] - fDy(x_now[0], x_now[1])*alfa;

		distanta = sqrt(pow((x_new[0] - x_now[0]), 2) + pow((x_new[1] - x_now[1]), 2));
        
		x_now[0] = x_new[0];
		x_now[1] = x_new[1];

        printf("k=%d  x = (%.3f, %.3f)  f(x1,x2)=%.3f\n", k,
                x_new[0], x_new[1],
                f(x_new[0], x_new[1]));
        fprintf(logfile3, "%f, %f, %f\n", x_new[0], x_new[1], f(x_new[0], x_new[1]));

        if (distanta < tol) {
            break;
        }
	}
}


/**************
***** MAIN ****
**************/
int main() {
    logfile3 = fopen("log3_grad.csv","w");
    if(logfile3==NULL) {
        cout << "Could not open one of files sorry.\n";
        return 1;
    }

    printf("\n*** GRADIENT DESCENT ***\n");
    alfa = 0.5;
    tol = 1e-8;
    x_start[0] = 2;
    x_start[1] = 0;
    // Gradient_Descent(f1_paraboloid, f1_paraboloid_Dx, f1_paraboloid_Dy, alfa, tol, x_start);
    //Gradient_Descent(f2_rosenbrock, f2_rosenbrock_Dx, f2_rosenbrock_Dy, alfa, tol, x_start);

    return 0;
}
