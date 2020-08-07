#include <stdio.h>
#include <stdlib.h>
#include <math.h>

FILE *logfile1, *logfile2, *logfile3;

double f2_rosenbrock(double x, double y) {
	return pow((1-x), 2) + 100*pow((y-x*x), 2); }
double f2_rosenbrock_Dx(double x, double y) {
	return 2*( 200*x*x*x - 200*x*y + x - 1 ); }
double f2_rosenbrock_Dy(double x, double y) {
	return 200*( y - x*x ); }

double f2_paraboloid_xy(double x, double y) {
	return x*x + y*y;
}
double f2_paraboloid_xy_Dx(double x, double y) {
	return 2*x;
}
double f2_paraboloid_xy_Dy(double x, double y) {
	return 2*y;
}

/******************************************
/* Functii ajutatoare de afisare (debug) */
void fprint(FILE *log, double vec[2]) {
	fprintf(log, "%.6f,%.5f\n", vec[0], vec[1]);
}
void print(double x) {
	printf("%.6f\n",x);
}
void printVec(int n, double vec[]) {
	printf("vec = [");
	for (int i = 0; i < n; i++) {
		printf(" %.6f", vec[i]);
	}
	printf(" ]\n");
}
void printMat(int n, int m, double mat[n][m]) {
	printf("mat = [");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			printf(" %.6f", mat[i][j]);
		}
		printf(";");
	}
	printf("]\n");
}

/************************************************
/* Implementarea functiilor de algebra liniara */
void Vec_negat(int n, double vec[n]) {
	for (int i = 0; i < n; i++) {
		vec[i] = -vec[i];
	}
}
void Vec_add_coefVec(int n,
				double coef,
				double a[n],
				double b[n],
				double out[n]) {
	/**a_scaleadd**/
	for (int i = 0; i < n; i++) {
		out[i] = a[i] + coef*b[i];
	}
}
void Vec_copy(int n, double in[n], double out[n]) {
	for (int i = 0; i < n; i++) {
		out[i] = in[i];
	}
}
void Mat_mult_Vec(int n, int m,
					double mat[n][m],
					double vec[m],
					double out[n]) {
	/**mat_vec: Inmultire matrice cu vector "coloana" (pt aflare gradient) **/
	// Initializam out cu zero
	for (int i = 0; i < m; i++) { out[i] = 0; }
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			out[i] += mat[i][j] * vec[j];
		}
	}
}
double Vec_Mat_Vec(int n, int m,
				   double va[n],
				   double mat[n][m],
				   double vb[m]) {
	/** Sum(a[i]*A[i][j]*b[j]) **/
	// Prima parte: va(2,1)*A(2,2)=tmp(1,2)
	double va_mat[n], res = 0;
	for (int i = 0; i < n; i++) { va_mat[i] = 0; }
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			va_mat[i] += va[j]*mat[j][i];
		}
	}
	// A doua parte: tmp(1,2)*vb(2,1)=result(1,1)
	for (int i = 0; i < n; i++) {
		res += va_mat[i] * vb[i];
	}
	return res;
}
double Vec_norm(int n, double vec[n]) {
	/**Norma = Radical din suma patratelor elementelor lui a**/
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += pow(vec[i], 2);
	}
	return sqrt(sum);
}
double Vec_sprod(int n, double a[n], double b[n]) {
	/**Produs scalar (suma din inmultirea punct la punct a elem) intre doi vectori
	<a,b> = a[0]*b[0] + a[1]*b[1];
	necesar pt gradient conjugat**/
	double result;
	for (int i = 0; i < n; i++) {
		result += a[i] * b[i];
	}
	return result;
}

/*******************************
/* Implementarea ALGORITMILOR */
/*******************/
/** METODA CAUCHY **/
void Cauchy(int n,
			double A[n][n], double b[n],
			double tol,
			double x[n]) {
	fprintf(logfile1, "x[0],x[1]\n");

	double tmp[n], rezid[n];
	double lambda_tmp1, lambda_tmp2, lambda;

	int nr_iter_max = 15;
	for (int k = 0; k < nr_iter_max; k++) {

		/* Calcul rezidual r = -(Ax-b) */
		Mat_mult_Vec(n, n, A, x, tmp);
		Vec_add_coefVec(2, -1, tmp, b, rezid);
		Vec_negat(n, rezid);

		/* Daca rezidul este prea mic, returneaza punctul curent */
		if (Vec_norm(n, rezid) < tol) {
			break;
		}
		printVec(n, x);
		fprint(logfile1, x);

		/* Alegere pas optim: arg min f(x+lambda*rezid) => pas_lambda=(rkT*rk)(rkT*A*rk)^-1 */
		lambda_tmp1 = Vec_sprod(n, rezid, rezid);
		lambda_tmp2 = 1 / Vec_Mat_Vec(n, n, rezid, A, rezid);
		lambda = lambda_tmp1 * lambda_tmp2;

		/* Adun pasul la coordonata (x1,x2): x = x + lambda*rezid */
		Vec_add_coefVec(n, lambda, x, rezid, x);
	}
}


/*****************************************************************************/
/** METODA Gradientilor Conjugati: Nr iteratii = Nr dimensiuni f(x1,x2,...) **/
void Gradienti_conjugati(int n,
						 double A[n][n], double b[n],
						 double tol,
						 double x[n]) {
	fprintf(logfile2, "x[0],x[1]\n");

	double lambda, Beta;
	double r[n], d[n], tmp[n];

	// Calcul rezidual
	Mat_mult_Vec(n, n, A, x, tmp);
	Vec_add_coefVec(2, -1, tmp, b, r);
	Vec_negat(n, r);

	// Directia initiala de cautare d0 = r0
	Vec_copy(n, r, d);

	for (int k = 0; k < 8; k++) {
		// Se alege pasul optim lambda = (r_k*d_k)(d_k*A*d_k)^-1
		double tmp1 = Vec_sprod(n, r, d);
		double tmp2 = Vec_Mat_Vec(n, n, d, A, d);
		lambda = tmp1 / tmp2;

		// Pas: x_k+1 = x_k + lambda*r_k
		Vec_add_coefVec(n, lambda, x, r, x);

		// Pas: r_k+1 = -gradient f(x_k)
		double r_tmp[n];
		Vec_copy(n, r, r_tmp);
		
		Mat_mult_Vec(n, n, A, x, tmp);
		Vec_add_coefVec(n, -1, tmp, b, r);
		Vec_negat(n, r);

		// Verifica daca rezidual prea mic
		if (Vec_norm(n, r) < tol) {
			break;
		}
		printVec(n, x);
		fprint(logfile2, x);

		// Alegem o directie conjugata:
		// Beta = norm(r_k+1)^2 * norm(r_k)^2
		tmp1 = pow(Vec_norm(n, r), 2);
		tmp2 = pow(Vec_norm(n, r_tmp),2);
		Beta = tmp1 / tmp2;
		// d_k+1 = r_k+1 + Beta*d_k
		Vec_add_coefVec(n, Beta, r, d, d);
	}
}

double goldsearch(int n, double f(double, double), double dir[], double x[], double tol, double xnew[]);
void a_scale(int n, double coef, double in[], double out[]);
void a_add(int n, double in[], double coef[], double out[]);
/*************************************************************/
/** METODA Gradientilor Conjugati pentru functii nepatratice */
void Gradienti_conjugati_nepatr(int n,
								double f(double, double),
								double fDx(double, double),
								double fDy(double, double),
								double tol,
								double x[n]) {
	fprintf(logfile3, "x[0],x[1]\n");
	double x_new[n];
	double Beta;

	// Calcul gradient al functiei nepatratice
	double grad[n];
	grad[0] = fDx(x[0], x[1]);
	grad[1] = fDy(x[0], x[1]);

	// Calcul directie initiala de cautare:
	double dir[n];
	Vec_copy(n, grad, dir);
	Vec_negat(n, dir);

	for (int k = 0; k < 100; k++) {
		// Pas: Calculam pasul optim lambda folosind o metoda line search
		printf("Iteratie %d\n", k);
		double lambda = goldsearch(n, f, dir, x, tol, x_new);
		printf("\tlambda: "); print(lambda);

		// Pas: x_k+1 = x_k + lambda*r_k
		Vec_add_coefVec(n, lambda, x, dir, x);
		printVec(n, x);
		fprint(logfile3, x);

		// Recalculam grad in functie de noua pozitie:
		double grad_tmp[n];
		Vec_copy(n, grad, grad_tmp);
		grad[0] = fDx(x[0], x[1]);
		grad[1] = fDy(x[0], x[1]);

		// Calculam Beta
		Beta = pow(Vec_norm(n, grad),2) / pow(Vec_norm(n, grad_tmp), 2);

		// Recalculam directia:
		Vec_negat(n, grad);
		Vec_add_coefVec(n, Beta, grad, dir, dir);
	}
}

/** Gold Search **/
double goldsearch(int n,
				  double f(double, double),
				  double dir[],
				  double x[],
				  double tol,
				  double xnew[]){

		double gr = (1 + sqrt(5)) / 2;
		//gr = 2.23;

		double step[n];
		double alphar =  1/gr;
		a_scale(2, alphar, dir, step);
		a_add(2,x,step,xnew);

		double fr = f(xnew[0], xnew[1]);
		double alphal =  1/(gr*gr);
		a_scale(2, alphal, dir, step);
		a_add(2,x,step,xnew);

		double fl = f(xnew[0], xnew[1]);
		double ngrad = Vec_norm(2, dir);
		double dist = (alphar-alphal);

		while (ngrad * dist > tol) {
			dist = (alphar - alphal) / gr;
			if (fr > fl) {
				alphar = alphal;
				fr = fl;
				alphal -= dist;
				a_scale(2, alphal, dir, step);
				a_add(2, x, step, xnew);
				fl = f(xnew[0], xnew[1]);
			}
			else {
				alphal = alphar;
				fl = fr;
				alphar += dist;
				a_scale(2, alphar, dir, step);
				a_add(2, x, step, xnew);
				fr = f(xnew[0], xnew[1]);
			}
		}
		//print((alphal+alphar)/2);
		return (alphal+alphar)/2;
}


/*************
**** MAIN ****
*************/
int main() {
	/** SCRIERE IN FISIERE **/
	logfile1 = fopen("log1.csv","w");
	logfile2 = fopen("log2.csv","w");
	logfile3 = fopen("log3.csv","w");
	if(logfile1==NULL || logfile2==NULL || logfile3==NULL) {
		printf("Could not open one of files, sorry!\n");
		return 1;
	}

	/** FUNCTIA PARABOLOID 1 **/
	// f(x,y) = pow(x,2) + pow(y,2) - x*y - x - y
	double A[2][2] = {{2,-1},{-1,2}};
	double b[2] = {1, 1};

	/** FUNCTIA PARABOLOID 2 **/
	// f(x,y) = 1.5*pow(x,2) + 3*pow(y,2) + 2*x*y - 2*x +8*y
	double A2[2][2] = {{3,2},{2,6}};
	double b2[2] = {2, -8};

	printf("/** METODA CAUCHY **/\n");
	int n = 2;
	double tol = 1e-7;
	double x_start[2] = {0, 1};
	// Cauchy(n, A, b, tol, x_start);

	printf("\n/** METODA Gradientilor Conjugati **/\n");
	x_start[0] = 0;
	x_start[1] = 1;
	// Gradienti_conjugati(n, A, b, tol, x_start);

	printf("\n/** METODA Gradientilor Conjugati pentru functii nepatratice **/\n");
	x_start[0] = 0;
	x_start[1] = 1;
	double x_new[2];
	// Gradienti_conjugati_nepatr(n,
	//                            f2_rosenbrock,
	//                            f2_rosenbrock_Dx,
	//                            f2_rosenbrock_Dy,
	//                            tol, x_start);

	//goldsearch(2, f2_rosenbrock, dir, x_start, tol, x_new);

	x_start[0] = 9;
	x_start[1] = 7;
	Gradienti_conjugati_nepatr(2, f2_paraboloid_xy, f2_paraboloid_xy_Dx, f2_paraboloid_xy_Dy, tol, x_start);


	//// TESTARE FUNCTII AJUTATOARE DE ALGEBRA LINIARA ///////
	/*////////////////////////////////////////////////////////
	double A[2][2] = {{2,3},{1,2}};
	double B[2] = {2,3};
	double C[2];
	Mat_mult_Vec(2,2,A,B,C);
	printVec(2,C);
	*/

	/*////////////////////////////////////////////////////////
	double A[2] = {3,4};
	printf("%.2f", Vec_norm(2, A));

	double a[2] = {2,2}, b[2] = {2,3};
	printf("%.2f", Vec_sprod(2,a,b));
	*/
	/*////////////////////////////////////////////////////////
	double a[2] = {1,2}, b[2] = {2,1}, A[2][2] = {{2,2},{2,2}};
	printf("%.2f", Vec_Mat_Vec(2,2,a,A,b));
	*/
	return 1;
}



void a_scale(int n, double coef, double in[], double out[]) {
	for (int i = 0; i < n; i++)
		out[i] = coef*in[i];
}
void a_add(int n, double in[], double step[], double out[n]) {
	for (int i = 0; i < n; i++)
		out[i] = in[i]+step[i];
}

// rezidual == Gradientul negativ = -(Ax - b)
// rezidual prea mic == Daca norma acestu rezidual este mai mica decat toleranta
// lambda = (rkT*rk)(rkT*A*rk)^-1 == Pasul optim (obtinut din prima analitic)
// Pentru metoda Cauchy, practic ne folosim de ce cunoastem acum in interiorul functiei (gradientul)
