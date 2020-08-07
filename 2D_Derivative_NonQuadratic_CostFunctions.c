/* TOP LAB 5 - Metode de tip Newton */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

FILE *logfile1, *logfile2;

/* Functii de test */
double f1_paraboloid(double x, double y) {
	return x*x + y*y - x*y - x - y; }
double f1_paraboloid_Dx(double x, double y) {
	return 2*x - y - 1; }
double f1_paraboloid_Dy(double x, double y) {
	return -x + 2*y - 1; }

double f2_rosenbrock(double x, double y) {
	return pow((1 - x), 2) + 100 * pow((y - x * x), 2); }
double f2_rosenbrock_Dx(double x, double y) {
	return 2 * (200 * x * x * x - 200 * x * y + x - 1); }
double f2_rosenbrock_Dy(double x, double y) {
	return 200 * (y - x * x); }

double f3_himenblau(double x, double y) {
	return pow((x * x + y - 11), 2) + pow((x + y * y - 7), 2); }
double f3_himenblau_Dx(double x, double y) {
	return 2 * (2 * x * (x * x + y - 11) + x + y * y - 7); }
double f3_himenblau_Dy(double x, double y) {
	return 2 * (x * x + 2 * y * (x + y * y - 7) + y - 11); }

/******************************************
/* Functii ajutatoare de afisare (debug) */
void fprint(FILE *log, double vec[2]);
void print(double x);
void printVec(int n, double vec[]);
void printMat(int n, int m, double mat[n][m]);

/*******************************
/* Functii de algebra liniara */
void Vec_add(int n, double in[], double step[], double out[n]);
void Vec_add_coefVec(int n, double a[n], double coef, double b[n], double out[n]);
void Vec_scale(int n, double coef, double in[n], double out[n]);
void Vec_mult_Vec(int n, double a[n], double b[n], double out[n][n]);
void Mat_mult_Vec(int n, int m, double mat[n][m], double vec[m], double out[n]);
void Mat_scale(int n, int m, double in[n][n], double coef, double out[n][n]);
void Mat_add_Mat(int n, int m, double A[n][n], double B[n][n], double OUT[n][n]);
double Vec_Mat_Vec(int n, int m, double va[n], double mat[n][m], double vb[m]);
double Vec_norm(int n, double vec[n]);
double Vec_sprod(int n, double a[n], double b[n]);

void Mat_inv(int n, double mat[n][n], double mat_inv[n][n]);

double goldsearch(int n, double f(double, double), double dir[], double x[], double tol, double xnew[]);

/*******************************
/* Implementarea ALGORITMILOR */
/*******************/
void Gradienti_Conjugati_Rez(int n, double A[n][n], double b[n], double out[n]);

/** METODA Newton Clasica **/
void Newton_Clasic(int n,
	double f(double, double),
	double fDx(double, double),
	double fDy(double, double),
	double tol,
	double x[n]) {

	fprintf(logfile1, "x[0],x[1],f(x0 x1)\n");
	double grad[n], pas[n], H[n][n];

	int nr_iter = 50;
	for (int k = 0; k < nr_iter; k++) {

		// Calculeaza gradientul
		grad[0] = fDx(x[0], x[1]);
		grad[1] = fDy(x[0], x[1]);

		/* ESTIMEAZA HESSIANA */
		double x_tmp[n], grad_tmp[n];
		Vec_scale(n, 1, x, x_tmp);
		// Parcurgere pentru cele 2 linii ale hessianei
		for (int i = 0; i < n; i++) {
			pas[i] = fabs(x[i]) * sqrt(tol);
			x_tmp[i] += pas[i];
			grad_tmp[0] = fDx(x_tmp[0], x_tmp[1]);
			grad_tmp[1] = fDy(x_tmp[0], x_tmp[1]);
			x_tmp[i] -= pas[i];

			Vec_add_coefVec(n, grad_tmp, -1, grad, H[i]);
			if (pas[i] != 0) { // evitam divizarea cu zero
				Vec_scale(n, 1/pas[i], H[i], H[i]);
			}
		}
		// printMat(n, n, H);

		/* Conditie stop daca gradientul este prea mic */
		if (Vec_norm(n, grad) < tol) {
			printf("Conditia stop norm(r) < tol dupa %d iteratii.\n", k);
			break;
		}

		/* Calculeaza H^-1 * grad ... */
		/*** MET 1: CALCULAM INVERSA ***/
		double H_inv[n][n], p[n];
		Mat_inv(n, H, H_inv);
		Mat_mult_Vec(n, n, H_inv, grad, p);

		/*** MET 2: Rezolvam cu Metoda Grad Conjugati ***/
		// /* OBS !!! In cazul Rosenbrock, Hessiana nu va fi pozitiv definita
		// => Metoda nu va converge catre solutia A_inv*b !!! */
		// double p[n];
		// Gradienti_Conjugati_Rez(n, H, grad, p);
		// printVec(n, p);

		// Mutam x
		Vec_add_coefVec(n, x, -1, p, x);

		fprintf(logfile1, "%.5f, %.5f, %.5f\n", x[0], x[1], f(x[0], x[1]));
		printf("%.5f, %.5f,		f(x0,x1)=%.5f\n", x[0], x[1], f(x[0], x[1]));
	}
}

/** METODA Quasi-Newton **/
void Quasi_Newton(int n,
	double f(double, double),
	double fDx(double, double),
	double fDy(double, double),
	double tol,
	double x[n]) {

	fprintf(logfile2, "x[0],x[1],f(x0 x1)\n");
	/* Este extensia in mai multe dimensiuni a metodei Secantei */
	// Incercam sa estimam inversa Hessienei H in mod direct, unde B va fi B = H
	//fprintf(logfile2, "x[0],x[1],f(x0 x1)\n");

	double s[n], q[n], p[n], lambda;

	// Matricea identitate
	double B[n][n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i == j)
				B[i][j] = 1;
			else
				B[i][j] = 0;

	// Se estimeaza gradient g
	double grad[n];
	grad[0] = fDx(x[0], x[1]);
	grad[1] = fDy(x[0], x[1]);

	int nr_iter = 10000;
	for (int k = 0; k < nr_iter; k++) {
		printf("\tk = %d\n", k);
		// Conditie Stop: norma gradientului prea mica
		if (Vec_norm(n, grad) < tol) {
			printf("Conditia stop norm(r) < tol dupa %d iteratii.\n", k);
			break;
		}

		// Pas 5: Se alege directia de inaintare s = -B*g (estimare echiv: -(H^-1)*g)
		Mat_mult_Vec(n, n, B, grad, s);
		//Vec_scale(n, -1/Vec_norm(n,s), s, s);
		Vec_scale(n, -1, s, s);

		// Pas 6: Se alege pasul de inaintare
		double xnew[n];
		double dir[n];
		Vec_scale(n, -1, grad, dir);
		lambda = goldsearch(n, f, dir, x, tol, xnew);
		printf("lambda: "); print(lambda);
		printf("s: "); printVec(n, s);

		// Pas 7: Salvam x curent si Mutam x
		double x_tmp[n];
		x_tmp[0] = x[0];
		x_tmp[1] = x[1];
		Vec_add_coefVec(n, x, lambda, s, x);

		// Pas 8: Salvam gradient curent si Se estimeaza gradientul g
		double grad_tmp[n];
		grad_tmp[0] = grad[0];
		grad_tmp[1] = grad[1];

		grad[0] = fDx(x[0], x[1]);
		grad[1] = fDy(x[0], x[1]);

		// Pas 9: Aflam q si p
		Vec_add_coefVec(n, grad, -1, grad_tmp, q);
		Vec_add_coefVec(n, x, -1, x_tmp, p);

		// Pas 10: Estimarea Davidon-Fletcher-Powell
		// al doilea termen
		double pp_T[n][n];
		Vec_mult_Vec(n, p, p, pp_T);
		Mat_scale(n, n, pp_T, 1/Vec_sprod(n, p, q), pp_T);

		// al treilea termen
		double Bq[n], q_TB[n], Bqq_TB[n][n], q_TBq;
		Mat_mult_Vec(n, n, B, q, Bq);

		Mat_mult_Vec(n, n, B, q, q_TB);
		Vec_mult_Vec(n, Bq, q_TB, Bqq_TB);
		q_TBq = Vec_Mat_Vec(n, n, q, B, q);
		Mat_scale(n, n, Bqq_TB, 1/q_TBq, Bqq_TB);

		// B = B + pp_T/p_Tq - Bqq_TB/q_TBq
		double tmp_B[n][n];
		Mat_add_Mat(n, n, B, pp_T, tmp_B);
		Mat_scale(n, n, Bqq_TB, -1, Bqq_TB);
		Mat_add_Mat(n, n, tmp_B, Bqq_TB, B);

		//printf("B: "); printMat(n,n,B);
		//printf("x: "); printVec(n, x);
		fprintf(logfile2, "%.5f, %.5f, %.5f\n", x[0], x[1], f(x[0], x[1]));
		printf("%.5f, %.5f,		f(x0,x1)=%.5f\n", x[0], x[1], f(x[0], x[1]));
	}
}

/********/
/* MAIN */
/********/
int main() {
	/** SCRIERE IN FISIERE **/
	logfile1 = fopen("log1_himenblau.csv","w");
	logfile2 = fopen("log2_himenblau.csv","w");
	if(logfile1==NULL || logfile2==NULL) {
		printf("Could not open one of files, sorry!\n");
		return 1;
	}

	printf("/** METODA Newton Clasica **/\n");
	int n = 2;
	double tol = 1e-3;
	double x_start[2] = {9, 7};
	Newton_Clasic(n, f1_paraboloid, f1_paraboloid_Dx, f1_paraboloid_Dy, tol, x_start);
	// Newton_Clasic(n, f2_rosenbrock, f2_rosenbrock_Dx, f2_rosenbrock_Dy, tol, x_start);
	// Newton_Clasic(n, f3_himenblau, f3_himenblau_Dx, f3_himenblau_Dy, tol, x_start);

	printf("\n** METODA Quasi-Newton **/\n");
	x_start[0] = -0.5;
	x_start[1] = -0.5;
	// Quasi_Newton(n, f1_paraboloid, f1_paraboloid_Dx, f1_paraboloid_Dy, tol, x_start);
	// Quasi_Newton(n, f2_rosenbrock, f2_rosenbrock_Dx, f2_rosenbrock_Dy, tol, x_start);
	// Quasi_Newton(n, f3_himenblau, f3_himenblau_Dx, f3_himenblau_Dy, tol, x_start);

	return 1;
}

/*****************************************************************************/
/** METODA Gradientilor Conjugati: Pentru rezolvarea ecuatiilor A^-1 * g(x) **/
void Gradienti_Conjugati_Rez(int n, double A[n][n], double b[n], double out[n]) {
	/* OBS: Metoda nu converge daca A nu este pozitiv definita*/
	double lambda, Beta, tol = 1e-4;
	double r[n], d[n], tmp[n];
	double x[n];
	x[0] = 0.1;
	x[1] = 0.1;

	// Calcul rezidual
	Mat_mult_Vec(n, n, A, x, tmp);
	Vec_add_coefVec(2, tmp, -1, b, r);
	Vec_scale(n, -1, r, r);

	// Directia initiala de cautare d0 = r0
	Vec_scale(n, -1, r, d);

	while (Vec_norm(n, r) > tol) {
		// Se alege pasul optim lambda = (r_k*d_k)(d_k*A*d_k)^-1
		double tmp1 = Vec_sprod(n, r, d);
		double tmp2 = Vec_Mat_Vec(n, n, d, A, d);
		lambda = tmp1 / tmp2;

		// Pas 5: x_k+1 = x_k + lambda*r_k
		Vec_add_coefVec(n, x, lambda, r, x);

		// Pas 6: r_k+1 = -gradient f(x_k), recalculam gradient
		double r_tmp[n];
		Vec_scale(n, 1, r, r_tmp);

		Mat_mult_Vec(n, n, A, x, tmp);
		Vec_add_coefVec(n, tmp, -1, b, r);
		Vec_scale(n, -1, r, r);

		// Pas 8: Alegem o directie conjugata:
		// Beta = norm(r_k+1)^2 * norm(r_k)^2
		tmp1 = pow(Vec_norm(n, r), 2);
		tmp2 = pow(Vec_norm(n, r_tmp),2);
		Beta = tmp1 / tmp2;
		// d_k+1 = r_k+1 + Beta*d_k
		Vec_add_coefVec(n, r, Beta, d, d);

		//printVec(n, x);
	}

	Vec_scale(n, 1, x, out);
}


/******************************************
/* Functii ajutatoare de afisare (debug) */
void fprint(FILE *log, double vec[2]) {
	fprintf(log, "%.4f,%.4f\n", vec[0], vec[1]);
}
void print(double x) {
	printf("%.4f\n", x);
}
void printVec(int n, double vec[]) {
	printf("vec = [");
	for (int i = 0; i < n; i++) {
		printf(" %.4f", vec[i]);
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
void Vec_add_coefVec(int n, double a[n], double coef, double b[n], double out[n]) {
	/**a_scaleadd**/
	for (int i = 0; i < n; i++) {
		out[i] = a[i] + coef * b[i];
	}
}

void Vec_scale(int n, double coef, double in[n], double out[n]) {
	for (int i = 0; i < n; i++) {
		out[i] = coef * in[i];
	}
}

void Mat_scale(int n, int m, double in[n][n], double coef, double out[n][n]) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			out[i][j] = coef * in[i][j];
		}
	}
}

void Mat_mult_Vec(int n, int m, double mat[n][m], double vec[m], double out[n]) {
	/**mat_vec: Inmultire matrice cu vector (pt aflare gradient) **/
	// Initializam out cu zero
	for (int i = 0; i < n; i++) {
		out[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			out[i] += mat[i][j] * vec[j];
		}
	}
}

void Vec_mult_Vec(int n, double a[n], double b[n], double out[n][n]){
	// inmultire vector coloana p cu vector linie pT pentru a obtine o matrice (Vector Outer Product):
	// p(n=2,m=1) * pT(n=1,m=2) = M(n=2,m=2)
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			out[i][j] = a[i] * b[j];
		}
	}
}

double Vec_Mat_Vec(int n, int m, double va[n], double mat[n][m], double vb[m]) {
	/** Sum(a[i]*A[i][j]*b[j]) **/
	// Prima parte: va(2,1)*A(2,2)=tmp(1,2)
	double va_mat[n], res = 0;
	for (int i = 0; i < n; i++) {
		va_mat[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			va_mat[i] += va[j] * mat[j][i];
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

void Mat_add_Mat(int n, int m, double A[n][n], double B[n][n], double OUT[n][n]) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			OUT[i][j] = A[i][j] + B[i][j];
		}
	}
}

void Mat_inv(int n, double mat[n][n], double mat_inv[n][n]) {
	//finding determinant
	double determinant;
	for(int i = 0; i < n; i++)
		determinant = determinant + (mat[0][i] * (mat[1][(i+1)%3] * mat[2][(i+2)%3] - mat[1][(i+2)%3] * mat[2][(i+1)%3]));

	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++)
			mat_inv[i][j] = ((mat[(j+1)%3][(i+1)%3] * mat[(j+2)%3][(i+2)%3]) - (mat[(j+1)%3][(i+2)%3] * mat[(j+2)%3][(i+1)%3]))/ determinant;
	}
}

void Vec_add(int n, double in[], double step[], double out[n]) {
	for (int i = 0; i < n; i++)
		out[i] = in[i]+step[i];
}


		// // Quasi-Newton
		// // Pas 6: Se alege pasul de inaintare folosind un algoritm 1D
		// // (ex: Secantei, avand derivata 1 si folosind-o pentru a aproxima deriv a 2a)
		// double x_tmp0 = x[0];
		// double x_tmp1 = x[1];
		// for (int i = 0; i < 10; i++) {
		// 	double x_tmp_next0 = x_tmp0 - fDx(x_tmp0, x_tmp1);
		// 	double x_tmp_next1 = x_tmp1 - fDy(x_tmp0, x_tmp1);

		// 	double fdd_x0 = (fDx(x_tmp_next0, x_tmp_next1) - fDx(x_tmp0, x_tmp1)) / (x_tmp_next0-x_tmp0);
		// 	double fdd_x1 = (fDy(x_tmp_next0, x_tmp_next1) - fDy(x_tmp0, x_tmp1)) / (x_tmp_next1-x_tmp1);

		// 	double lambda0 = fDx(x_tmp_next0, x_tmp_next1) / fdd_x0;
		// 	double lambda1 = fDx(x_tmp_next0, x_tmp_next1) / fdd_x1;

		// 	lambda[0] = lambda0;
		// 	lambda[1] = lambda1;

		// 	printVec(n, lambda);

		// 	x_tmp0 = x_tmp_next0;
		// 	x_tmp1 = x_tmp_next1;
		// }

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
		Vec_scale(2, alphar, dir, step);
		Vec_add(2,x,step,xnew);

		double fr = f(xnew[0], xnew[1]);
		double alphal =  1/(gr*gr);
		Vec_scale(2, alphal, dir, step);
		Vec_add(2,x,step,xnew);

		double fl = f(xnew[0], xnew[1]);
		double ngrad = Vec_norm(2, dir);
		double dist = (alphar-alphal);

		while(ngrad*dist>tol){
			dist = (alphar-alphal)/gr;
			if(fr > fl){
				alphar = alphal;
				fr = fl;
				alphal -= dist;
				Vec_scale(2, alphal, dir, step);
				Vec_add(2,x,step,xnew);
				fl = f(xnew[0], xnew[1]);
			}
			else{
				alphal = alphar;
				fl = fr;
				alphar += dist;
				Vec_scale(2, alphar, dir, step);
				Vec_add(2,x,step,xnew);
				fr = f(xnew[0], xnew[1]);
			}
		}
		//print((alphal+alphar)/2);
		return (alphal+alphar)/2;
}
