#include <math.h>
#include <stdio.h>
#include "rkf45.h"
#include "quanc8.h"
#include "zeroin.h"
#include "float.h"

using namespace std;

double t1 = 0, t2 = 20, h = 0.5, A = 1.2, B = 0, C = -0.40218305, D, M = 1 / 82.45, M1 = 1 - M;
int N = 4; //кол-во уравнений в системе

double integrandC(double x) {
	return (1 / sqrt(((pow(x, 2) + 1) * (3 * pow(x, 2) + 4))));
}

double calculationD(double x) {
	return 2 * log10(x) + 1 - (x / 2);
}

void systemForRKFy(double t, double *y, double *dy) {
	dy[0] = y[1]; //y[1] = y'
	dy[1] = -2 * y[3] + y[0] - M1 * y[0] / pow(pow(pow(y[2] + M, 2) + pow(y[0], 2), 0.5), 3) - M * y[0] / pow(pow(pow(y[2] - M1, 2) + pow(y[0], 2), 0.5), 3); //y1'
	dy[2] = y[3]; //y[3] = x'
	dy[3] = 2 * y[1] + y[2] - M1 * (y[2] + M) / pow(pow(pow(y[2] + M, 2) + pow(y[0], 2), 0.5), 3) - M * (y[2] - M1) / pow(pow(pow(y[2] - M1, 2) + pow(y[0], 2), 0.5), 3); //y2'
}

int main()
{
	//вычисляется значение C
	double a1 = 0, b1 = 1;
	const double ABSERR = 1.0e-14, RELERR = 1.0e-14;
	double errest, flag;
	int nofun, k = 0;
	double q[1];
	quanc8(integrandC, a1, b1, ABSERR, RELERR, &q[0], &errest, &nofun, &flag); //the value of the integral in C
	
	C += q[k]; //прибавка результата интеграла к остальной части выражения
	C = pow(C, 4); //получение C
	printf("\nC = %.64f", C);
	printf(";      errest C = %.32f\n", errest);

	//вычисляется значение D
	double a2 = 1, tol = 1.0e-14;
	double b2 = DBL_MAX;
	double resultD = Zeroin(calculationD, a2, b2, tol);
	D = -2.639594 * ((double)resultD);
	printf("\nD = %.14f", D);
	printf(";       tol D = %.32f\n", tol);

	//calculating system
	//t - 
	double *t = new double(t1), *tout = new double(t1), *re = new double(1e-13), *ae = new double(1e-13),
           *t5 = new double(t1), *tout5 = new double(t1), *re5 = new double(1e-13), *ae5 = new double(1e-13),
           *t10 = new double(t1), *tout10 = new double(t1), *re10 = new double(1e-13), *ae10 = new double(1e-13),
           *t25 = new double(t1), *tout25 = new double(t1), *re25 = new double(1e-13), *ae25 = new double(1e-13);

	int flagY = 1, flagY5 = 1, flagY10 = 1, flagY25 = 1;
	int iwork[10], iwork5[10], iwork10[10], iwork25[10];
	double work[30], work5[30], work10[30], work25[30];

	double y[4] = {1.2, C, D, 0};
	double y5[4] = {1.2/100*95, C/100*95, D/100*95, 0};
	double y10[4] = {1.2/100*90, C/100*90, D/100*90, 0};
	double y25[4] = {1.2/100*75, C/100*75, D/100*75, 0};

	int i = 0;
	*tout = 0;
	printf("\nN      T           Y                     X\n");
	
	for (*tout = t1; *tout < t2 + h; *tout += h) {

		RKF45(systemForRKFy, N, y, t, tout, re, ae, &flagY, work, iwork);

		printf("%3d", i);
		printf("%10.2f", *tout);
		printf("%22.14f", y[0]);
		printf("%22.14f\n", y[2]);
		i++;
	}
	
	i = 0;
	*tout = 0;
	printf("\nN      T        Y error 5%%            X error 5%%\n");
	
	for (*tout = t1; *tout < t2 + h; *tout += h) {

		RKF45(systemForRKFy, N, y5, t5, tout, re5, ae5, &flagY5, work5, iwork5);

		printf("%3d", i);
		printf("%10.2f", *tout);
		printf("%22.14f", y5[0]);
		printf("%22.14f\n", y5[2]);
		i++;
	}
	
	i = 0;
	*tout = 0;
	printf("\nN      T        Y error 10%%            X error 10%%\n");	
	for (*tout = t1; *tout < t2 + h; *tout += h) {

		RKF45(systemForRKFy, N, y10, t10, tout, re10, ae10, &flagY10, work10, iwork10);

		printf("%3d", i);
		printf("%10.2f", *tout);
		printf("%22.14f", y10[0]);
		printf("%22.14f\n", y10[2]);
		i++;
	}
	
	i = 0;
	*tout = 0;
	printf("\nN      T        Y error 25%%            X error 25%%\n");	
	for (*tout = t1; *tout < t2 + h; *tout += h) {

		RKF45(systemForRKFy, N, y25, t25, tout, re25, ae25, &flagY25, work25, iwork25);

		printf("%3d", i);
		printf("%10.2f", *tout);
		printf("%22.14f", y25[0]);
		printf("%22.14f\n", y25[2]);
		i++;
	}

	return 0;
}


