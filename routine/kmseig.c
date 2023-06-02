#include "mex.h"
#include "math.h"

/* Trench, Spectral decomposition of Kac-Murdock-Szego matrices */
void kmseig(mwSize N, double rho, double *E, double *V)
{
	double theta, theta0, theta1, g, g0, g1, eta, value, norm, degree;
	int sign = 1;
        int i;
	  for ( i = 0; i < N; i++) {
		theta0 = i * M_PI / (N + 1);
		theta1 = (i + 1) * M_PI / (N + 1);
		g0 = sin((N + 1)*theta0) - 2 * rho * sin(N * theta0) + rho*rho * sin((N - 1)*theta0);
		g1 = sin((N + 1)*theta1) - 2 * rho * sin(N * theta1) + rho*rho * sin((N - 1)*theta1);
		while (theta1 - theta0 > 0.000000000001) {
			theta = (theta0 + theta1) / 2;
			g = sin((N + 1)*theta) - 2 * rho * sin(N * theta) + rho*rho * sin((N - 1)*theta);
			if (g*g1 <= 0) {
				theta0 = theta;
				g0 = g;
			}
			else {
				theta1 = theta;
				g1 = g;
			}
		}
		theta = (theta0 + theta1) / 2;
		E[i] = (1 - rho*rho) / (1 - 2 * rho * cos(theta) + rho*rho);
		eta = 2 * atan(rho*sin(theta) / (1 - rho*cos(theta)));
		norm = 0;
		degree = (i + 1)*M_PI / (N + 1) - eta / (N + 1);
                int j;
		for (j = 0; j < (N + 1) / 2; j++) {
			value = sin((j + 1) * degree + eta / 2);
			V[j + i*N] = value;
			norm += 2 * value * value;
		}
		if (N % 2 == 1) {
			norm -= value * value;
			norm = sqrt(norm);
			V[(N - 1) / 2 + i*N] /= norm;
		}
		else
			norm = sqrt(norm);

		for (j = 0; j < N / 2; j++) {
			V[j + i*N] /= norm;
			V[(i + 1)*N - j - 1] = sign * V[j + i*N];
		}
		sign = -sign;
	}
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mwSize N;    // size of matrix 
	double rho;    // input value  
	double *outEigValue;    // vector of eigenvalues
	double *outEigVec;    // matrix of eigenvectors

	if (nrhs != 2)
		mexErrMsgTxt("Must have exactly two input arguments.");
	if (nlhs > 2)
		mexErrMsgTxt("Too many output arguments.");

	N = *(mxGetPr(prhs[0]));
	rho = *(mxGetPr(prhs[1]));
	plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);    // create the output eigenvalue vector
	plhs[1] = mxCreateDoubleMatrix(N, N, mxREAL);    // create the output eigenvector matrix
	outEigValue = mxGetPr(plhs[0]);    // Nx1 output matrix
	outEigVec = mxGetPr(plhs[1]);    // NxN output matrix
	kmseig(N, rho, outEigValue, outEigVec);    // call the computational routine
}
