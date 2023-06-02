#include "mex.h"

/* Michael Stewart,
Cholesky factorization of semi-definite Toeplitz matrices, Linear Algebra and its Applications,
Volume 254, pages 497-525, 1997. */
void toep_chol(mwSize n, double *t, double *L)
{
	double *g = (double *)malloc(2 * n * sizeof (double));
	double div, g1j, g2j, rho, t0;
	int i, j;
    
    t0 = sqrt(t[0]);
	for (j = 0; j < n; j++)
	{
		L[j] = t[j] / t0;
        g[j * 2] = t[j];
        g[1 + j * 2] = t[j];
	}
	g[1] = 0.0;
	for (j = n - 1; 1 <= j; j--)
	{
		g[j * 2] = g[(j - 1) * 2];
	}
    g[0] = 0.0;
	for (i = 1; i < n; i++)
	{
		rho = -g[1 + i * 2] / g[i * 2];
		div = sqrt((1.0 - rho) * (1.0 + rho));
		for (j = i; j < n; j++)
		{
			g1j = g[0 + j * 2];
			g2j = g[1 + j * 2];
			g[0 + j * 2] = (g1j + rho * g2j) / div;
			g[1 + j * 2] = (rho * g1j + g2j) / div;
            L[j + i*n] = g[j * 2] / t0;
		}
		for (j = n - 1; i < j; j--)
		{
			g[j * 2] = g[(j - 1) * 2];
		}
		g[i * 2] = 0.0;
	}
	free(g);
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],	int nrhs, const mxArray *prhs[])
{
	mwSize nr, nc;						/* size of matrix */
	double *inMatrix;               /* 1xN input matrix */
	double *outMatrix;              /* output matrix */

	/* get dimensions of the input matrix */
    nr = mxGetM(prhs[0]);
	nc = mxGetN(prhs[0]);
	if (nr == 0 || nc == 0 )
		mexErrMsgIdAndTxt("MyToolbox:toep_cholesky", "Input matrix cannot be empty.");
    if (nr < nc)
        nr = nc;
    else
        nc = nr;
    
	/* create the output matrix */
	plhs[0] = mxCreateDoubleMatrix(nc, nc, mxREAL);

	/* create a pointer to the real data in the input matrix  */
	inMatrix = mxGetPr(prhs[0]);

	/* get a pointer to the real data in the output matrix */
	outMatrix = mxGetPr(plhs[0]);

	/* call the computational routine */
	toep_chol(nc, inMatrix, outMatrix);
}