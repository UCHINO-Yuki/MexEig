#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *jobvl, char *jobvr, double *A, double *B, double *C, double *D, double *VL, double *VR, double *ar, double *ai, double *b, mwSignedIndex *n);
void dggev(char *JOBVL, char *JOBVR,
           mwSignedIndex *N, double *A, mwSignedIndex *LDA, double *B, mwSignedIndex *LDB,
           double *ALPHAR, double *ALPHAI, double *BETA, 
           double *VL, mwSignedIndex *LDVL, double *VR, mwSignedIndex *LDVR,
           double *WORK, mwSignedIndex *LWORK, mwSignedIndex *INFO);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 2)
        mexErrMsgIdAndTxt("dggev:inputNotEnough", "Input argument is not enough.");

    // Scalars
	mwSignedIndex n, typelen, cflag, i, mflag, Dindex;
    char jobl, jobr;
    char *type;
    
    // Arrays
	double *A, *B, *C, *D, *VL, *VR, *eigenvalue, *ar, *ai, *b;

    // Arguments
	n = (mwSignedIndex)mxGetN(prhs[0]);
	A = mxGetPr(prhs[0]);
	B = mxGetPr(prhs[1]);
    C = MALLOC(double, n*n);
    D = MALLOC(double, n*n);
    ar = MALLOC(double, n);
    ai = MALLOC(double, n);
    b = MALLOC(double, n);

    if (nlhs == 1) {
        //-----
        // D = dggev(A,B,...)
        //-----
        jobl = 'N';
        jobr = 'N';
        mflag = 0;
        Dindex = 0;

    } else if (nlhs == 2) {
        //-----
        // [VL,D] = dggev(A,B,...)
        //-----
        jobl = 'V';
        jobr = 'N';
        plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
        VL = mxGetPr(plhs[0]);
        mflag = 1;
        Dindex = 1;

    } else if (nlhs == 3) {
        //-----
        // [VL,D,VR] = dggev(A,B,...)
        //-----
        jobl = 'V';
        jobr = 'V';
        plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
        VL = mxGetPr(plhs[0]);
        plhs[2] = mxCreateDoubleMatrix(n, n, mxREAL);
        VR = mxGetPr(plhs[2]);
        mflag = 1;
        Dindex = 1;
        
    } else {
        free(C);
		free(D);
		free(ar);
		free(ai);
		free(b);
        mexErrMsgIdAndTxt("dggev:TooManyOutput", "Too many output.");
    }

    // Run
    run_vec(&jobl, &jobr, A, B, C, D, VL, VR, ar, ai, b, &n);

    // check complex
    cflag = 0;
    for (i=n-1; i>=0; --i) {
        if (ai[i] != 0) {
            cflag = 1;
            break;
        }
    }
    
    if (nrhs == 3) {
        //-----
        // [...] = dggev(A,B,'...')
        //-----
        if (mxIsChar(prhs[2])) {
            typelen = mxGetN(prhs[2]) + 1;
            type = mxMalloc(typelen);
            if (mxGetString(prhs[2], type, (mwSize) typelen) == 0) {
                if (strcmp(type, "matrix") == 0) mflag = 1;
                else if (strcmp(type, "vector") == 0) mflag = 0;
            }
        }
    }

    // Return values
    if (mflag == 0) {
        if (cflag == 0) { // all eigenvalues are not complex number
            double *eigenvalue;
            plhs[Dindex] = mxCreateDoubleMatrix(n, 1, mxREAL);
            eigenvalue = mxGetPr(plhs[Dindex]);
            for (i=n-1; i>=0; --i) eigenvalue[i] = ar[i]/b[i];
        } else {
            mxComplexDouble *eigenvalue;
            plhs[Dindex] = mxCreateDoubleMatrix(n, 1, mxCOMPLEX);
            eigenvalue = mxGetComplexDoubles(plhs[Dindex]);
            for (i=n-1; i>=0; --i) {
                eigenvalue[i].real = ar[i]/b[i];
                eigenvalue[i].imag = ai[i]/b[i];
            }
        }
    } else {
        if (cflag == 0) { // all eigenvalues are not complex number
            double *eigenvalue;
            plhs[Dindex] = mxCreateDoubleMatrix(n, n, mxREAL);
			eigenvalue = mxGetPr(plhs[Dindex]);
			for (i=n*n-1; i>=0; --i) eigenvalue[i] = 0;
			for (i=n-1; i>=0; --i) eigenvalue[i*n + i] = ar[i]/b[i];
		}
		else {
            mxComplexDouble *eigenvalue;
			plhs[Dindex] = mxCreateDoubleMatrix(n, n, mxCOMPLEX);
			eigenvalue = mxGetComplexDoubles(plhs[Dindex]);
			for (i=n*n-1; i>=0; --i) {
				eigenvalue[i].real = 0;
				eigenvalue[i].imag = 0;
			}
			for (i=n-1; i>=0; --i) {
                eigenvalue[i*n + i].real = ar[i]/b[i];
                eigenvalue[i*n + i].imag = ai[i]/b[i];
			}
		}
    } 

    free(C);
    free(D);
    free(ar);
    free(ai);
    free(b);
}

void run_vec(char *jobvl, char *jobvr, double *A, double *B, double *C, double *D, double *VL, double *VR, double *ar, double *ai, double *b, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, info=0;
    
    // Arrays
	double *work;
    double dammy[1];

    for (i=(*n)*(*n)-1; i>=0; --i) {
        C[i] = A[i];
        D[i] = B[i];
    }

    // Get the length
    dggev(jobvl, jobvr, n, C, n, D, n, ar, ai, b, VL, n, VR, n, dammy, &lwork, &info);
    
    // Set the workspace
    lwork=(mwSignedIndex)dammy[0];
    work = MALLOC(double, lwork);

    // Run
    dggev(jobvl, jobvr, n, C, n, D, n, ar, ai, b, VL, n, VR, n, work, &lwork, &info);
    
    // Release the workspace
    free(work);
}
