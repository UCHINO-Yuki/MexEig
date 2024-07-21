#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *jobl, char *jobr, double *A, double *B, double *VL, double *VR, double *eigr, double *eigi, mwSignedIndex *n);
void dgeev(char *JOBVL, char *JOBVR, mwSignedIndex *N, double *A, mwSignedIndex *LDA, double *WR, double *WI, 
            double *VL, mwSignedIndex *LDVL, double *VR, mwSignedIndex *LDVR, double *WORK, mwSignedIndex *LWORK, mwSignedIndex *INFO);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Scalars
	mwSignedIndex n, typelen, cflag, i, mflag, Dindex;
    char jobl, jobr;
    char *type;

    // Arrays
	double *A, *B, *VL, *VR, *eigr, *eigi;

    // Arguments
	n = mxGetN(prhs[0]);
	A = mxGetPr(prhs[0]);
    B = MALLOC(double, n*n);
    eigr = MALLOC(double, n);
    eigi = MALLOC(double, n);

    if (nlhs == 1) {
        //-----
        // D = dgeev(A,...)
        //-----
        jobl = 'N';
        jobr = 'N';
        mflag = 0;
        Dindex = 0;

    } else if (nlhs == 2) {
        //-----
        // [VL,D] = dgeev(A,...)
        //-----
        jobl = 'V';
        jobr = 'N';
        plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
        VL = mxGetPr(plhs[0]);
        mflag = 1;
        Dindex = 1;

    } else if (nlhs == 3) {
        //-----
        // [VL,D,VR] = dgeev(A,...)
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
		free(B);
		free(eigr);
		free(eigi);
        mexErrMsgIdAndTxt("dgeev:TooManyOutput", "Too many output.");
    }

    // Run
    run_vec(&jobl, &jobr, A, B, VL, VR, eigr, eigi, &n);
    
    // check complex
    cflag = 0;
    for (i=n-1; i>=0; --i) {
        if (eigi[i] != 0) {
            cflag = 1;
            break;
        }
    }

    if (nrhs > 1) {
        if (mxIsChar(prhs[1])) {
            typelen = mxGetN(prhs[1]) + 1;
            type = mxMalloc(typelen);
            if (mxGetString(prhs[1], type, (mwSize) typelen) == 0) {
                if (strcmp(type, "matrix") == 0) mflag = 1;
                else mflag = 0;
            }
        }
    }

    // Return values
    if (mflag == 0) {
        if (cflag == 0) { // all eigenvalues are not complex number
            double *eigenvalue;
            plhs[Dindex] = mxCreateDoubleMatrix(n, 1, mxREAL);
            eigenvalue = mxGetPr(plhs[Dindex]);
            for (i=n-1; i>=0; --i) eigenvalue[i] = eigr[i];
        } else {
            mxComplexDouble *eigenvalue;
            plhs[Dindex] = mxCreateDoubleMatrix(n, 1, mxCOMPLEX);
            eigenvalue = mxGetComplexDoubles(plhs[Dindex]);
            for (i=n-1; i>=0; --i) {
                eigenvalue[i].real = eigr[i];
                eigenvalue[i].imag = eigi[i];
            }
        }
    } else {
        if (cflag == 0) { // all eigenvalues are not complex number
            double *eigenvalue;
            plhs[Dindex] = mxCreateDoubleMatrix(n, n, mxREAL);
			eigenvalue = mxGetPr(plhs[Dindex]);
			for (i=n*n-1; i>=0; --i) eigenvalue[i] = 0;
			for (i=n-1; i>=0; --i) eigenvalue[i*n + i] = eigr[i];
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
				eigenvalue[i*n + i].real = eigr[i];
				eigenvalue[i*n + i].imag = eigi[i];
			}
		}
    } 

	free(B);
	free(eigr);
	free(eigi);
}

void run_vec(char *jobl, char *jobr, double *A, double *B, double *VL, double *VR, double *eigr, double *eigi, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, info;
    
    // Arrays
	double *work;
    double dammy[1];

    for (i=0; i<(*n)*(*n); i++) B[i] = A[i];

    // Get the length
	dgeev(jobl, jobr, n, B, n, eigr, eigi, VL, n, VR, n, dammy, &lwork, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(double, lwork);

    // Run
	dgeev(jobl, jobr, n, B, n, eigr, eigi, VL, n, VR, n, work, &lwork, &info);
    
    // Release the workspace
    free(work);
}
