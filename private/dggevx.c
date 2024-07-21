#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *BALANC, char *jobl, char *jobr, double *A, double *B, double *C, double *D, double *VL, double *VR, double *ar, double *ai, double *b, mwSignedIndex *n);
void dggevx(char *BALANC, char *JOBVL, char *JOBVR, char *SENSE,
           mwSignedIndex *N, double *A, mwSignedIndex *LDA, double *B, mwSignedIndex *LDB,
           double *ALPHAR, double *ALPHAI, double *BETA, 
           double *VL, mwSignedIndex *LDVL, double *VR, mwSignedIndex *LDVR,
           mwSignedIndex *ILO, mwSignedIndex *IHI, double *LSCALE, double *RSCALE, 
           double *ABNRM, double *BBNRM, double *RCONDE, double *RCONDV, 
           double *WORK, mwSignedIndex *LWORK, mwSignedIndex *IWORK, mxLogical *BWORK, mwSignedIndex *INFO);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Scalars
	mwSignedIndex n, typelen, cflag, i, mflag, Dindex;
    char jobl, jobr, NN = 'N';
    char *type, *BALANC;
    BALANC = &NN;

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
        // D = dggevx(A,...)
        //-----
        jobl = 'N';
        jobr = 'N';
        mflag = 0;
        Dindex = 0;

    } else if (nlhs == 2) {
        //-----
        // [VL,D] = dggevx(A,...)
        //-----
        jobl = 'V';
        jobr = 'N';
        plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
        VL = mxGetPr(plhs[0]);
        mflag = 1;
        Dindex = 1;

    } else if (nlhs == 3) {
        //-----
        // [VL,D,VR] = dggevx(A,...)
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
        mexErrMsgIdAndTxt("dggevx:TooManyOutput", "Too many output.");
    }

    if (nrhs == 3) {
        //-----
        // [...] = dggevx(A,B,'...')
        //-----
        if (mxIsChar(prhs[2])) {
            typelen = mxGetN(prhs[2]) + 1;
            type = mxMalloc(typelen);
            if (mxGetString(prhs[2], type, (mwSize) typelen) == 0) {
                if (strcmp(type, "matrix") == 0) mflag = 1;
                else if (strcmp(type, "vector") == 0) mflag = 0;
                else {
                    BALANC = mxMalloc(typelen);
                    mxGetString(prhs[2], BALANC, (mwSize) typelen);
                }
            }
        }
    }
    else if (nrhs == 4) {
        //-----
        // [...] = dggevx(A,B,'...','...')
        //-----
        if (mxIsChar(prhs[2])) {
            typelen = mxGetN(prhs[2]) + 1;
            type = mxMalloc(typelen);
            if (mxGetString(prhs[2], type, (mwSize) typelen) == 0) {
                if (strcmp(type, "matrix") == 0) mflag = 1;
                else if (strcmp(type, "vector") == 0) mflag = 0;
            }
        }
        if (mxIsChar(prhs[3])) {
            typelen = mxGetN(prhs[3]) + 1;
            BALANC = mxMalloc(typelen);
            mxGetString(prhs[3], BALANC, (mwSize) typelen);
        }
    }

    // Run
    run_vec(BALANC, &jobl, &jobr, A, B, C, D, VL, VR, ar, ai, b, &n);
    
    // check complex
    cflag = 0;
    for (i=n-1; i>=0; --i) {
        if (ai[i] != 0) {
            cflag = 1;
            break;
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

void run_vec(char *BALANC, char *jobl, char *jobr, double *A, double *B, double *C, double *D, double *VL, double *VR, double *ar, double *ai, double *b, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, info, ILO, IHI;
    char SENSE = 'N';
    double ABNRM, BBNRM;
    mxLogical BWORK;
    
    // Arrays
    mwSignedIndex *iwork;
	double *work, *LSCALE, *RSCALE, *RCONDE, *RCONDV;
    double dammy[1];

    iwork = MALLOC(mwSignedIndex, *n + 6);
    LSCALE = MALLOC(double, *n);
    RSCALE = MALLOC(double, *n);
    RCONDE = MALLOC(double, *n);
    RCONDV = MALLOC(double, *n);

    for (i=(*n)*(*n)-1; i>=0; --i) {
        C[i] = A[i];
        D[i] = B[i];
    }

    // Get the length
    dggevx(BALANC, jobl, jobr, &SENSE, n, C, n, D, n, ar, ai, b, VL, n, VR, n,
            &ILO, &IHI, LSCALE, RSCALE, &ABNRM, &BBNRM, RCONDE, RCONDV, 
            dammy, &lwork, iwork, &BWORK, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(double, lwork);

    // Run
    dggevx(BALANC, jobl, jobr, &SENSE, n, C, n, D, n, ar, ai, b, VL, n, VR, n,
            &ILO, &IHI, LSCALE, RSCALE, &ABNRM, &BBNRM, RCONDE, RCONDV, 
            work, &lwork, iwork, &BWORK, &info);
    
    // Release the workspace
    free(work);
    free(LSCALE);
    free(RSCALE);
    free(iwork);
    free(RCONDE);
    free(RCONDV);
}
