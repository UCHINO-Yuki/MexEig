#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *BALANC, char *jobl, char *jobr, double *A, double *B, double *VL, double *VR, double *eigr, double *eigi, mwSignedIndex *n);
void dgeevx(char *BALANC, char *JOBVL, char *JOBVR, char *SENSE, mwSignedIndex *N, double *A, mwSignedIndex *LDA, double *WR, double *WI,
            double *VL, mwSignedIndex *LDVL, double *VR, mwSignedIndex *LDVR, mwSignedIndex *ILO, mwSignedIndex *IHI, double *SCALE, double *ABNRM, 
            double *RCONDE, double *RCONDV, double *WORK, mwSignedIndex *LWORK, mwSignedIndex *IWORK, mwSignedIndex *INFO);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Scalars
	mwSignedIndex n, typelen, cflag, i, mflag, Dindex;
    char jobl, jobr, NN = 'N';
    char *type, *BALANC;
    BALANC = &NN;

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
        // D = dgeevx(A,...)
        //-----
        jobl = 'N';
        jobr = 'N';
        mflag = 0;
        Dindex = 0;

    } else if (nlhs == 2) {
        //-----
        // [VL,D] = dgeevx(A,...)
        //-----
        jobl = 'V';
        jobr = 'N';
        plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
        VL = mxGetPr(plhs[0]);
        mflag = 1;
        Dindex = 1;

    } else if (nlhs == 3) {
        //-----
        // [VL,D,VR] = dgeevx(A,...)
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
        mexErrMsgIdAndTxt("dgeevx:TooManyOutput", "Too many output.");
    }

    if (nrhs == 2) {
        if (mxIsChar(prhs[1])) {
            typelen = mxGetN(prhs[1]) + 1;
            type = mxMalloc(typelen);
            if (mxGetString(prhs[1], type, (mwSize) typelen) == 0) {
                if (strcmp(type, "matrix") == 0) mflag = 1;
                else if (strcmp(type, "vector") == 0) mflag = 0;
                else {
                    BALANC = mxMalloc(typelen);
                    mxGetString(prhs[1], BALANC, (mwSize) typelen);
                }
            }
        }
    }
    else if (nrhs == 3) {
        if (mxIsChar(prhs[1])) {
            typelen = mxGetN(prhs[1]) + 1;
            type = mxMalloc(typelen);
            if (mxGetString(prhs[1], type, (mwSize) typelen) == 0) {
                if (strcmp(type, "matrix") == 0) mflag = 1;
                else if (strcmp(type, "vector") == 0) mflag = 0;
            }
        }
        if (mxIsChar(prhs[2])) {
            typelen = mxGetN(prhs[2]) + 1;
            BALANC = mxMalloc(typelen);
            mxGetString(prhs[2], BALANC, (mwSize) typelen);
        }
    }

    // Run
    run_vec(BALANC, &jobl, &jobr, A, B, VL, VR, eigr, eigi, &n);
    
    // check complex
    cflag = 0;
    for (i=n-1; i>=0; --i) {
        if (eigi[i] != 0) {
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

void run_vec(char *BALANC, char *jobl, char *jobr, double *A, double *B, double *VL, double *VR, double *eigr, double *eigi, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, info, ILO, IHI;
    char SENSE = 'N';
    double ABNRM;
    
    // Arrays
    mwSignedIndex *iwork;
	double *work, *SCALE, *RCONDE, *RCONDV;
    double dammy[1];

    iwork = MALLOC(mwSignedIndex, 2*(*n)-1);
    SCALE = MALLOC(double, *n);
    RCONDE = MALLOC(double, *n);
    RCONDV = MALLOC(double, *n);
    for (i=0; i<(*n)*(*n); i++) B[i] = A[i];

    // Get the length
    dgeevx(BALANC, jobl, jobr, &SENSE, n, B, n, eigr, eigi, VL, n, VR, n, 
            &ILO, &IHI, SCALE, &ABNRM, RCONDE, RCONDV, dammy, &lwork, iwork, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(double, lwork);

    // Run
    dgeevx(BALANC, jobl, jobr, &SENSE, n, B, n, eigr, eigi, VL, n, VR, n, 
            &ILO, &IHI, SCALE, &ABNRM, RCONDE, RCONDV, work, &lwork, iwork, &info);
    
    // Release the workspace
    free(work);
    free(SCALE);
    free(iwork);
    free(RCONDE);
    free(RCONDV);
}
