#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *jobvl, char *jobvr, float *A, float *B, float *C, float *D, float *VL, float *VR, float *ar, float *ai, float *b, mwSignedIndex *n);
void sggev3(char *JOBVL, char *JOBVR,
           mwSignedIndex *N, float *A, mwSignedIndex *LDA, float *B, mwSignedIndex *LDB,
           float *ALPHAR, float *ALPHAI, float *BETA, 
           float *VL, mwSignedIndex *LDVL, float *VR, mwSignedIndex *LDVR,
           float *WORK, mwSignedIndex *LWORK, mwSignedIndex *INFO);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 2)
        mexErrMsgIdAndTxt("sggev3:inputNotEnough", "Input argument is not enough.");

    // Scalars
	mwSignedIndex n, typelen, cflag, i, mflag, Dindex;
    char jobl, jobr;
    char *type;
    
    // Arrays
	float *A, *B, *C, *D, *VL, *VR, *eigenvalue, *ar, *ai, *b;

    // Arguments
	n = (mwSignedIndex)mxGetN(prhs[0]);
	A = (float *)mxGetPr(prhs[0]);
	B = (float *)mxGetPr(prhs[1]);
    C = MALLOC(float, n*n);
    D = MALLOC(float, n*n);
    ar = MALLOC(float, n);
    ai = MALLOC(float, n);
    b = MALLOC(float, n);

    if (nlhs == 1) {
        //-----
        // D = sggev3(A,B,...)
        //-----
        jobl = 'N';
        jobr = 'N';
        mflag = 0;
        Dindex = 0;

    } else if (nlhs == 2) {
        //-----
        // [VL,D] = sggev3(A,B,...)
        //-----
        jobl = 'V';
        jobr = 'N';
        plhs[0] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
        VL = (float *)mxGetPr(plhs[0]);
        mflag = 1;
        Dindex = 1;

    } else if (nlhs == 3) {
        //-----
        // [VL,D,VR] = sggev3(A,B,...)
        //-----
        jobl = 'V';
        jobr = 'V';
        plhs[0] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
        VL = (float *)mxGetPr(plhs[0]);
        plhs[2] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
        VR = (float *)mxGetPr(plhs[2]);
        mflag = 1;
        Dindex = 1;
        
    } else {
        free(C);
		free(D);
		free(ar);
		free(ai);
		free(b);
        mexErrMsgIdAndTxt("sggev3:TooManyOutput", "Too many output.");
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
        // [...] = sggev3(A,B,'...')
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
            float *eigenvalue;
            plhs[Dindex] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
            eigenvalue = (float *)mxGetPr(plhs[Dindex]);
            for (i=n-1; i>=0; --i) eigenvalue[i] = ar[i]/b[i];
        } else {
            mxComplexSingle *eigenvalue;
            plhs[Dindex] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxCOMPLEX);
            eigenvalue = mxGetComplexSingles(plhs[Dindex]);
            for (i=n-1; i>=0; --i) {
                eigenvalue[i].real = ar[i]/b[i];
                eigenvalue[i].imag = ai[i]/b[i];
            }
        }
    } else {
        if (cflag == 0) { // all eigenvalues are not complex number
            float *eigenvalue;
            plhs[Dindex] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
			eigenvalue = (float *)mxGetPr(plhs[Dindex]);
			for (i=n*n-1; i>=0; --i) eigenvalue[i] = 0;
			for (i=n-1; i>=0; --i) eigenvalue[i*n + i] = ar[i]/b[i];
		}
		else {
            mxComplexSingle *eigenvalue;
			plhs[Dindex] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxCOMPLEX);
			eigenvalue = mxGetComplexSingles(plhs[Dindex]);
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

void run_vec(char *jobvl, char *jobvr, float *A, float *B, float *C, float *D, float *VL, float *VR, float *ar, float *ai, float *b, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, info=0;
    
    // Arrays
	float *work;
    float dammy[1];

    for (i=(*n)*(*n)-1; i>=0; --i) {
        C[i] = A[i];
        D[i] = B[i];
    }

    // Get the length
    sggev3(jobvl, jobvr, n, C, n, D, n, ar, ai, b, VL, n, VR, n, dammy, &lwork, &info);
    
    // Set the workspace
    lwork=(mwSignedIndex)dammy[0];
    work = MALLOC(float, lwork);

    // Run
    sggev3(jobvl, jobvr, n, C, n, D, n, ar, ai, b, VL, n, VR, n, work, &lwork, &info);
    
    // Release the workspace
    free(work);
}
