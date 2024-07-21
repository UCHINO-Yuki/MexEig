#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *BALANC, char *jobl, char *jobr, float *A, float *B, float *C, float *D, float *VL, float *VR, float *ar, float *ai, float *b, mwSignedIndex *n);
void sggevx(char *BALANC, char *JOBVL, char *JOBVR, char *SENSE,
           mwSignedIndex *N, float *A, mwSignedIndex *LDA, float *B, mwSignedIndex *LDB,
           float *ALPHAR, float *ALPHAI, float *BETA, 
           float *VL, mwSignedIndex *LDVL, float *VR, mwSignedIndex *LDVR,
           mwSignedIndex *ILO, mwSignedIndex *IHI, float *LSCALE, float *RSCALE, 
           float *ABNRM, float *BBNRM, float *RCONDE, float *RCONDV, 
           float *WORK, mwSignedIndex *LWORK, mwSignedIndex *IWORK, mxLogical *BWORK, mwSignedIndex *INFO);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Scalars
	mwSignedIndex n, typelen, cflag, i, mflag, Dindex;
    char jobl, jobr, NN = 'N';
    char *type, *BALANC;
    BALANC = &NN;

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
        // D = sggevx(A,...)
        //-----
        jobl = 'N';
        jobr = 'N';
        mflag = 0;
        Dindex = 0;

    } else if (nlhs == 2) {
        //-----
        // [VL,D] = sggevx(A,...)
        //-----
        jobl = 'V';
        jobr = 'N';
        plhs[0] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
        VL = (float *)mxGetPr(plhs[0]);
        mflag = 1;
        Dindex = 1;

    } else if (nlhs == 3) {
        //-----
        // [VL,D,VR] = sggevx(A,...)
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
        mexErrMsgIdAndTxt("sggevx:TooManyOutput", "Too many output.");
    }

    if (nrhs == 3) {
        //-----
        // [...] = sggevx(A,B,'...')
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
        // [...] = sggevx(A,B,'...','...')
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

void run_vec(char *BALANC, char *jobl, char *jobr, float *A, float *B, float *C, float *D, float *VL, float *VR, float *ar, float *ai, float *b, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, info, ILO, IHI;
    char SENSE = 'N';
    float ABNRM, BBNRM;
    mxLogical BWORK;
    
    // Arrays
    mwSignedIndex *iwork;
	float *work, *LSCALE, *RSCALE, *RCONDE, *RCONDV;
    float dammy[1];

    iwork = MALLOC(mwSignedIndex, *n + 6);
    LSCALE = MALLOC(float, *n);
    RSCALE = MALLOC(float, *n);
    RCONDE = MALLOC(float, *n);
    RCONDV = MALLOC(float, *n);

    for (i=(*n)*(*n)-1; i>=0; --i) {
        C[i] = A[i];
        D[i] = B[i];
    }

    // Get the length
    sggevx(BALANC, jobl, jobr, &SENSE, n, C, n, D, n, ar, ai, b, VL, n, VR, n,
            &ILO, &IHI, LSCALE, RSCALE, &ABNRM, &BBNRM, RCONDE, RCONDV, 
            dammy, &lwork, iwork, &BWORK, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(float, lwork);

    // Run
    sggevx(BALANC, jobl, jobr, &SENSE, n, C, n, D, n, ar, ai, b, VL, n, VR, n,
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
