#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *BALANC, char *jobl, char *jobr, float *A, float *B, float *VL, float *VR, float *eigr, float *eigi, mwSignedIndex *n);
void sgeevx(char *BALANC, char *JOBVL, char *JOBVR, char *SENSE, mwSignedIndex *N, float *A, mwSignedIndex *LDA, float *WR, float *WI,
            float *VL, mwSignedIndex *LDVL, float *VR, mwSignedIndex *LDVR, mwSignedIndex *ILO, mwSignedIndex *IHI, float *SCALE, float *ABNRM, 
            float *RCONDE, float *RCONDV, float *WORK, mwSignedIndex *LWORK, mwSignedIndex *IWORK, mwSignedIndex *INFO);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Scalars
	mwSignedIndex n, typelen, cflag, i, mflag, Dindex;
    char jobl, jobr, NN = 'N';
    char *type, *BALANC;
    BALANC = &NN;

    // Arrays
	float *A, *B, *VL, *VR, *eigr, *eigi;

    // Arguments
	n = mxGetN(prhs[0]);
	A = (float *)mxGetPr(prhs[0]);
    B = MALLOC(float, n*n);
    eigr = MALLOC(float, n);
    eigi = MALLOC(float, n);

    if (nlhs == 1) {
        //-----
        // D = sgeevx(A,...)
        //-----
        jobl = 'N';
        jobr = 'N';
        mflag = 0;
        Dindex = 0;

    } else if (nlhs == 2) {
        //-----
        // [VL,D] = sgeevx(A,...)
        //-----
        jobl = 'V';
        jobr = 'N';
        plhs[0] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
        VL = (float *)mxGetPr(plhs[0]);
        mflag = 1;
        Dindex = 1;

    } else if (nlhs == 3) {
        //-----
        // [VL,D,VR] = sgeevx(A,...)
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
		free(B);
		free(eigr);
		free(eigi);
        mexErrMsgIdAndTxt("sgeevx:TooManyOutput", "Too many output.");
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
            float *eigenvalue;
            plhs[Dindex] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
            eigenvalue = (float *)mxGetPr(plhs[Dindex]);
            for (i=n-1; i>=0; --i) eigenvalue[i] = eigr[i];
        } else {
            mxComplexSingle *eigenvalue;
            plhs[Dindex] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxCOMPLEX);
            eigenvalue = mxGetComplexSingles(plhs[Dindex]);
            for (i=n-1; i>=0; --i) {
                eigenvalue[i].real = eigr[i];
                eigenvalue[i].imag = eigi[i];
            }
        }
    } else {
        if (cflag == 0) { // all eigenvalues are not complex number
            float *eigenvalue;
            plhs[Dindex] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
			eigenvalue = (float *)mxGetPr(plhs[Dindex]);
			for (i=n*n-1; i>=0; --i) eigenvalue[i] = 0;
			for (i=n-1; i>=0; --i) eigenvalue[i*n + i] = eigr[i];
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
				eigenvalue[i*n + i].real = eigr[i];
				eigenvalue[i*n + i].imag = eigi[i];
			}
		}
    } 

	free(B);
	free(eigr);
	free(eigi);
}

void run_vec(char *BALANC, char *jobl, char *jobr, float *A, float *B, float *VL, float *VR, float *eigr, float *eigi, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, info, ILO, IHI;
    char SENSE = 'N';
    float ABNRM;
    
    // Arrays
    mwSignedIndex *iwork;
	float *work, *SCALE, *RCONDE, *RCONDV;
    float dammy[1];

    iwork = MALLOC(mwSignedIndex, 2*(*n)-1);
    SCALE = MALLOC(float, *n);
    RCONDE = MALLOC(float, *n);
    RCONDV = MALLOC(float, *n);
    for (i=0; i<(*n)*(*n); i++) B[i] = A[i];

    // Get the length
    sgeevx(BALANC, jobl, jobr, &SENSE, n, B, n, eigr, eigi, VL, n, VR, n, 
            &ILO, &IHI, SCALE, &ABNRM, RCONDE, RCONDV, dammy, &lwork, iwork, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(float, lwork);

    // Run
    sgeevx(BALANC, jobl, jobr, &SENSE, n, B, n, eigr, eigi, VL, n, VR, n, 
            &ILO, &IHI, SCALE, &ABNRM, RCONDE, RCONDV, work, &lwork, iwork, &info);
    
    // Release the workspace
    free(work);
    free(SCALE);
    free(iwork);
    free(RCONDE);
    free(RCONDV);
}
