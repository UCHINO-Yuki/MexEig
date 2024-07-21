#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *jobs, float *A, float *B, float *eigenvalue, float *abstol, mwSignedIndex *n);
void run_mat(char *jobs, float *A, float *B, float *eigenvalue, float *abstol, mwSignedIndex *n);
void ssyevx(char *jobs, char *range, char *uplo, mwSignedIndex *n, float *A, mwSignedIndex *lda, float *VL, float *VU, mwSignedIndex *IL, mwSignedIndex *IU, float *ABSTOL, mwSignedIndex *m, float *W, float *Z, mwSignedIndex *ldz, float *work, mwSignedIndex *lwork, mwSignedIndex *iwork, mwSignedIndex *ifail, mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Scalars
	mwSignedIndex n, typelen;
    char jobs;
    char *type;
    float abstol;

    // Arrays
	float *A, *B, *eigenvalue;

    // Arguments
	n = mxGetN(prhs[0]);
	A = (float *)mxGetPr(prhs[0]);
    
    if (nlhs == 1) {
        jobs = 'N';
        B = MALLOC(float, n*n);
        
        if (nrhs == 1) {
            //-----
            // D = ssyevx(A)
            //-----
            // Return values
            plhs[0] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
            eigenvalue = (float *)mxGetPr(plhs[0]);

            // Run
            abstol = -1;
            run_vec(&jobs, A, B, eigenvalue, &abstol, &n);
        }
        else if (nrhs == 2) {
            if (mxIsChar(prhs[1])){
                typelen = mxGetN(prhs[1]) + 1;
                type = mxMalloc(typelen);
                mxGetString(prhs[1], type, (mwSize) typelen);

                if (strcmp(type, "matrix") == 0){
                    //-----
                    // D = ssyevx(A,'matrix')
                    //-----
                    // Return values
                    plhs[0] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
                    eigenvalue = (float *)mxGetPr(plhs[0]);

                    // Run
                    abstol = -1;
                    run_mat(&jobs, A, B, eigenvalue, &abstol, &n);
                }
                else{
                    mexErrMsgIdAndTxt("ssyevx:inputNotChar", "Input argument is not appropriate.");
                }
            }
            else{
                //-----
                // D = ssyevx(A,abstol)
                //-----
                // Return values
                plhs[0] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
                eigenvalue = (float *)mxGetPr(plhs[0]);

                // Run
                abstol = (float)mxGetScalar(prhs[1]);
                run_vec(&jobs, A, B, eigenvalue, &abstol, &n);
            }
        }
        else {
            if (mxIsChar(prhs[2])){
                typelen = mxGetN(prhs[2]) + 1;
                type = mxMalloc(typelen);
                mxGetString(prhs[2], type, (mwSize) typelen);

                if (strcmp(type, "matrix") == 0){
                    //-----
                    // D = ssyevx(A,abstol,'matrix')
                    //-----
                    // Return values
                    plhs[0] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
                    eigenvalue = (float *)mxGetPr(plhs[0]);

                    // Run
                    abstol = (float)mxGetScalar(prhs[1]);
                    run_mat(&jobs, A, B, eigenvalue, &abstol, &n);
                }
                else{
                    mexErrMsgIdAndTxt("ssyevx:inputNotChar", "Input argument is not appropriate.");
                }
            }
            else{
                mexErrMsgIdAndTxt("ssyevx:inputNotChar", "Input argument must be of type char.");
            }
        }
        free(B);
    }
    else {
        jobs = 'V';
        
        // Return values
 	    plhs[0] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
        B = (float *)mxGetPr(plhs[0]);

        if (nrhs == 1){
            //-----
            // [X,D] = ssyevx(A)
            //-----
            // Return values
            plhs[1] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
            eigenvalue = (float *)mxGetPr(plhs[1]);

            // Run
            abstol = -1;
            run_mat(&jobs, A, B, eigenvalue, &abstol, &n);
        }
        else if (nrhs == 2) {
            if (mxIsChar(prhs[1])){
                typelen = mxGetN(prhs[1]) + 1;
                type = mxMalloc(typelen);
                mxGetString(prhs[1], type, (mwSize) typelen);

                if (strcmp(type, "vector") == 0){
                    //-----
                    // [X,D] = ssyevx(A,'vector')
                    //-----
                    // Return values
                    plhs[1] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
                    eigenvalue = (float *)mxGetPr(plhs[1]);

                    // Run
                    abstol = -1;
                    run_vec(&jobs, A, B, eigenvalue, &abstol, &n);
                }
                else{
                    mexErrMsgIdAndTxt("ssyevx:inputNotChar", "Input argument is not appropriate.");
                }
            }
            else{
                //-----
                // [X,D] = ssyevx(A,abstol)
                //-----
                // Return values
                plhs[1] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
                eigenvalue = mxGetPr(plhs[1]);

                // Run
                abstol = (float)mxGetScalar(prhs[1]);
                run_mat(&jobs, A, B, eigenvalue, &abstol, &n);
            }
        }
        else {
            if (mxIsChar(prhs[2])){
                typelen = mxGetN(prhs[2]) + 1;
                type = mxMalloc(typelen);
                mxGetString(prhs[2], type, (mwSize) typelen);

                if (strcmp(type, "vector") == 0){
                    //-----
                    // [X,D] = ssyevx(A,abstol,'vector')
                    //-----
                    // Return values
                    plhs[1] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
                    eigenvalue = mxGetPr(plhs[1]);

                    // Run
                    abstol = (float)mxGetScalar(prhs[1]);
                    run_vec(&jobs, A, B, eigenvalue, &abstol, &n);
                }
                else{
                    mexErrMsgIdAndTxt("ssyevx:inputNotChar", "Input argument is not appropriate.");
                }
            }
            else{
                mexErrMsgIdAndTxt("ssyevx:inputNotChar", "Input argument must be of type char.");
            }
        }
    }
}


void run_vec(char *jobs, float *A, float *B, float *eigenvalue, float *abstol, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, m, IL, IU, lwork=-1, info;
    float VL, VU;
    char UPLO = 'U', RANGE = 'A';
    
    // Arrays
	float *work, *C;
    float dammy[1];
    mwSignedIndex *iwork, *ifail;
    
    // Set the workspace
    C = MALLOC(float, (*n)*(*n));
    ifail = MALLOC(mwSignedIndex, *n);
    iwork = MALLOC(mwSignedIndex, 5*(*n));

    for (i=0; i<(*n)*(*n); i++) C[i] = A[i];

    // Get the length
    ssyevx(jobs, &RANGE, &UPLO, n, C, n, &VL, &VU, &IL, &IU, 
            abstol, &m, eigenvalue, B, n, dammy, &lwork, iwork, ifail, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(float, lwork);

    // Run
    ssyevx(jobs, &RANGE, &UPLO, n, C, n, &VL, &VU, &IL, &IU, 
            abstol, &m, eigenvalue, B, n, work, &lwork, iwork, ifail, &info);

    // Release the workspace
    free(work);
    free(iwork);
    free(C);
    free(ifail);
}

void run_mat(char *jobs, float *A, float *B, float *eigenvalue, float *abstol, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, m, IL, IU, lwork=-1, info;
    float VL, VU;
    char UPLO = 'U', RANGE = 'A';
    
    // Arrays
	float *work, *C, *eigtemp;
    float dammy[1];
    mwSignedIndex *iwork, *ifail;
    
    // Set the workspace
    C = MALLOC(float, (*n)*(*n));
    ifail = MALLOC(mwSignedIndex, *n);
    iwork = MALLOC(mwSignedIndex, 5*(*n));
    eigtemp = MALLOC(float, *n);

    for (i=0; i<(*n)*(*n); i++) {
        C[i] = A[i];
        eigenvalue[i] = 0;
    }

    // Get the length
    ssyevx(jobs, &RANGE, &UPLO, n, C, n, &VL, &VU, &IL, &IU, 
            abstol, &m, eigtemp, B, n, dammy, &lwork, iwork, ifail, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(float, lwork);

    // Run
    ssyevx(jobs, &RANGE, &UPLO, n, C, n, &VL, &VU, &IL, &IU, 
            abstol, &m, eigtemp, B, n, work, &lwork, iwork, ifail, &info);
    
    // Set the return values
    for (i = 0; i < *n; i++) eigenvalue[i*(*n) + i] = eigtemp[i];

    // Release the workspace
    free(work);
    free(iwork);
    free(C);
    free(ifail);
    free(eigtemp);
}
