#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *jobs, float *A, float *B, float *eigenvalue, mwSignedIndex *n);
void run_mat(char *jobs, float *A, float *B, float *eigenvalue, mwSignedIndex *n);
void ssyevd(char *jobs, char *uplo, mwSignedIndex *n, float *A, mwSignedIndex *lda, float *W, float *work, mwSignedIndex *lwork, mwSignedIndex *iwork, mwSignedIndex *liwork, mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Scalars
	mwSignedIndex n, typelen;
    char jobs;
    char *type;

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
            // D = ssyevd(A)
            //-----
            // Return values
            plhs[0] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
            eigenvalue = (float *)mxGetPr(plhs[0]);

            // Run
            run_vec(&jobs, A, B, eigenvalue, &n);
        }
        else {
            if (mxIsChar(prhs[1])){
                typelen = mxGetN(prhs[1]) + 1;
                type = mxMalloc(typelen);

                if (mxGetString(prhs[1], type, (mwSize) typelen) == 0){
                    if (strcmp(type, "matrix") == 0){
                        //-----
                        // D = ssyevd(A,'matrix')
                        //-----
                        // Return values
                        plhs[0] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
                        eigenvalue = (float *)mxGetPr(plhs[0]);

                        // Run
                        run_mat(&jobs, A, B, eigenvalue, &n);
                    }
                    else{
                        mexErrMsgIdAndTxt("ssyevd:inputNotChar", "Input argument is not appropriate.");
                    }
                }
            }
            else{
                mexErrMsgIdAndTxt("ssyevd:inputNotChar", "Input argument must be of type char.");
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
            // [X,D] = dsyevd(A)
            //-----
            // Return values
            plhs[1] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
            eigenvalue = (float *)mxGetPr(plhs[1]);

            // Run
            run_mat(&jobs, A, B, eigenvalue, &n);
        }
        else {
            if (mxIsChar(prhs[1])){
                typelen = mxGetN(prhs[1]) + 1;
                type = mxMalloc(typelen);

                if (mxGetString(prhs[1], type, (mwSize) typelen) == 0){
                    if (strcmp(type, "vector") == 0){
                        //-----
                        // D = ssyevd(A,'vector')
                        //-----
                        // Return values
                        plhs[1] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
                        eigenvalue = (float *)mxGetPr(plhs[1]);

                        // Run
                        run_vec(&jobs, A, B, eigenvalue, &n);
                    }
                    else{
                        mexErrMsgIdAndTxt("ssyevd:inputNotChar", "Input argument is not appropriate.");
                    }
                }
            }
            else{
                mexErrMsgIdAndTxt("ssyevd:inputNotChar", "Input argument must be of type char.");
            }
        }
    }
}

void run_vec(char *jobs, float *A, float *B, float *eigenvalue, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, liwork=-1, info;
    char UPLO = 'U';
    
    // Arrays
	float *work;
    float dammy[1];
    mwSignedIndex *iwork;
    mwSignedIndex idammy[1];

    for (i=0; i<(*n)*(*n); i++) B[i] = A[i];

    // Get the length
    ssyevd(jobs, &UPLO, n, B, n, eigenvalue, dammy, &lwork, idammy, &liwork, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(float, lwork);
    liwork = idammy[0];
    iwork = MALLOC(mwSignedIndex, liwork);

    // Run
    ssyevd(jobs, &UPLO, n, B, n, eigenvalue, work, &lwork, iwork, &liwork, &info);
    
    // Release the workspace
    free(work);
    free(iwork);
}

void run_mat(char *jobs, float *A, float *B, float *eigenvalue, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, liwork=-1, info;
    char UPLO = 'U';
    
    // Arrays
	float *work, *eigtemp;
    float dammy[1];
    mwSignedIndex *iwork;
    mwSignedIndex idammy[1];

    eigtemp = MALLOC(float, *n);

    for (i=0; i<(*n)*(*n); i++) {
        B[i] = A[i];
        eigenvalue[i] = 0;
    }
    
    // Get the length
    ssyevd(jobs, &UPLO, n, B, n, eigtemp, dammy, &lwork, idammy, &liwork, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(float, lwork);
    liwork = idammy[0];
    iwork = MALLOC(mwSignedIndex, liwork);

    // Run
    ssyevd(jobs, &UPLO, n, B, n, eigtemp, work, &lwork, iwork, &liwork, &info);
    
    // Set the return values
    for (i = 0; i < *n; i++) eigenvalue[i*(*n) + i] = eigtemp[i];

    // Release the workspace
    free(eigtemp);
    free(work);
    free(iwork);
}