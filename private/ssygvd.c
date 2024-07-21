#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *jobs, float *A, float *B, float *C, float *D, float *eigenvalue, mwSignedIndex *n);
void run_mat(char *jobs, float *A, float *B, float *C, float *D, float *eigenvalue, mwSignedIndex *n);
void ssygvd(mwSignedIndex *Itype, char *jobs, char *uplo, mwSignedIndex *n, float *A, mwSignedIndex *lda, float *B, mwSignedIndex *ldb, float *W, float *work, mwSignedIndex *lwork, mwSignedIndex *iwork, mwSignedIndex *liwork, mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 2)
        mexErrMsgIdAndTxt("ssygvd:inputNotEnough", "Input argument is not enough.");

    // Scalars
	mwSignedIndex n, typelen;
    char jobs;
    char *type;
    
    // Arrays
	float *A, *B, *C, *D, *eigenvalue;

    // Arguments
	n = (mwSignedIndex)mxGetN(prhs[0]);
	A = (float *)mxGetPr(prhs[0]);
	B = (float *)mxGetPr(prhs[1]);
    D = MALLOC(float, n*n);

    if (nlhs == 1) {
        jobs = 'N';
        C = MALLOC(float, n*n);

        if (nrhs == 2) {
            //-----
            // D = ssygvd(A,B)
            //-----
            // Return values
            plhs[0] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
            eigenvalue = (float *)mxGetPr(plhs[0]);
            
            // Run
            run_vec(&jobs, A, B, C, D, eigenvalue, &n);
        }
        else {
            if (mxIsChar(prhs[2])){
                typelen = mxGetN(prhs[2]) + 1;
                type = mxMalloc(typelen);

                if (mxGetString(prhs[2], type, (mwSize) typelen) == 0){
                    if (strcmp(type, "matrix") == 0){
                        //-----
                        // D = ssygvd(A,B,'matrix')
                        //-----
                        // Return values
                        plhs[0] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
                        eigenvalue = (float *)mxGetPr(plhs[0]);

                        // Run
                        run_mat(&jobs, A, B, C, D, eigenvalue, &n);
                    }
                    else{
                        free(C);
                        free(D);
                        mexErrMsgIdAndTxt("ssygvd:inputNotChar", "Input argument is not appropriate.");
                    }
                }
            }
            else{
                free(C);
                free(D);
                mexErrMsgIdAndTxt("ssygvd:inputNotChar", "Input argument must be of type char.");
            }
        }
        free(C);
    }
    else {
        jobs = 'V';
        
        // Return values
 	    plhs[0] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
        C = (float *)mxGetPr(plhs[0]);

        if (nrhs == 2){
            //-----
            // [X,D] = ssygvd(A,B)
            //-----
            // Return values
            plhs[1] = mxCreateNumericMatrix(n, n, mxSINGLE_CLASS, mxREAL);
            eigenvalue = (float *)mxGetPr(plhs[1]);

            // Run
            run_mat(&jobs, A, B, C, D, eigenvalue, &n);
        }
        else {
            if (mxIsChar(prhs[2])){
                typelen = mxGetN(prhs[2]) + 1;
                type = mxMalloc(typelen);

                if (mxGetString(prhs[2], type, (mwSize) typelen) == 0){
                    if (strcmp(type, "vector") == 0){
                        //-----
                        // [X,D] = ssygvd(A,B,'vector')
                        //-----
                        // Return values
                        plhs[1] = mxCreateNumericArray(1, &n, mxSINGLE_CLASS, mxREAL);
                        eigenvalue = (float *)mxGetPr(plhs[1]);

                        // Run
                        run_vec(&jobs, A, B, C, D, eigenvalue, &n);
                    }
                    else{
                        free(D);
                        mexErrMsgIdAndTxt("ssygvd:inputNotChar", "Input argument is not appropriate.");
                    }
                }
            }
            else{
                free(D);
                mexErrMsgIdAndTxt("ssygvd:inputNotChar", "Input argument must be of type char.");
            }
        }
    }
    free(D);
}

void run_vec(char *jobs, float *A, float *B, float *C, float *D, float *eigenvalue, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, liwork=-1, info=0, Itype=1;
    char UPLO = 'U';
    
    // Arrays
	float *work;
    float dammy[1];
    mwSignedIndex *iwork;
    mwSignedIndex idammy[1];

    for (i=(*n)*(*n)-1; i>=0; --i) {
        C[i] = A[i];
        D[i] = B[i];
    }

    // Get the length
    ssygvd(&Itype, jobs, &UPLO, n, C, n, D, n, eigenvalue, dammy, &lwork, idammy, &liwork, &info);
    
    // Set the workspace
    lwork=(mwSignedIndex)dammy[0];
    work = MALLOC(float, lwork);
    liwork = idammy[0];
    iwork = MALLOC(mwSignedIndex, liwork);

    // Run
    ssygvd(&Itype, jobs, &UPLO, n, C, n, D, n, eigenvalue, work, &lwork, iwork, &liwork, &info);
    
    // Release the workspace
    free(work);
    free(iwork);
}

void run_mat(char *jobs, float *A, float *B, float *C, float *D, float *eigenvalue, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, liwork=-1, info, Itype=1;
    char UPLO = 'U';
    
    // Arrays
	float *work, *eigtemp;
    float dammy[1];
    mwSignedIndex *iwork;
    mwSignedIndex idammy[1];

    eigtemp = MALLOC(float, *n);

    for (i=(*n)*(*n)-1; i>=0; i--) {
        C[i] = A[i];
        D[i] = B[i];
        eigenvalue[i] = 0;
    }
    
    // Get the length
    ssygvd(&Itype, jobs, &UPLO, n, C, n, D, n, eigtemp, dammy, &lwork, idammy, &liwork, &info);
    
    // Set the workspace
    lwork=(mwSignedIndex)dammy[0];
    work = MALLOC(float, lwork);
    liwork = idammy[0];
    iwork = MALLOC(mwSignedIndex, liwork);

    // Run
    ssygvd(&Itype, jobs, &UPLO, n, C, n, D, n, eigtemp, work, &lwork, iwork, &liwork, &info);

    // Set the return values
    for (i = 0; i < *n; i++) eigenvalue[i*(*n) + i] = eigtemp[i];

    // Release the workspace
    free(eigtemp);
    free(work);
    free(iwork);
}
