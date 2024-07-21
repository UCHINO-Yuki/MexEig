#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *jobs, double *A, double *B, double *eigenvalue, mwSignedIndex *n);
void run_mat(char *jobs, double *A, double *B, double *eigenvalue, mwSignedIndex *n);
void dsyevd(char *jobs, char *uplo, mwSignedIndex *n, double *A, mwSignedIndex *lda, double *W, double *work, mwSignedIndex *lwork, mwSignedIndex *iwork, mwSignedIndex *liwork, mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Scalars
	mwSignedIndex n, typelen;
    char jobs;
    char *type;

    // Arrays
	double *A, *B, *eigenvalue;

    // Arguments
	n = mxGetN(prhs[0]);
	A = mxGetPr(prhs[0]);

    if (nlhs == 1) {
        jobs = 'N';
        B = MALLOC(double, n*n);

        if (nrhs == 1) {
            //-----
            // D = dsyevd(A)
            //-----
            // Return values
            plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
            eigenvalue = mxGetPr(plhs[0]);

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
                        // D = dsyevd(A,'matrix')
                        //-----
                        // Return values
                        plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
                        eigenvalue = mxGetPr(plhs[0]);

                        // Run
                        run_mat(&jobs, A, B, eigenvalue, &n);
                    }
                    else{
                        mexErrMsgIdAndTxt("dsyevd:inputNotChar", "Input argument is not appropriate.");
                    }
                }
            }
            else{
                mexErrMsgIdAndTxt("dsyevd:inputNotChar", "Input argument must be of type char.");
            }
        }

        free(B);
    }
    else {
        jobs = 'V';
        
        // Return values
 	    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
        B = mxGetPr(plhs[0]);

        if (nrhs == 1){
            //-----
            // [X,D] = dsyevd(A)
            //-----
            // Return values
            plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);
            eigenvalue = mxGetPr(plhs[1]);

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
                        // D = dsyevd(A,'vector')
                        //-----
                        // Return values
                        plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
                        eigenvalue = mxGetPr(plhs[1]);

                        // Run
                        run_vec(&jobs, A, B, eigenvalue, &n);
                    }
                    else{
                        mexErrMsgIdAndTxt("dsyevd:inputNotChar", "Input argument is not appropriate.");
                    }
                }
            }
            else{
                mexErrMsgIdAndTxt("dsyevd:inputNotChar", "Input argument must be of type char.");
            }
        }
    }
}

void run_vec(char *jobs, double *A, double *B, double *eigenvalue, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, liwork=-1, info;
    char UPLO = 'U';
    
    // Arrays
	double *work;
    double dammy[1];
    mwSignedIndex *iwork;
    mwSignedIndex idammy[1];

    for (i=0; i<(*n)*(*n); i++) B[i] = A[i];

    // Get the length
    dsyevd(jobs, &UPLO, n, B, n, eigenvalue, dammy, &lwork, idammy, &liwork, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(double, lwork);
    liwork = idammy[0];
    iwork = MALLOC(mwSignedIndex, liwork);

    // Run
    dsyevd(jobs, &UPLO, n, B, n, eigenvalue, work, &lwork, iwork, &liwork, &info);
    
    // Release the workspace
    free(work);
    free(iwork);
}

void run_mat(char *jobs, double *A, double *B, double *eigenvalue, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, liwork=-1, info;
    char UPLO = 'U';
    
    // Arrays
	double *work, *eigtemp;
    double dammy[1];
    mwSignedIndex *iwork;
    mwSignedIndex idammy[1];

    eigtemp = MALLOC(double, *n);

    for (i=0; i<(*n)*(*n); i++) {
        B[i] = A[i];
        eigenvalue[i] = 0;
    }
    
    // Get the length
    dsyevd(jobs, &UPLO, n, B, n, eigtemp, dammy, &lwork, idammy, &liwork, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(double, lwork);
    liwork = idammy[0];
    iwork = MALLOC(mwSignedIndex, liwork);

    // Run
    dsyevd(jobs, &UPLO, n, B, n, eigtemp, work, &lwork, iwork, &liwork, &info);
    
    // Set the return values
    for (i = 0; i < *n; i++) eigenvalue[i*(*n) + i] = eigtemp[i];

    // Release the workspace
    free(eigtemp);
    free(work);
    free(iwork);
}