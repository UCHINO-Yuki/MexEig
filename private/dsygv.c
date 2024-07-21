#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *jobs, double *A, double *B, double *C, double *D, double *eigenvalue, mwSignedIndex *n);
void run_mat(char *jobs, double *A, double *B, double *C, double *D, double *eigenvalue, mwSignedIndex *n);
void dsygv(mwSignedIndex *Itype, char *jobs, char *uplo, mwSignedIndex *n, double *A, mwSignedIndex *lda, double *B, mwSignedIndex *ldb, double *W, double *work, mwSignedIndex *lwork, mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 2)
        mexErrMsgIdAndTxt("dsygv:inputNotEnough", "Input argument is not enough.");

    // Scalars
	mwSignedIndex n, typelen;
    char jobs;
    char *type;
    
    // Arrays
	double *A, *B, *C, *D, *eigenvalue;

    // Arguments
	n = (mwSignedIndex)mxGetN(prhs[0]);
	A = mxGetPr(prhs[0]);
	B = mxGetPr(prhs[1]);
    D = MALLOC(double, n*n);

    if (nlhs == 1) {
        jobs = 'N';
        C = MALLOC(double, n*n);

        if (nrhs == 2) {
            //-----
            // D = dsygv(A,B)
            //-----
            // Return values
            plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
            eigenvalue = mxGetPr(plhs[0]);
            
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
                        // D = dsygv(A,B,'matrix')
                        //-----
                        // Return values
                        plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
                        eigenvalue = mxGetPr(plhs[0]);

                        // Run
                        run_mat(&jobs, A, B, C, D, eigenvalue, &n);
                    }
                    else{
                        free(C);
                        free(D);
                        mexErrMsgIdAndTxt("dsygv:inputNotChar", "Input argument is not appropriate.");
                    }
                }
            }
            else{
                free(C);
                free(D);
                mexErrMsgIdAndTxt("dsygv:inputNotChar", "Input argument must be of type char.");
            }
        }
        free(C);
    }
    else {
        jobs = 'V';
        
        // Return values
 	    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
        C = mxGetPr(plhs[0]);

        if (nrhs == 2){
            //-----
            // [X,D] = dsygv(A,B)
            //-----
            // Return values
            plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);
            eigenvalue = mxGetPr(plhs[1]);

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
                        // [X,D] = dsygv(A,B,'vector')
                        //-----
                        // Return values
                        plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
                        eigenvalue = mxGetPr(plhs[1]);

                        // Run
                        run_vec(&jobs, A, B, C, D, eigenvalue, &n);
                    }
                    else{
                        free(D);
                        mexErrMsgIdAndTxt("dsygv:inputNotChar", "Input argument is not appropriate.");
                    }
                }
            }
            else{
                free(D);
                mexErrMsgIdAndTxt("dsygv:inputNotChar", "Input argument must be of type char.");
            }
        }
    }
    free(D);
}

void run_vec(char *jobs, double *A, double *B, double *C, double *D, double *eigenvalue, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, info=0, Itype=1;
    char UPLO = 'U';
    
    // Arrays
	double *work;
    double dammy[1];

    for (i=(*n)*(*n)-1; i>=0; --i) {
        C[i] = A[i];
        D[i] = B[i];
    }

    // Get the length
    dsygv(&Itype, jobs, &UPLO, n, C, n, D, n, eigenvalue, dammy, &lwork, &info);
    
    // Set the workspace
    lwork=(mwSignedIndex)dammy[0];
    work = MALLOC(double, lwork);

    // Run
    dsygv(&Itype, jobs, &UPLO, n, C, n, D, n, eigenvalue, work, &lwork, &info);
    
    // Release the workspace
    free(work);
}

void run_mat(char *jobs, double *A, double *B, double *C, double *D, double *eigenvalue, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, lwork=-1, info, Itype=1;
    char UPLO = 'U';
    
    // Arrays
	double *work, *eigtemp;
    double dammy[1];

    eigtemp = MALLOC(double, *n);

    for (i=(*n)*(*n)-1; i>=0; i--) {
        C[i] = A[i];
        D[i] = B[i];
        eigenvalue[i] = 0;
    }
    
    // Get the length
    dsygv(&Itype, jobs, &UPLO, n, C, n, D, n, eigtemp, dammy, &lwork, &info);
    
    // Set the workspace
    lwork=(mwSignedIndex)dammy[0];
    work = MALLOC(double, lwork);

    // Run
    dsygv(&Itype, jobs, &UPLO, n, C, n, D, n, eigtemp, work, &lwork, &info);

    // Set the return values
    for (i = 0; i < *n; i++) eigenvalue[i*(*n) + i] = eigtemp[i];

    // Release the workspace
    free(eigtemp);
    free(work);
}
