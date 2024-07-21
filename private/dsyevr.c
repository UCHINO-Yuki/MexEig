#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *jobs, double *A, double *B, double *eigenvalue, double *abstol, mwSignedIndex *n);
void run_mat(char *jobs, double *A, double *B, double *eigenvalue, double *abstol, mwSignedIndex *n);
void dsyevr(char *jobs, char *range, char *uplo, mwSignedIndex *n, double *A, mwSignedIndex *lda, double *VL, double *VU, mwSignedIndex *IL, mwSignedIndex *IU, double *ABSTOL, mwSignedIndex *m, double *W, double *Z, mwSignedIndex *ldz, mwSignedIndex *isuppz, double *work, mwSignedIndex *lwork, mwSignedIndex *iwork, mwSignedIndex *liwork, mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Scalars
	mwSignedIndex n, typelen;
    char jobs;
    char *type;
    double abstol;

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
            // D = dsyevr(A)
            //-----
            // Return values
            plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
            eigenvalue = mxGetPr(plhs[0]);

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
                    // D = dsyevr(A,'matrix')
                    //-----
                    // Return values
                    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
                    eigenvalue = mxGetPr(plhs[0]);

                    // Run
                    abstol = -1;
                    run_mat(&jobs, A, B, eigenvalue, &abstol, &n);
                }
                else{
                    mexErrMsgIdAndTxt("dsyevr:inputNotChar", "Input argument is not appropriate.");
                }
            }
            else{
                //-----
                // D = dsyevr(A,abstol)
                //-----
                // Return values
                plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
                eigenvalue = mxGetPr(plhs[0]);

                // Run
                abstol = mxGetScalar(prhs[1]);
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
                    // D = dsyevr(A,abstol,'matrix')
                    //-----
                    // Return values
                    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
                    eigenvalue = mxGetPr(plhs[0]);

                    // Run
                    abstol = mxGetScalar(prhs[1]);
                    run_mat(&jobs, A, B, eigenvalue, &abstol, &n);
                }
                else{
                    mexErrMsgIdAndTxt("dsyevr:inputNotChar", "Input argument is not appropriate.");
                }
            }
            else{
                mexErrMsgIdAndTxt("dsyevr:inputNotChar", "Input argument must be of type char.");
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
            // [X,D] = dsyevx(A)
            //-----
            // Return values
            plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);
            eigenvalue = mxGetPr(plhs[1]);

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
                    // [X,D] = dsyevr(A,'vector')
                    //-----
                    // Return values
                    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
                    eigenvalue = mxGetPr(plhs[1]);

                    // Run
                    abstol = -1;
                    run_vec(&jobs, A, B, eigenvalue, &abstol, &n);
                }
                else{
                    mexErrMsgIdAndTxt("dsyevr:inputNotChar", "Input argument is not appropriate.");
                }
            }
            else{
                //-----
                // [X,D] = dsyevr(A,abstol)
                //-----
                // Return values
                plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);
                eigenvalue = mxGetPr(plhs[1]);

                // Run
                abstol = mxGetScalar(prhs[1]);
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
                    // [X,D] = dsyevr(A,abstol,'vector')
                    //-----
                    // Return values
                    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
                    eigenvalue = mxGetPr(plhs[1]);

                    // Run
                    abstol = mxGetScalar(prhs[1]);
                    run_vec(&jobs, A, B, eigenvalue, &abstol, &n);
                }
                else{
                    mexErrMsgIdAndTxt("dsyevr:inputNotChar", "Input argument is not appropriate.");
                }
            }
            else{
                mexErrMsgIdAndTxt("dsyevr:inputNotChar", "Input argument must be of type char.");
            }
        }
    }
}

void run_vec(char *jobs, double *A, double *B, double *eigenvalue, double *abstol, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, m, IL, IU, lwork=-1, liwork=-1, info;
    double VL, VU;
    char UPLO = 'U', RANGE = 'A';
    
    // Arrays
	double *work, *C;
    double dammy[1];
    mwSignedIndex *iwork, *isuppz;
    mwSignedIndex idammy[1];
    
    // Set the workspace
    C = MALLOC(double, (*n)*(*n));
    isuppz = MALLOC(mwSignedIndex, 2*(*n));

    for (i=0; i<(*n)*(*n); i++) C[i] = A[i];

    // Get the length
    dsyevr(jobs, &RANGE, &UPLO, n, C, n, &VL, &VU, &IL, &IU, 
             abstol, &m, eigenvalue, B, n, isuppz, 
            dammy, &lwork, idammy, &liwork, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    liwork = idammy[0];
    work = MALLOC(double, lwork);
    iwork = MALLOC(mwSignedIndex, liwork);

    // Run
    dsyevr(jobs, &RANGE, &UPLO, n, C, n, &VL, &VU, &IL, &IU, 
             abstol, &m, eigenvalue, B, n, isuppz, 
            work, &lwork, iwork, &liwork, &info);

    // Release the workspace
    free(work);
    free(iwork);
    free(C);
    free(isuppz);
}

void run_mat(char *jobs, double *A, double *B, double *eigenvalue, double *abstol, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, m, IL, IU, lwork=-1, liwork=-1, info;
    double VL, VU;
    char UPLO = 'U', RANGE = 'A';
    
    // Arrays
	double *work, *C, *eigtemp;
    double dammy[1];
    mwSignedIndex *iwork, *isuppz;
    mwSignedIndex idammy[1];
    
    // Set the workspace
    C = MALLOC(double, (*n)*(*n));
    isuppz = MALLOC(mwSignedIndex, 2*(*n));
    eigtemp = MALLOC(double, *n);

    for (i=0; i<(*n)*(*n); i++) {
        C[i] = A[i];
        eigenvalue[i] = 0;
    }

    // Get the length
    dsyevr(jobs, &RANGE, &UPLO, n, C, n, &VL, &VU, &IL, &IU, 
             abstol, &m, eigtemp, B, n, isuppz, 
            dammy, &lwork, idammy, &liwork, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    liwork = idammy[0];
    work = MALLOC(double, lwork);
    iwork = MALLOC(mwSignedIndex, liwork);

    // Run
    dsyevr(jobs, &RANGE, &UPLO, n, C, n, &VL, &VU, &IL, &IU, 
             abstol, &m, eigtemp, B, n, isuppz, 
            work, &lwork, iwork, &liwork, &info);
    
    // Set the return values
    for (i = 0; i < *n; i++) eigenvalue[i*(*n) + i] = eigtemp[i];

    // Release the workspace
    free(work);
    free(iwork);
    free(C);
    free(isuppz);
    free(eigtemp);
}
