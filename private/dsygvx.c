#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define MALLOC(type, n) ((type *) malloc(sizeof(type) * n))

void run_vec(char *jobs, double *A, double *B, double *C, double *D, double *eigenvalue, double *abstol, mwSignedIndex *n);
void run_mat(char *jobs, double *A, double *B, double *C, double *D, double *eigenvalue, double *abstol, mwSignedIndex *n);
void dsygvx(mwSignedIndex *Itype, char *jobs, char *range, char *uplo, 
            mwSignedIndex *n, double *A, mwSignedIndex *lda, double *B, mwSignedIndex *ldb, 
            double *VL, double *VU, mwSignedIndex *IL, mwSignedIndex *IU, 
            double *ABSTOL, mwSignedIndex *m, double *W, double *Z, mwSignedIndex *ldz, 
            double *work, mwSignedIndex *lwork, mwSignedIndex *iwork, mwSignedIndex *ifail, mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 2)
        mexErrMsgIdAndTxt("dsygvx:inputNotEnough", "Input argument is not enough.");

    // Scalars
	mwSignedIndex n, typelen;
    char jobs;
    char *type;
    double abstol;
    
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
            // D = dsygvx(A,B)
            //-----
            // Return values
            plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
            eigenvalue = mxGetPr(plhs[0]);
            
            // Run
            abstol = -1;
            run_vec(&jobs, A, B, C, D, eigenvalue, &abstol, &n);
        }
        else if (nrhs == 3){
            if (mxIsChar(prhs[2])){
                typelen = mxGetN(prhs[2]) + 1;
                type = mxMalloc(typelen);

                if (mxGetString(prhs[2], type, (mwSize) typelen) == 0){
                    if (strcmp(type, "matrix") == 0){
                        //-----
                        // D = dsygvx(A,B,'matrix')
                        //-----
                        // Return values
                        plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
                        eigenvalue = mxGetPr(plhs[0]);

                        // Run
                        abstol = -1;
                        run_mat(&jobs, A, B, C, D, eigenvalue, &abstol, &n);
                    }
                    else{
                        free(C);
                        free(D);
                        mexErrMsgIdAndTxt("dsygvx:inputNotChar", "Input argument is not appropriate.");
                    }
                }
            }
            else{
                //-----
                // D = dsygvx(A,B,abstol)
                //-----
                // Return values
                plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
                eigenvalue = mxGetPr(plhs[0]);
                
                // Run
                abstol = mxGetScalar(prhs[2]);
                run_vec(&jobs, A, B, C, D, eigenvalue, &abstol, &n);
            }
        }
        else{
            if (mxIsChar(prhs[3])){
                typelen = mxGetN(prhs[3]) + 1;
                type = mxMalloc(typelen);

                if (mxGetString(prhs[3], type, (mwSize) typelen) == 0){
                    if (strcmp(type, "matrix") == 0){
                        //-----
                        // D = dsygvx(A,B,abstol,'matrix')
                        //-----
                        // Return values
                        plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
                        eigenvalue = mxGetPr(plhs[0]);

                        // Run
                        abstol = mxGetScalar(prhs[2]);
                        run_mat(&jobs, A, B, C, D, eigenvalue, &abstol, &n);
                    }
                    else{
                        free(C);
                        free(D);
                        mexErrMsgIdAndTxt("dsygvx:inputNotChar", "Input argument is not appropriate.");
                    }
                }
            }
            else{
                free(C);
                free(D);
                mexErrMsgIdAndTxt("dsygvx:inputNotChar", "Input argument is not appropriate.");
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
            // [X,D] = dsygvx(A,B)
            //-----
            // Return values
            plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);
            eigenvalue = mxGetPr(plhs[1]);

            // Run
            abstol = -1;
            run_mat(&jobs, A, B, C, D, eigenvalue, &abstol, &n);
        }
        else if (nrhs == 3){
            if (mxIsChar(prhs[2])){
                typelen = mxGetN(prhs[2]) + 1;
                type = mxMalloc(typelen);

                if (mxGetString(prhs[2], type, (mwSize) typelen) == 0){
                    if (strcmp(type, "vector") == 0){
                        //-----
                        // [X,D] = dsygvx(A,B,'vector')
                        //-----
                        // Return values
                        plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
                        eigenvalue = mxGetPr(plhs[1]);

                        // Run
                        abstol = -1;
                        run_vec(&jobs, A, B, C, D, eigenvalue, &abstol, &n);
                    }
                    else{
                        free(D);
                        mexErrMsgIdAndTxt("dsygvx:inputNotChar", "Input argument is not appropriate.");
                    }
                }
            }
            else{
                //-----
                // [X,D] = dsygvx(A,B,abstol)
                //-----
                // Return values
                plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);
                eigenvalue = mxGetPr(plhs[1]);
                
                // Run
                abstol = mxGetScalar(prhs[2]);
                run_mat(&jobs, A, B, C, D, eigenvalue, &abstol, &n);
            }
        }
        else{
            if (mxIsChar(prhs[3])){
                typelen = mxGetN(prhs[3]) + 1;
                type = mxMalloc(typelen);

                if (mxGetString(prhs[3], type, (mwSize) typelen) == 0){
                    if (strcmp(type, "vector") == 0){
                        //-----
                        // [X,D] = dsygvx(A,B,abstol,'vector')
                        //-----
                        // Return values
                        plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
                        eigenvalue = mxGetPr(plhs[1]);

                        // Run
                        abstol = mxGetScalar(prhs[2]);
                        run_vec(&jobs, A, B, C, D, eigenvalue, &abstol, &n);
                    }
                    else{
                        free(C);
                        free(D);
                        mexErrMsgIdAndTxt("dsygvx:inputNotChar", "Input argument is not appropriate.");
                    }
                }
            }
            else{
                free(C);
                free(D);
                mexErrMsgIdAndTxt("dsygvx:inputNotChar", "Input argument is not appropriate.");
            }
        }
    }
    free(D);
}

void run_vec(char *jobs, double *A, double *B, double *C, double *D, double *eigenvalue, double *abstol, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, m, IL, IU, lwork=-1, info=0, Itype=1;
    double VL, VU;
    char UPLO = 'U', RANGE = 'A';
    
    // Arrays
	double *work, *E;
    double dammy[1];
    mwSignedIndex *iwork, *ifail;
    ifail = MALLOC(mwSignedIndex, *n);
    iwork = MALLOC(mwSignedIndex, 5*(*n));
    E = MALLOC(double, (*n)*(*n));

    for (i=(*n)*(*n)-1; i>=0; --i) {
        E[i] = A[i];
        D[i] = B[i];
    }

    // Get the length
    dsygvx(&Itype, jobs, &RANGE, &UPLO, n, E, n, D, n, &VL, &VU, &IL, &IU, 
            abstol, &m, eigenvalue, C, n, dammy, &lwork, iwork, ifail, &info);
    
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(double, lwork);

    // Run
    dsygvx(&Itype, jobs, &RANGE, &UPLO, n, E, n, D, n, &VL, &VU, &IL, &IU, 
            abstol, &m, eigenvalue, C, n, work, &lwork, iwork, ifail, &info);
    
    // Release the workspace
    free(work);
    free(iwork);
    free(E);
}

void run_mat(char *jobs, double *A, double *B, double *C, double *D, double *eigenvalue, double *abstol, mwSignedIndex *n) {
    // Scalars
    mwSignedIndex i, m, IL, IU, lwork=-1, info=0, Itype=1;
    double VL, VU;
    char UPLO = 'U', RANGE = 'A';
    
    // Arrays
	double *work, *eigtemp, *E;
    mwSignedIndex *iwork, *ifail;
    double dammy[1];
    ifail = MALLOC(mwSignedIndex, *n);
    iwork = MALLOC(mwSignedIndex, 5*(*n));
    E = MALLOC(double, (*n)*(*n));

    eigtemp = MALLOC(double, *n);

    for (i=(*n)*(*n)-1; i>=0; i--) {
        E[i] = A[i];
        D[i] = B[i];
        eigenvalue[i] = 0;
    }
    
    // Get the length
    dsygvx(&Itype, jobs, &RANGE, &UPLO, n, E, n, D, n, &VL, &VU, &IL, &IU, 
            abstol, &m, eigtemp, C, n, dammy, &lwork, iwork, ifail, &info);
     
    // Set the workspace
    lwork = (mwSignedIndex)dammy[0];
    work = MALLOC(double, lwork);

    // Run
    dsygvx(&Itype, jobs, &RANGE, &UPLO, n, E, n, D, n, &VL, &VU, &IL, &IU, 
            abstol, &m, eigtemp, C, n, work, &lwork, iwork, ifail, &info);
    
    // Set the return values
    for (i = 0; i < *n; i++) eigenvalue[i*(*n) + i] = eigtemp[i];

    // Release the workspace
    free(eigtemp);
    free(work);
    free(iwork);
    free(E);
}
