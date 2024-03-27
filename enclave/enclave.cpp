#include "enclave_t.h"
#include "string.h"
#include <math.h>
#include <vector>
#include <sstream> //istringstream
#include <iostream> // cout
#include <fstream> // ifstream

#include "enclave.h"
#include "myEigenFunctions.h"


double ecall_DotProduct_aa(double **A, double **B, int n, int m,  size_t len)
{
	//
	//	This is a function that takes two matrices A and B of identical dimensions (n*m) and 
	//  calculates and returns their dot product.
	//
	double dot = 0.0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			dot += A[i][j]*B[i][j];
		}
	}

	return dot;
}

double ecall_DotProduct_vv(double *A, double *B, int n, size_t len) 
{
	//
	//	This is a function that takes two vectors A and B of identical length (n) and 
	//  calculates and returns their dot product.
	//

	double dot = 0.0;

	for (int i = 0; i < n; i++) {
		dot += A[i]*B[i];
	}

	return dot;
}

double* ecall_DotProduct_av(double **A, double *v, int n, size_t len1, size_t len2)
{
	//
	//  This is a function that takes a nxn-matrix A and an n-dimensional vector v stores
	//  the product A.v at the original location of v
	//

    double *result = (double *) malloc(sizeof(double)*n); // pointer to result vector
  
    for (int i = 0; i < n; i++) {
      	result[i] = 0.0; // initialize ith element of result v
      	for (int j = 0; j < n; j++) {
        	result[i] += A[i][j] * v[j]; 
      	}	
    }

    return result;
}

double** ecall_CenterMatrix(double **A, int n, int m, size_t len)
{
	//
	//  This is a function that takes a nxm-matrix A and subtracts the emperical mean of each column.
	//  The function returns the point to the centered matrix
	//
      
	// initialise new vars
    double **result =  (double **) malloc(sizeof(double*)*n);
	double mean;
      
	// allocate memory
	for (int row = 0; row < n; row++) {
		result[row] =  (double) malloc(sizeof(double)*n);
	}

	// loop over entire matrix
	for (int col = 0; col < m; col++) {
		mean = 0;
		// get sum of column
		for (int row = 0; row < n; row++) {
			mean += A[row][col];
		}
		// calculate mean
		mean = mean / n;
		// subtract mean from result
		for (int row = 0; row < n; row++) {
			result[row][col] = A[row][col] - mean;
		}
	}
	  
    return result;
}

double** ecall_CovarianceMatrix(double **A, int n, int m, size_t len)
{
	//
	//  This is a function that takes a nxm-matrix A and computes its mxm covariance matrix.
	//
      
    double **cov =  (double *) malloc(sizeof(double*)*m);
    double **cA = ecall_CenterMatrix(A, n, m);
      
	for(int i = 0; i < m; i++)
	cov[i] =  (double) malloc(sizeof(double)*n);

	// loop over covariance matrix
	for (int row = 0; row < m; row++) {
		for (int col = 0; col < m; col++) {
			// initialise matrix entry
			cov[row][col] = 0;
			for (int i = 0; i < n; i++) {
				// calculate the entry from multiplying cA with its transpose
				cov[row][col] += cA[i][row] * cA[i][col];
			}
			cov[row][col] /= n - 1.0; // divide by n-1
		}
	}
      
    return cov;
}

Eigenpair ecall_power_method(double **A, double *v, int n, double tol, size_t len1, size_t len2)
{
    //
    // This function computes the largest eigenvalue and the corresponding eigenvector
    //
    // v: initial eigenvector iterate

    Eigenpair eigenpair(n);
    double lambda;

    for (int i = 0; i < n; i++) {
	    eigenpair.vector[i] = v[i];
    }
	eigenpair.normalize();

    do {
	    lambda = eigenpair.value;
	    eigenpair.vector = ecall_DotProduct_av(A, eigenpair.vector, eigenpair.length);
	    eigenpair.normalize();
	  
  	} while (abs(eigenpair.value - lambda)/abs(eigenpair.value) > tol);
  
  	return eigenpair;
}
