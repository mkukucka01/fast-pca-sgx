/*****************************************************************************************************
	myEigenFunctions.h

  header information about the functions

====================================================================================================*/

#include <cmath>
#include <string>
#include <vector>
#include <sstream> //istringstream
#include <iostream> // cout
#include <fstream> // ifstream

#include "fast_pca_u.h"

using namespace std;

struct Eigenpair {
  double value; //Eigenvalue
  double *vector; //Eigenvector
  int length; //Length of eigenvector
  void normalize(oe_enclave_t* enclave, oe_result_t* result) {
    // Set eigenvalue to norm of vector and normalize eigenvector
    *result = enclave_DotProduct_vv(enclave, &value, vector, vector, length, length);
    value = sqrt(value);
    for (int i=0; i<length; i++)
      vector[i]/= value;
  }; //

  void print() {
    std::cout << value << ": \t";
    for (int i=0; i<length; i++)
      std::cout << vector[i] << "\t";
    std::cout << "\n";
  }
  // Constructor
  // Attribute value is set to 0.0 and attribute vector to an array of doubles with length n
Eigenpair(const int n) : value(0.0), length(n), vector(new double[n]) {} //Constructor

};

double** ReadData(string inputFileName, int n, int m);

double** CenterMatrix(double **A, int n, int m);

double** CovarianceMatrix(double **A, int n, int m);

Eigenpair power_method(oe_enclave_t* enclave, oe_result_t* result, double **A, double *v, int n, double tol);

void deflate(oe_enclave_t* enclave, oe_result_t* result, double **A, Eigenpair eigenpair);

void print_matrix(double **A, int n, int m);
void print_matrix(double *v, int n);