/*****************************************************************************************************

 Source information for Eigenvalue decomposition functions

====================================================================================================*/

#include <string>

#include "myEigenFunctions.h"
#include "fast_pca_u.h"

using namespace std;

double** ReadData(string inputFileName, int n, int m)
    {
          //
          //  This is a function that returns a pointer to a 56x286 matrix, which contains the reflectance spectra.
          //  The first 29 rows contain spectra from Arabica samples, the following 27 rows contain spectra from Robusta samples.
          // The code for parsing the CSV file has been adapted from Bernardo Trindade
          // https://waterprogramming.wordpress.com/2017/08/20/reading-csv-files-in-c

      double **spectra = new double* [n];

     
      for(int sample = 0; sample < n; sample++)
	spectra[sample] = new double [m];

    vector<vector<double> > data;
    cout << inputFileName << endl;
    ifstream inputFile(inputFileName);
    int l = 0;
      
    while (inputFile) {
        string s;
        if (!getline(inputFile, s)) break;
	// cout << s << endl;
        if (l!=0) // ignore first line, which contains wavelengths
	  {
            istringstream ss(s);
            vector<double> record;

	    int column = 0;
            while (ss) {
	      // cout << "Row " << l << " " << ss << endl;
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
		  // cout << "Row " << l << " Column " << column << " line " << line << endl; 
                    spectra[l-1][column] = stod(line);
		    column++;
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
		}

        }
	l++;
    }

    return spectra;
  }

void print_matrix(double **A, int n, int m)
{
  //
  // This procedure prints an nxm matrix.
  //
  for (int i=0; i < n; i++)
    {
    for (int j=0; j < m; j++)
      std::cout << A[i][j] << '\t';
    std::cout << '\n';
    }
}

void print_matrix(double* v, int n)
{
	//
	// This procedure prints an nx1 vector.
	//
	for (int i = 0; i < n; i++)
	{
		std::cout << v[i] << '\t';
	}
}

Eigenpair power_method(oe_enclave_t* enclave, oe_result_t* result, double **A, double *v, int n, double tol)
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
	eigenpair.normalize(enclave, result);
    do {
	    lambda = eigenpair.value;
	    *result = enclave_DotProduct_av(enclave, &eigenpair.vector, A, eigenpair.vector, eigenpair.length, eigenpair.length*eigenpair.length, eigenpair.length);
	    eigenpair.normalize(enclave, result);
	  
  } while (abs(eigenpair.value - lambda)/abs(eigenpair.value) > tol);
  
      
  return eigenpair;
}

void deflate(oe_enclave_t* enclave, oe_result_t* result, double **A, Eigenpair eigenpair)
{
    //
    // This procedure removes eigenpair.vector from transformation matrix A in place
    //

    // store the eigenvalue before normalising
	double lambda = eigenpair.value;
	double deflator;
	eigenpair.normalize(enclave, result);

	for (int i = 0; i < eigenpair.length; i++) {
		for (int j = 0; j < eigenpair.length; j++) {
			// calculate value to deflate each entry in the original matrix by
			deflator = lambda * eigenpair.vector[i] * eigenpair.vector[j];
			A[i][j] = A[i][j] - deflator;
		}
	}
}

double** CenterMatrix(double **A, int n, int m)
    {
          //
          //  This is a function that takes a nxm-matrix A and subtracts the emperical mean of each column.
          //  The function returns the point to the centered matrix
          //
      
	// initialise new vars
      double **result = new double* [n];
	  double mean;
      
	  // allocate memory
      for(int row = 0; row < n; row++)
	result[row] = new double [m];

	  // loop over entire matrix
	  for (int col = 0; col < m; col++) {
		  mean = 0;
		  // get sum of column
		  for (int row = 0; row < n; row++)
		mean += A[row][col];
		  // calculate mean
		  mean = mean / n;
		  // subtract mean from result
		  for (int row = 0; row < n; row++)
		result[row][col] = A[row][col] - mean;
	  }
    return result;
  }

double** CovarianceMatrix(double **A, int n, int m)
    {
          //
          //  This is a function that takes a nxm-matrix A and computes its mxm covariance matrix.
          //
      
      double **cov = new double* [m];
      double **cA = CenterMatrix(A, n, m);
      
      for(int i = 0; i < m; i++)
	cov[i] = new double [m];

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
