#include "enclave_u.h"
#include "sgx_urts.h"

#include <fstream>

#include "app.h"
#include "utils.h"
#include "myEigenFunctions.h"
#include "enclave.h"

using namespace std;

#define PI 3.14159265358979323846
#define TOL 1E-6

bool AE(double a, double b) {
	return (abs(a - b) < TOL); // Compare two values with respect to tolerance
}

// main function implementations
void deflate(double **A, Eigenpair eigenpair)
{
    //
    // This procedure removes eigenpair.vector from transformation matrix A in place
    //

    // store the eigenvalue before normalising
	double lambda = eigenpair.value;
	double deflator;
	eigenpair.normalize();

	for (int i = 0; i < eigenpair.length; i++) {
		for (int j = 0; j < eigenpair.length; j++) {
			// calculate value to deflate each entry in the original matrix by
			deflator = lambda * eigenpair.vector[i] * eigenpair.vector[j];
			A[i][j] = A[i][j] - deflator;
		}
	}
}

double** ReadData(string inputFileName, int n, int m)
{
	//
	//  This is a function that returns a pointer to a 56x286 matrix, which contains the reflectance spectra.
	//  The first 29 rows contain spectra from Arabica samples, the following 27 rows contain spectra from Robusta samples.
	// The code for parsing the CSV file has been adapted from Bernardo Trindade
	// https://waterprogramming.wordpress.com/2017/08/20/reading-csv-files-in-c

    double **spectra = new double* [n];
     
    for (int sample = 0; sample < n; sample++) {
	    spectra[sample] = new double [m];
    }

    vector<vector<double>> data;
    cout << inputFileName << endl;
    ifstream inputFile(inputFileName);
    int l = 0;
      
    while (inputFile) {
        string s;
        if (!getline(inputFile, s)) break;
		// cout << s << endl;
        if (l!=0) { // ignore first line, which contains wavelengths
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
  	for (int i=0; i < n; i++) {
    	for (int j=0; j < m; j++) {
      		std::cout << A[i][j] << '\t';
		}
    	std::cout << '\n';
    }
}

void print_matrix(double* v, int n)
{
	//
	// This procedure prints an nx1 vector.
	//
	for (int i = 0; i < n; i++) {
		std::cout << v[i] << '\t';
	}
}

int main(int argc, char** argv) {

    sgx_enclave_id_t eid = 0;
    sgx_launch_token_t token = {0};
    int updated, ret;
    sgx_status_t ecall_status, enclave_status;

    enclave_status = sgx_create_enclave(ENCLAVE_FILE, SGX_DEBUG_FLAG, &token, &updated, &eid, NULL);
    if (enclave_status != SGX_SUCCESS) {
        error_print("Fail to initialize enclave."); 
        return -1;
    }
    info_print("Enclave successfully initilised.");

    // ----------------------------------------------------------------------------------------------
	//
	// PART 1: Initialisation
	//
	// ----------------------------------------------------------------------------------------------

	// Defining local variables to be used:

	int n = 2, nc = 3;
	double* v, * result = NULL;

	double** A = NULL;
	double** C = NULL;
	// Allocating memory for the 1D arrays - these are the number of masses, n, long:
	v = new double[n];

	// struct comprising eigenvalue, eigenvector pair
	Eigenpair eigenpair(2), expected_eigenpair(2), deflate_eigenpair(3);

	// Allocating memory for the 2D arrays - these have dimensions of n*n:
	A = new double* [n];
	C = new double* [nc];
	for (int i = 0; i < n; i++) {
		A[i] = new double[n];
	}
	for (int i = 0; i < nc; i++) {
		C[i] = new double[nc];
	}

	if (eigenpair.value == 0.0) {
		cout << "OK: Initialization of eigenpair with eigenvalue 0.\n";
	}
	else {
		cout << "ERROR: Initialization of eigenpair with eigenvalue 0.\n";
	}

	eigenpair.vector[0] = 2.0;
	eigenpair.vector[1] = 0.0;
	eigenpair.normalize();
	if (eigenpair.value == 2.0 &&
		eigenpair.vector[0] == 1.0 &&
		eigenpair.vector[1] == 0.0) {
		cout << "OK: Normation of eigenpair.\n";
	}
	else {
		cout << "ERROR: Normation of eigenpair. ";
	}

	A[0][0] = 4.0;
	A[0][1] = 6.0;
	A[1][0] = 1.0;
	A[1][1] = 3.0;

	v[0] = 3.0;
	v[1] = 1.0;

	ecall_status = ecall_DotProduct_av(eid, &result, A, v, 2, 2*2, 2);
	if (ecall_status != SGX_SUCCESS)) {
		error_print("Failed calculation");
	}
	if (result[0] == 18.0 && result[1] == 6.0) {
		cout << "OK: DotProduct of Matrix and Vector\n";
	}
	else {
		cout << "ERROR: DotProduct of Matrix and Vector\n";
	}

	v[0] = 1.0;
	v[1] = 1.0;
	ecall_status = ecall_power_method(eid, &eignepair, A, v, 2, TOL, 2*2, 2);
	if (ecall_status != SGX_SUCCESS)) {
		error_print("Failed calculation");
	}
	if (abs(eigenpair.value - 6.0) / 6.0 < TOL) {
		cout << "OK: Power Method computing largest eigenvalue\n";
	}
	else {
		cout << "ERROR: Power Method computing largest eigenvalue\n";
	}

	if (abs(eigenpair.vector[0] - 0.94868417) / 0.94868417 < sqrt(TOL) &&
		abs(eigenpair.vector[1] - 0.31622515) / 0.31622515 < sqrt(TOL)) {
		cout << "OK: Power Method computing eigenvector of largest eigenvalue\n";
	}
	else {
		cout << "ERROR: Power Method computing eigenvector of largest eigenvalue\n";
		cout << "(" << eigenpair.vector[0] << "," << eigenpair.vector[1] << ")\n";
	}


	ifstream infile;

	infile.open("C.txt");
	for (int i = 0; i < nc; i++) {
		for (int j = 0; j < nc; j++) {
			infile >> C[i][j];
		}
	}
	cout << "Test matrix" << endl;
	print_matrix(C, 3, 3);

	double** cov;
	ecall_status = ecall_CovarianceMatrix(eid, &cov, 3, 3, 3*3);
	if (ecall_status != SGX_SUCCESS)) {
		error_print("Failed calculation");
	}
	//  cout << "CovarianceMatrix" << endl;
	// print_matrix(cov, 3, 3);
	if (AE(cov[0][0], 1.0) && AE(cov[0][1], 0.5) && AE(cov[0][2], -1.0) &&
		AE(cov[1][0], 0.5) && AE(cov[1][1], 1.0) && AE(cov[1][2], -1.0) &&
		AE(cov[2][0], -1.0) && AE(cov[2][1], -1.0) && AE(cov[2][2], 4.0 / 3)
		) {
		cout << "OK: covariance matrix computed correctly\n";
	}
	else {
		cout << "ERROR: in computation of covariance matrix\n";
		print_matrix(C, deflate_eigenpair.length, deflate_eigenpair.length);
	}


	deflate_eigenpair.value = 3.0;
	deflate_eigenpair.vector[0] = 1.0;
	deflate_eigenpair.vector[1] = 1.0;
	deflate_eigenpair.vector[2] = 0.0;

	deflate(C, deflate_eigenpair);
	if (AE(C[0][0], -0.5) && AE(C[0][1], 0.5) && AE(C[0][2], 0.0) &&
		AE(C[1][0], 0.5) && AE(C[1][1], -0.5) && AE(C[1][2], 0.0) &&
		AE(C[2][0], 0.0) && AE(C[2][1], 0.0) && AE(C[2][2], 2.0)
		) {
		cout << "OK: deflate method computes deflated matrix correctly\n";
	}
	else {
		cout << "ERROR: deflate method does not compute deflated matrix correctly\n";
		print_matrix(C, deflate_eigenpair.length, deflate_eigenpair.length);
	}

	// Test reading of CSV file
	double** spectra;

	const int N = 56;   // rows
	const int M = 286;  // columns

	spectra = ReadData("DS19hH2_dk0_FTIR_Spectra_instant_coffee.csv", N, M);
	if (AE(spectra[0][0], 21.227620) && AE(spectra[N - 1][M - 1], 1.574679)) {
		cout << "OK: Reading of CSV file" << endl;
	}
	else {
		cout << "ERROR: reading of CSV file" << spectra[0][0] << " " << spectra[N - 1][M - 1] << endl;
	}

	//print_matrix(spectra, N, M);
	ecall_status = ecall_CovarianceMatrix(eid, &cov, spectra, N, M, N*M);
	if (ecall_status != SGX_SUCCESS)) {
		error_print("Failed calculation");
	}

	// ----------------------------------------------------------------------------------------------
	// PART 2: Compute principal Components
	// ----------------------------------------------------------------------------------------------
	// Eigenvectors will be of length M

	const int n_pc = 56; // Number of principal components to compute
	double** Eigenvectors = NULL;        // 2D array where the principal components will be stored, row for each PC
	double* Eigenvalues = new double[n_pc]; // 1D array where the eigenvalues will be stored, col for each PC
	Eigenpair dominant(M);     // Eigenpair that stores the most recently computed eigenpair
	double* v_init = new double[M]; // pointer to initial eigenvector estimate

	// memory allocation
	Eigenvectors = new double* [n_pc];
	for (int i = 0; i < N; i++) {
		Eigenvectors[i] = new double[M];
	}

	// construct initial eigenvector guess
	v_init[0] = 1.0;
	for (int i = 1; i < M; i++) {
		v_init[i] = 1;
	}

	// Compute first n_pc principal components
	for (int i = 0; i < n_pc; i++) {
		ecall_status = ecall_power_method(eid, &dominant, cov, v_init, M, TOL, M*M, M); // compute dominant eigenpair (power method)
		if (ecall_status != SGX_SUCCESS)) {
		error_print("Failed calculation");
	}
		// store principal component
		Eigenvectors[i] = dominant.vector;
		Eigenvalues[i] = dominant.value;

		// deflate matrix with dominant eigenpair
		deflate(cov, dominant);
	}

	// ----------------------------------------------------------------------------------------------
	// PART 3: Exporting Data 
	// ----------------------------------------------------------------------------------------------
	
	std::ofstream output_vectors("coffeeAI_principal_components.csv");

	for (int i = 0; i < n_pc; i++) {
		for (int j = 0; j < M; j++) {
			// each row corresponds to a principal component
			output_vectors << Eigenvectors[i][j] << ",";
		}
		output_vectors << "\n";
	}
	output_vectors.close();

	std::ofstream output_values("coffeeAI_eigenvalues.csv");
	for (int i = 0; i < n_pc; i++) {
		// each row corresponds to a principal component
		output_values << Eigenvalues[i] << "\n";
	}
	output_values.close();
    
    // destroy enclave
    enclave_status = sgx_destroy_enclave(eid);
    if (enclave_status != SGX_SUCCESS) {
        error_print("Fail to destroy enclave."); 
        return -1;
    }
    info_print("Enclave successfully destroyed.");
    info_print("Program exit success.");
    return 0;
}
