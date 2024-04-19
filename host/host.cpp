


#include <fstream>
#include <openenclave/host.h>
#include "myEigenFunctions.h"
#include <sys/time.h>

// Include the untrusted helloworld header that is generated
// during the build. This file is generated by calling the
// sdk tool oeedger8r against the helloworld.edl file.
#include "fast_pca_u.h"

using namespace std;

bool check_simulate_opt(int* argc, const char* argv[])
{
    for (int i = 0; i < *argc; i++) {
        if (strcmp(argv[i], "--simulate") == 0) {
            fprintf(stdout, "Running in simulation mode\n");
            memmove(&argv[i], &argv[i + 1], (*argc - i) * sizeof(char*));
            (*argc)--;
            return true;
        }
    }
    return false;
}

#define PI 3.14159265358979323846
#define TOL 1E-6

bool AE(double a, double b) {
	return (abs(a - b) < TOL); // Compare two values with respect to tolerance
}

int main(int argc, const char* argv[]) {
 	// ----------------------------------------------------------------------------------------------
	//
	// PART 1: Initialisation and check function functionality
	//
	// ----------------------------------------------------------------------------------------------

	// Defining local variables to be used:
    oe_result_t result;
    int ret = 1;
    oe_enclave_t* enclave = NULL;

	uint64_t start; // for benchmarking
	uint64_t end;
	uint64_t elapsed;
	struct timeval tv;
	

    uint32_t flags = OE_ENCLAVE_FLAG_DEBUG;
    if (check_simulate_opt(&argc, argv)) {
        flags |= OE_ENCLAVE_FLAG_SIMULATE;
    }

    if (argc != 2) {
        fprintf(stderr, "Usage: %s enclave_image_path [ --simulate  ]\n", argv[0]);
        return 0;
    }
    // Create the enclave

	/* get start timestamp */
 	
    	gettimeofday(&tv,NULL);
    	start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;

    result = oe_create_fast_pca_enclave(
        argv[1], OE_ENCLAVE_TYPE_AUTO, flags, NULL, 0, &enclave);
    if (result != OE_OK) {
        fprintf(
            stderr,
            "oe_create_fast_pca_enclave(): result=%u (%s)\n",
            result,
            oe_result_str(result));
        return 0;
    }

	/* get elapsed time */
    	gettimeofday(&tv,NULL);
    	end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	elapsed = end - start;
	printf("@@@ Enclave creation elapsed time (usec): %lu\n", elapsed);

	int n = 2, nc = 3;
	// double* v, * result = NULL;
	double* v;

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
	/* get start timestamp */
 	
    	gettimeofday(&tv,NULL);
    	start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
	eigenpair.normalize(enclave, &result);
	/* get elapsed time */
    	gettimeofday(&tv,NULL);
    	end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	elapsed = end - start;
	printf("@@@ normalize() elapsed time (usec): %lu\n", elapsed);
	if (result != OE_OK) {
		fprintf(
			stderr,
			"error in normalize: %s\n",
			oe_result_str(result));
		return 0;
  	}
	if (eigenpair.value == 2.0 &&
		eigenpair.vector[0] == 1.0 &&
		eigenpair.vector[1] == 0.0) {
		cout << "OK: Normalization of eigenpair.\n";
	}
	else {
		cout << "ERROR: Normalization of eigenpair. ";
	}

	A[0][0] = 4.0;
	A[0][1] = 6.0;
	A[1][0] = 1.0;
	A[1][1] = 3.0;

	v[0] = 3.0;
	v[1] = 1.0;
	
	double* dot_product_rst = new double[2];
	/* get start timestamp */
 	
    	gettimeofday(&tv,NULL);
    	start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
	result = enclave_DotProduct_av(enclave, dot_product_rst, A, v, 2, 2*2, 2);
	/* get elapsed time */
    	gettimeofday(&tv,NULL);
    	end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	elapsed = end - start;
	printf("@@@ DotProduct_av() elapsed time (usec): %lu\n", elapsed);
	if (result != OE_OK) {
		fprintf(
			stderr,
			"error in DotProduct_av: %s\n",
			oe_result_str(result));
		return 0;
  	}
	if (dot_product_rst[0] == 18.0 && dot_product_rst[1] == 6.0) {
		cout << "OK: DotProduct of Matrix and Vector\n";
	}
	else {
		cout << "ERROR: DotProduct of Matrix and Vector\n";
		cout << dot_product_rst[0] << endl;
	}

	v[0] = 1.0;
	v[1] = 1.0;
	/* get start timestamp */
 	
    	gettimeofday(&tv,NULL);
    	start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
	eigenpair = power_method(enclave, &result, A, v, 2, TOL);
	/* get elapsed time */
    	gettimeofday(&tv,NULL);
    	end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	elapsed = end - start;
	printf("@@@ power_method() elapsed time (usec): %lu\n", elapsed);
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
	infile.close();

	double** cov = new double* [3];
	for (int i = 0; i < 3; i++) {
    	cov[i] = new double [3];
	}
	/* get start timestamp */
 	
    	gettimeofday(&tv,NULL);
    	start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
	result = enclave_CovarianceMatrix(enclave, cov, C, 3, 3, 3*3);
	/* get elapsed time */
    	gettimeofday(&tv,NULL);
    	end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	elapsed = end - start;
	printf("@@@ CovarianceMatrix elapsed time (usec): %lu\n", elapsed);
	if (result != OE_OK) {
		fprintf(
			stderr,
			"error in covariance matrix: %s\n",
			oe_result_str(result));
		return 0;
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

	/* get start timestamp */
 	
    	gettimeofday(&tv,NULL);
    	start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
	deflate(enclave, &result, C, deflate_eigenpair);
	/* get elapsed time */
    	gettimeofday(&tv,NULL);
    	end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	elapsed = end - start;
	printf("@@@ deflate() elapsed time (usec): %lu\n", elapsed);
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

	// ----------------------------------------------------------------------------------------------
	// PART 2: CSV Testbenches: Coffee Spectral Data
	// ----------------------------------------------------------------------------------------------
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
	/* get start timestamp */
 	
    	gettimeofday(&tv,NULL);
    	start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;

	//print_matrix(spectra, N, M);
	delete[] cov; // free memory
	cov = new double* [M];
	for (int i = 0; i < M; i++) {
    	cov[i] = new double [M];
	}
	result = enclave_CovarianceMatrix(enclave, cov, spectra, N, M, N*M);
	if (result != OE_OK) {
		fprintf(
			stderr,
			"error in covariance matrix: %s\n",
			oe_result_str(result));
		return 0;
  	}

	// ----------------------------------------------------------------------------------------------
	// Compute principal Components
	// ----------------------------------------------------------------------------------------------
	// Eigenvectors will be of length M

	const int n_pc = 40; // Number of principal components to compute
	double** Eigenvectors = new double* [n_pc];      // 2D array where the principal components will be stored, row for each PC
	double* Eigenvalues = new double[n_pc]; // 1D array where the eigenvalues will be stored, col for each PC
	Eigenpair dominant(M);     // Eigenpair that stores the most recently computed eigenpair
	double* v_init = new double[M]; // pointer to initial eigenvector estimate

	// memory allocation
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
		dominant = power_method(enclave, &result, cov, v_init, M, TOL); // compute dominant eigenpair (power method)
		// store principal component
		Eigenvectors[i] = dominant.vector;
		Eigenvalues[i] = dominant.value;

		// deflate matrix with dominant eigenpair
		deflate(enclave, &result, cov, dominant);
	}
	/* get elapsed time */
    	gettimeofday(&tv,NULL);
    	end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	elapsed = end - start;
	printf("@@@ Coffee spectral elapsed time (usec): %lu\n", elapsed);

	// ----------------------------------------------------------------------------------------------
	// Exporting Data 
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

	// // ----------------------------------------------------------------------------------------------
	// // PART 2: CSV Testbench: Milk Spectral Data
	// // ----------------------------------------------------------------------------------------------
	// // Test reading of CSV file
	
	// double** spectra;

	// const int N = 431;   // rows
	// const int M = 531;  // columns
	
	// spectra = ReadData("Milk_MIR_Traits_data_2023.csv", N, M);

	// if (AE(spectra[0][0], -0.065454953) && AE(spectra[N - 1][M - 1], 0.015999551)) {
	// 	cout << "OK: Reading of CSV file" << endl;
	// }
	// else {
	// 	cout << "ERROR: reading of CSV file " << spectra[0][0] << " " << spectra[N - 1][M - 1] << endl;
	// }
	// /* get start timestamp */
 	
    // 	gettimeofday(&tv,NULL);
    // 	start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
	
	// // print_matrix(spectra, N, M);
	// delete[] cov; // free memory
	// cov = new double* [M];
	// for (int i = 0; i < M; i++) {
    // 	cov[i] = new double [M];
	// }
	// result = enclave_CovarianceMatrix(enclave, cov, spectra, N, M, N*M);
	// if (result != OE_OK) {
	// 	fprintf(
	// 		stderr,
	// 		"error in covariance matrix: %s\n",
	// 		oe_result_str(result));
	// 	return 0;
  	// }

	// // ----------------------------------------------------------------------------------------------
	// // Compute principal Components
	// // ----------------------------------------------------------------------------------------------
	// // Eigenvectors will be of length M

	// const int n_pc = 10; // Number of principal components to compute
	// double** Eigenvectors = new double* [n_pc];        // 2D array where the principal components will be stored, row for each PC
	// double* Eigenvalues = new double[n_pc]; // 1D array where the eigenvalues will be stored, col for each PC
	// Eigenpair dominant(M);     // Eigenpair that stores the most recently computed eigenpair
	// double* v_init = new double[M]; // pointer to initial eigenvector estimate

	// // memory allocation
	// for (int i = 0; i < N; i++) {
	// 	Eigenvectors[i] = new double[M];
	// }

	// // construct initial eigenvector guess
	// v_init[0] = 1.0;
	// for (int i = 1; i < M; i++) {
	// 	v_init[i] = 1;
	// }

	// // Compute first n_pc principal components
	// for (int i = 0; i < n_pc; i++) {
	// 	dominant = power_method(enclave, &result, cov, v_init, M, TOL); // compute dominant eigenpair (power method)
	// 	if (result != OE_OK) {
	// 		fprintf(
	// 			stderr,
	// 			"error in power method: %s\n",
	// 			oe_result_str(result));
	// 		return 0;
  	// 	}
	// 	// store principal component
	// 	Eigenvectors[i] = dominant.vector;
	// 	Eigenvalues[i] = dominant.value;

	// 	// deflate matrix with dominant eigenpair
	// 	deflate(enclave, &result, cov, dominant);
	// 	if (result != OE_OK) {
	// 		fprintf(
	// 			stderr,
	// 			"error in deflate matrix: %s\n",
	// 			oe_result_str(result));
	// 		return 0;
  	// 	}
	// }
	// /* get elapsed time */
    // 	gettimeofday(&tv,NULL);
    // 	end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    // 	elapsed = end - start;
	// printf("@@@ Milk spectral testbench elapsed time (usec): %lu\n", elapsed);
	// // ----------------------------------------------------------------------------------------------
	// // Exporting Data 
	// // ----------------------------------------------------------------------------------------------
	
	// ofstream output_vectors("milk_principal_components.csv");

	// for (int i = 0; i < n_pc; i++) {
	// 	for (int j = 0; j < M; j++) {
	// 		// each row corresponds to a principal component
	// 		// cout << Eigenvectors[i][j] << ",";
	// 		output_vectors << Eigenvectors[i][j] << ",";
	// 	}
	// 	output_vectors << "\n";
	// 	// cout << "\n";
	// }
	// output_vectors.close();

	// std::ofstream output_values("milk_eigenvalues.csv");

	// for (int i = 0; i < n_pc; i++) {
	// 	// each row corresponds to a principal component
	// 	output_values << Eigenvalues[i] << "\n";
	// }
	// output_values.close();
	
	
	// // ----------------------------------------------------------------------------------------------
	// // PART 2: CSV Testbench: Grape Spectral Data
	// // ----------------------------------------------------------------------------------------------
	// // Test reading of CSV file
	
	// double** spectra;

	// const int N = 274;   // rows
	// const int M = 205;  // columns
	
	// spectra = ReadData("grape_maturity_dataset.csv", N, M);

	// if (AE(spectra[0][0], 144.738) && AE(spectra[N - 1][M - 1], 0.36499172)) {
	// 	cout << "OK: Reading of CSV file" << endl;
	// }
	// else {
	// 	cout << "ERROR: reading of CSV file " << spectra[0][0] << " " << spectra[N - 1][M - 1] << endl;
	// }
	// // print_matrix(spectra, N, M);
	// /* get start timestamp */
 	
    // 	gettimeofday(&tv,NULL);
    // 	start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
	
	// // print_matrix(spectra, N, M);
	// delete[] cov; // free memory
	// cov = new double* [M];
	// for (int i = 0; i < M; i++) {
    // 	cov[i] = new double [M];
	// }
	// result = enclave_CovarianceMatrix(enclave, cov, spectra, N, M, N*M);
	// if (result != OE_OK) {
	// 	fprintf(
	// 		stderr,
	// 		"error in covariance matrix: %s\n",
	// 		oe_result_str(result));
	// 	return 0;
  	// }

	// // ----------------------------------------------------------------------------------------------
	// // Compute principal Components
	// // ----------------------------------------------------------------------------------------------
	// // Eigenvectors will be of length M

	// const int n_pc = 10; // Number of principal components to compute
	// double** Eigenvectors = new double* [n_pc];        // 2D array where the principal components will be stored, row for each PC
	// double* Eigenvalues = new double[n_pc]; // 1D array where the eigenvalues will be stored, col for each PC
	// Eigenpair dominant(M);     // Eigenpair that stores the most recently computed eigenpair
	// double* v_init = new double[M]; // pointer to initial eigenvector estimate

	// // memory allocation
	// for (int i = 0; i < N; i++) {
	// 	Eigenvectors[i] = new double[M];
	// }

	// // construct initial eigenvector guess
	// v_init[0] = 1.0;
	// for (int i = 1; i < M; i++) {
	// 	v_init[i] = 1;
	// }

	// // Compute first n_pc principal components
	// for (int i = 0; i < n_pc; i++) {
	// 	dominant = power_method(enclave, &result, cov, v_init, M, TOL); // compute dominant eigenpair (power method)
	// 	if (result != OE_OK) {
	// 		fprintf(
	// 			stderr,
	// 			"error in power method: %s\n",
	// 			oe_result_str(result));
	// 		return 0;
  	// 	}
	// 	// store principal component
	// 	Eigenvectors[i] = dominant.vector;
	// 	Eigenvalues[i] = dominant.value;

	// 	// deflate matrix with dominant eigenpair
	// 	deflate(enclave, &result, cov, dominant);
	// 	if (result != OE_OK) {
	// 		fprintf(
	// 			stderr,
	// 			"error in deflate matrix: %s\n",
	// 			oe_result_str(result));
	// 		return 0;
  	// 	}
	// }
	// /* get elapsed time */
    // 	gettimeofday(&tv,NULL);
    // 	end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    // 	elapsed = end - start;
	// printf("@@@ Grape spectral testbench elapsed time (usec): %lu\n", elapsed);
	// // ----------------------------------------------------------------------------------------------
	// // Exporting Data 
	// // ----------------------------------------------------------------------------------------------
	
	// std::ofstream output_vectors("grape_principal_components.csv");

	// for (int i = 0; i < n_pc; i++) {
	// 	for (int j = 0; j < M; j++) {
	// 		// each row corresponds to a principal component
	// 		// cout << Eigenvectors[i][j] << ",";
	// 		output_vectors << Eigenvectors[i][j] << ",";
	// 	}
	// 	output_vectors << "\n";
	// 	// cout << "\n";
	// }
	// output_vectors.close();

	// std::ofstream output_values("grape_eigenvalues.csv");

	// for (int i = 0; i < n_pc; i++) {
	// 	// each row corresponds to a principal component
	// 	output_values << Eigenvalues[i] << "\n";
	// }
	// output_values.close();

	
	delete[] Eigenvalues;
	delete[] Eigenvectors;
	// Clean up the enclave if we created one
	if (enclave) {
		oe_terminate_enclave(enclave);
	}
	return 0;
} // end main
