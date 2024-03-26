#ifndef MYEIGENFUNCTIONS_H_
#define MYEIGENFUNCTIONS_H_

#include <math.h>
#include <string.h>
#include <vector>
#include <sstream> //istringstream
#include <iostream> // cout
#include <fstream> // ifstream

using namespace std;

struct Eigenpair {
  double value; //Eigenvalue
  double *vector; //Eigenvector
  int length; //Length of eigenvector
  void normalize() {
    // Set eigenvalue to norm of vector and normalize eigenvector
    value = sqrt(ecall_DotProduct_vv(vector, vector, length));
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

#endif // MYEIGENFUNCTIONS_H_
