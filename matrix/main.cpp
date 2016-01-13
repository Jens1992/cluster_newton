#include <iostream>
#include <stdlib.h>

#include <f2c.h>

#include <copasi/utilities/CMatrix.h>

extern "C"
{

#include <clapack.h>
#include <lapack.h>
#include <blaswrap.h>

}
using namespace std;


double callFrobenius(CMatrix<double>& matrix)
{
  integer rows = 3;
  integer cols = 3;
  char norm = 'F';

 
  doublereal* result = NULL;
  // doublereal dlange_(char *norm, integer *m, integer *n, doublereal *a, integer 
  // *lda, doublereal *work);

  return  dlange_(&norm, &rows, &cols, matrix.array(), &rows, result);


}

int main()
{

  // create matrix
  CMatrix<double> matrix(3, 3);
  matrix(0, 0) = 1;
  matrix(0, 1) = 2;
  matrix(0, 2) = 1;

  matrix(1, 0) = -1;
  matrix(1, 1) = 2;
  matrix(1, 2) = -3;

  matrix(2, 0) = 0;
  matrix(2, 1) = 1;
  matrix(2, 2) = -2;


  cout << "Norm is: " << callFrobenius(matrix) << endl;
  return 0;
}