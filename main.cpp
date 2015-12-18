#include <iostream>
#include <stdlib.h>

#include <f2c.h>
#include <clapack.h>
#include <lapack.h>
#include <blaswrap.h>

using namespace std;


double callFrobenius()
{
  integer rows = 3;
  integer cols = 3;
  char norm = 'F';
  doublereal* matrix = (doublereal*)malloc(sizeof(doublereal)*rows*cols);
  
  matrix[0] = 1;
  matrix[1] = 2;
  matrix[2] = 1;
  
  matrix[3] = -1;
  matrix[4] = 2;
  matrix[5] = -3;

  matrix[6] = 0;
  matrix[7] = 1;
  matrix[8] = -2;

  doublereal* result = NULL;
  // doublereal dlange_(char *norm, integer *m, integer *n, doublereal *a, integer 
  // *lda, doublereal *work);
  
  return  dlange_(&norm, &rows, &cols, matrix, &rows, result);

  

}

int main ()
{
  cout << "Norm is" << callFrobenius() << endl;
  return 0;
}