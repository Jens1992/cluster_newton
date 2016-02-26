//#include "stdafx.h"

#include <iostream>
#include <stdlib.h>

#include <f2c.h>

#include <vector>
#include <string>
#include <sstream>
#include <stdio.h>
#include <time.h>

#include<math.h>
#ifdef WIN32
#include<conio.h>
#endif
#include <cstdio>


extern "C" 
{
	#include <clapack.h>
	#include <lapack.h>
	#include <blaswrap.h>
}

using namespace std;

const int number_of_random_points = 10;
const int number_of_parameters_from_copasi = 5;
const int K1 = 10; // number of iterations is K1 + 1 in the first Stage
const int K2 = K1 + 20; // starts at K1+1, number of iterations is K2+1 in the second Stage 
/// K2=K1+2 is just temporarily that K2 is bigger than K1 because the iteration of K2 starts at K1+1
const int number_of_equations = 6;

const double eta = 0.1;

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
   doublereal dlange_(char *norm, integer *m, integer *n, doublereal *a, integer 
   *lda, doublereal *work);
  
  return  dlange_(&norm, &rows, &cols, matrix, &rows, result);
  
}

vector<vector<double > > func_vector_to_diag_matrix(vector<vector<double > > vector_input, int nRows, int nCols, string row_or_col_vec)
{
	/**************************************************************************************************************
	* this function creates from a vector given as a matrix with either only 1 column or 1 row a diagonal matrix  *
	* with the vector's elements as the diagonal elements of the new matrix.									  *
	***************************************************************************************************************/

	int nRowCol;
	if (row_or_col_vec == "row")
	{
		nRowCol = nCols;
	}
	else if (row_or_col_vec == "col")
	{
		nRowCol = nRows;
	}

	vector<vector<double > > matrix(nRowCol, vector<double>(nRowCol));
	
	for (int i = 0; i < nRowCol; i++)
	{
		for (int j = 0; j < nRowCol; j++)
		{
			if (i == j)
			{
				if (row_or_col_vec == "row")
				{
					matrix[i][j] = vector_input[0][i];
				}
				else if (row_or_col_vec == "col")
				{
					matrix[i][j] = vector_input[i][0];
				}
			}
			else
			{
				matrix[i][j] = 0;
			}
		}
	}

	return matrix;

	/**************************************************************************************************************
	* END of function creating a diagonal matrix from a vector given as a matrix.								  *
	***************************************************************************************************************/
}

void matrix_output(vector<vector<double > > matrix, int rows, int cols, string text_before, string text_after)
{
	cout << text_before + "\n";

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			printf("%.20f", matrix[i][j]);
			cout << "  ";
		}
		printf("\n");
	}

	cout << text_after + "\n\n";
}

vector<vector<double > > func_double_star_to_vector_matrix(double* double_array, int nRows, int nCols)

{
	/****************************************************************************************************
	* This function converts a double* array to a vector<vector<double > > matrix.						*
	*****************************************************************************************************/

	vector<vector<double > > matrix(nRows, vector<double>(nCols));
	int array_length;
	int memindex;

	array_length = nRows * nCols;
	memindex = 0;

	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			matrix[i][j] = double_array[memindex];
			memindex++;
		}
	}

	return matrix;

	/****************************************************************************************************
	* END of function converting double* array to vector<vector<double>>								*
	*****************************************************************************************************/

}

vector<vector<double > > func_double_star_star_to_vector_matrix(double** double_array, int nRows, int nCols)
{
	/****************************************************************************************************
	* This function converts a double** array to a vector<vector<double > > matrix.						*
	*****************************************************************************************************/

	vector<vector<double > > matrix(nRows, vector<double>(nCols));
	int array_length;
	int memindex;

	array_length = nRows * nCols;

	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			matrix[i][j] = double_array[i][j];
		}
	}

	return matrix;

	/****************************************************************************************************
	* END of function converting double** array to vector<vector<double>>								*
	*****************************************************************************************************/
}

double* func_vector_matrix_to_double_star(vector<vector<double > > matrix, int nRows, int nCols)
{
	/****************************************************************************************************
	* This function converts a vector<vector<double > > matrix into a double* array i.e. to use it in		*
	* the function inverse to create inverse matrices.													*
	*****************************************************************************************************/

	double* double_array = new double[nRows*nCols];
	int memindex = 0;

	for (int i = 0; i < nRows; ++i)
	{
		for (int j = 0; j < nCols; j++)
		{
			double_array[memindex] = matrix[i][j];
			++memindex;
		}
	}

	return double_array;

	/****************************************************************************************************
	* END of function converting vector matrix to double* array											*
	*****************************************************************************************************/
}

void inverse(double* A, integer N)
{
	/****************************************************************************************
	* This function calculates the inverse of a given matrix A								*
	*****************************************************************************************/

	integer *IPIV = new integer[N + 1];
	integer LWORK = N*N;
	double *WORK = new double[LWORK];
	integer INFO;

	// LU decomoposition of a general matrix
	//void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

	// generate inverse of a matrix given its LU decomposition
	//void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

	dgetrf_(&N, &N, A, &N, IPIV, &INFO);
	dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);

	delete IPIV;
	delete WORK;

	/****************************************************************************************
	* END of function calculating a given matrix A											*
	*****************************************************************************************/
}



void LGS(complex *matrix)
{

	/****************************************************************************************
	* This function should have calculate an LGS but doesn'T work yet						*
	*****************************************************************************************/

	integer rows = 3;
	integer cols = 3;
	integer *IPIV = new integer[rows];
	const int x = 9;

	integer cols2 = 1;

	char norm = 'F';

	//complex matrix[x];
	complex b[3];
	integer info;

  matrix = new complex[9];//();
	

//	doublereal* matrix = (doublereal*)malloc(sizeof(doublereal)*rows*cols);

//	doublereal* b = (doublereal*)malloc(sizeof(doublereal)*rows*cols2);
	
	matrix[0] = { 3,0 };

	matrix[1] = { 2, 0 };
	matrix[2] = { -1, 0 };

	matrix[3] = { 2, 0 };
	matrix[4] = { -2, 0 };
	matrix[5] = { 4, 0 };

	matrix[6] = { -1, 0 };
	matrix[7] = { 0.5, 0 };
	matrix[8] = { -1, 0 };
	
	b[0] = { 1, 0 };
	b[1] = { -2, 0 };
	b[2] = { 0, 0 };
	

	doublereal* result = NULL;
	//doublereal CGESV(char *norm, integer *m, integer *n, doublereal *a, integer
	//	*lda, doublereal *work);

	 //cgesv_(&rows,&cols2,matrix,&rows,IPIV,b,&rows,&info);

	 
	 /***********************************************************************************************
	 * END of function which should calculate a linear system of linear equations.					*
	 ************************************************************************************************/
}

// functions used in general

vector<vector<double > > func_matrix_addition(vector< vector<double > > matrixA, vector< vector<double > > matrixB, int nZeilen, int nSpalten)
{
	vector< vector<double > > matrixC(nZeilen, vector<double>(nSpalten));

	for (int i = 0; i < nZeilen; i++)
	{
		for (int j = 0; j < nSpalten; j++)
		{
			matrixC[i][j] = matrixA[i][j] + matrixB[i][j];
		}
	}

	return matrixC;
}

vector<vector<double > > func_matrix_difference(vector< vector<double > > matrixA, vector< vector<double > > matrixB, int nZeilen, int nSpalten)
{
	vector< vector<double > > matrixC(nZeilen, vector<double>(nSpalten));

	for (int i = 0; i < nZeilen; i++)
	{
		for (int j = 0; j < nSpalten; j++)
		{
			matrixC[i][j] = matrixA[i][j] - matrixB[i][j];
		}
	}

	return matrixC;
}


/************************************************************************************************************
* calculates the 2 norm of a vector																			*
*************************************************************************************************************/

double func_2_norm(vector<vector<double > > vector_matrix, int rows, int cols, string row_or_col)
{
	double sum_2_norm = 0;

	if (row_or_col == "row")
	{
		for (int i = 0; i < cols; i++)
		{
			sum_2_norm = sum_2_norm + vector_matrix[0][i] * vector_matrix[0][i];
		}
		sum_2_norm = sqrt(sum_2_norm);
	}
	else if (row_or_col == "col")
	{
		for (int i = 0; i < rows; i++)
		{
			sum_2_norm = sum_2_norm + vector_matrix[i][0] * vector_matrix[i][0];
		}
		sum_2_norm = sqrt(sum_2_norm);
	}

	return sum_2_norm;

}

/************************************************************************************************************
* calculates the 2 norm of a vector																			*
*************************************************************************************************************/

vector<double> func_get_diag(vector<vector<double > > matrix, int nZeilen, int nSpalten)
{
	vector<double> diagonal_elements;

	if (nZeilen != nSpalten)
	{
		cout << "Fehler, Zeilen und Spalten müssen bei der Ermittlung der Diagonalelemente gleich groß sein!";
	}
	else
	{
		for (int i = 0; i < nZeilen; i++)
		{
			diagonal_elements.push_back(matrix[i][i]);
		}
	}

	return diagonal_elements;
}

vector<vector<double > > func_take_vector_of_matrix(vector<vector<double > > matrix, int rows, int cols, int vector_row_or_col_nr, string row_or_col)
{
	vector<vector<double > > row_vector(1, vector<double>(cols));
	vector<vector<double > > col_vector(rows, vector<double>(1));

	if (row_or_col == "row")
	{
		for (int j = 0; j < cols; j++)
		{
			row_vector[0][j] = matrix[vector_row_or_col_nr][j];
		}

		return row_vector;
	}
	else if (row_or_col == "col")
	{
		for (int i = 0; i < rows; i++)
		{
			col_vector[i][0] = matrix[i][vector_row_or_col_nr];
		}
		return col_vector;
	}
	else
	{
		cout << "Es muss im letzten Parameter der Funktion func_take_vector_of_matrix festgelegt entweder \"row\", wenn ein Zeilenvektor oder \"col\", wenn ein Spaltenvektor erstellt werden soll.";
	}	
}
//func_add_vector_to_matrix(Xinverse_mult_X_mult_Y_transpose, 1, m + 1, i, "row", A);
vector<vector<double > > func_add_vector_to_matrix(vector<vector<double > > add_vector, int rows, int cols, int vector_row_or_col_nr, string row_or_col, vector<vector<double > > matrix)
{
	if (row_or_col == "row")
	{
		for (int i = 0; i < cols; i++)
		{
			matrix[vector_row_or_col_nr][i] = add_vector[rows-1][i];
		}
	}
	else if (row_or_col == "col")
	{
		for (int i = 0; i < rows; i++)
		{
			matrix[i][vector_row_or_col_nr] = add_vector[i][cols-1];
		}
	}
	else
	{
		cout << "Es muss im vorletzten Parameter der Funktion func_add_vector_to_matrix festgelegt entweder \"row\", wenn ein Zeilenvektor oder \"col\", wenn ein Spaltenvektor zur Matrix hinzugefügt werden soll.";
	}
	return matrix;
}

vector<vector<double > > func_create_x_tilde(vector<vector<double > > X, int m, int l)
{
	vector<vector<double > > X_tilde(m+1,vector<double>(l));

	for (int i = 0; i <= m; i++)
	{
		for (int j = 0; j < l; j++)
		{

			if (i < m)
			{
				X_tilde[i][j] = X[i][j];
			}
			else
			{
				X_tilde[i][j] = 1;				
			}
		}
	}

	return X_tilde;
}

vector<vector<double > > func_matrix_transpose(vector<vector<double > > matrix, int nZeilen, int nSpalten)
{
	vector<vector<double > > transposed_matrix;

	for (int col_nr = 0; col_nr < nSpalten; col_nr++)
	{
		vector<double> row_elements;
		for (int row_nr = 0; row_nr < nZeilen; row_nr++)
		{
			row_elements.push_back(matrix[row_nr][col_nr]);
		}
		transposed_matrix.push_back(row_elements);
	}

	return transposed_matrix;
}

vector<vector<double > > func_matrix_multiplication(vector<vector<double > > matrixA
	, vector<vector<double > > matrixB
	, int nZeilen_matrixA
	, int nSpalten_matrixA
	, int nZeilen_matrixB
	, int nSpalten_matrixB
	)
{
	vector<vector<double > > matrixC(nZeilen_matrixA, vector<double>(nSpalten_matrixB));
	double new_element = 0;

	if (nSpalten_matrixA != nZeilen_matrixB)
	{
		cout << "\n";
		cout << "Bei der Matrixmultiplikation stimmt die Spaltenzahl der ersten Matrixn nicht mit der Zeilenzahl der zweiten Matrix überein!";
		cout << "\n";
	}
	else
	{
		// multiply each row of matrix A: 
		for (int j = 0; j < nZeilen_matrixA; j++)
		{
			// ... with each col of matrix B:
			for (int k = 0; k < nSpalten_matrixB; k++)
			{
				new_element = 0;
				// ... therefor use each element of row j: 
				for (int i = 0; i < nSpalten_matrixA; i++)
				{
					// ... and multiply this element ji of matrixA with the element ik of 
					// matrixB and add all of them
					// to the new element jk of matrixC 
					new_element += matrixA[j][i] * matrixB[i][k];
				}
				matrixC[j][k] = new_element;
			}
		}
	}
	return matrixC;
}

vector<vector<double > > func_matrix_scalar_multiplication(double scalar
	, vector<vector<double > > matrix
	, int nZeilen
	, int nSpalten
	)
{
	for (int i = 0; i < nZeilen; i++)
	{
		for (int j = 0; j < nSpalten; j++)
		{
			matrix[i][j] = matrix[i][j] * scalar;
		}
	}

	return matrix;
}

/*

* C++ Program to Find the Determinant of a Given Matrix

*/

double det(int n, vector<vector<double > > matrix)
{
	int c, subi, i, j, subj;
	vector<vector<double > > submatrix(n, vector<double>(n));
	double determinante = 0;
	if (n == 2)
	{
		return((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
	}
	else
	{
		for (c = 0; c < n; c++)
		{
			subi = 0;
			for (i = 1; i < n; i++)
			{
				subj = 0;
				for (j = 0; j < n; j++)
				{
					if (j == c)
					{
						continue;
					}
					submatrix[subi][subj] = matrix[i][j];
					subj++;
				}
				subi++;
			}
			determinante = determinante + (pow(-1, c) * matrix[0][c] * det(n - 1, submatrix));
		}
	}
	return determinante;
}

/// Stage 1

// 1-1 -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/********************************************************************************************
* This part of the program is there to produce an cluster of random points based on (2.1)   *
*********************************************************************************************/

/********************************************************************************************
* This function checks whether the function of Stage 1  1: 1-1: (2.1) is fulfilled.			*
* for a given random value xij(0), here for the value value_in_range						*
*********************************************************************************************/

bool func_check_if_value_fulfills_random_condition(double value_in_range
	, double xDach
	, double max_range
	, double min_range
	)

{
	double left_side_of_condition;
	double right_side_of_condition;
	bool bool_condition_fulfilled = false;

	left_side_of_condition = abs((value_in_range - xDach) / xDach);
	right_side_of_condition = max_range - min_range;

	if (left_side_of_condition < right_side_of_condition)
	{
		bool_condition_fulfilled = true;
	}
	else
	{
		bool_condition_fulfilled = false;
	}

	return bool_condition_fulfilled;
}

/****************************************************************************************************
* This functions creates a random value between min_range and max_range which has to fulfill		*
* the condition (2.1)																				*
*****************************************************************************************************/

double fRand(double max_range, double min_range, double xDach)
{
	double f;
	double value_in_range = 0;
	int memindex = 0;

	bool bool_random_erfuellt_bedingung = false;

	while (bool_random_erfuellt_bedingung == false && memindex <= 10)
	{
		f = (double)rand() / RAND_MAX;
		value_in_range = min_range + f * (max_range - min_range);
		bool_random_erfuellt_bedingung = func_check_if_value_fulfills_random_condition(value_in_range
			, xDach
			, max_range
			, min_range
			);
		++memindex;
	}
	return value_in_range;
}

/*******************************************************************************************************
* This function writes the randomly created values of fRand in an Array. This array has a specific     *
* length given by the constant number_of_random_points.												   *
********************************************************************************************************/

double** func_get_random_points(double* min_range
	, double* max_range
	, double* xDach
	)
{
	//double* random_points = new double[number_of_random_points][number_of_parameters_from_copasi];
	double** random_points = 0;
	double random_number;
	random_points = new double*[number_of_parameters_from_copasi]; // rows of array

	for (int i = 0; i < number_of_parameters_from_copasi; i++)
	{
		random_points[i] = new double[number_of_random_points]; // cols of array
		//random_points[i] = { 0,0,0,0,0 };
		for (int j = 0; j < number_of_random_points; j++)
		{
			random_number = fRand(max_range[i], min_range[i], xDach[i]);
			random_points[i][j] = random_number;
		}
	}

	return random_points;
}

/********************************************************************************************************
* END of program part which produces a cluster of random points (two dimensional array)					*
*********************************************************************************************************/

// END of 1-1 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// 1-2 -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/********************************************************************************************************
* program part, which produces a random vector near y* with condition at								*
* maximum abs(yij -yi* / yi*) < 0.1																		*
*********************************************************************************************************/


vector<vector<double > > func_get_random_vector_near_ystar(vector<double> ystar)
{
	vector<vector<double > > random_points(number_of_equations, vector<double>(number_of_random_points));
	double random_number;

	vector<double> min_range(number_of_equations);
	vector<double> max_range(number_of_equations);

	/****************************************************************************************
	* set min_Range and max_range considering the condition from 1-2: (2.2)					*
	*****************************************************************************************/

	for (int i = 0; i < number_of_equations; i++)
	{
		if (ystar[i] != 0)
		{
			min_range[i] = ystar[i] - eta * ystar[i] + 0.0000000001;
			//cout << "\n min range: ";
			//cout << min_range[i];
			//cout << "\n";
			max_range[i] = eta * ystar[i] + ystar[i] - 0.0000000001;
			//cout << "\n max range: ";
			//cout << max_range[i];
			//cout << "\n";
		}
		else
		{
			min_range[i] = 0;
			max_range[i] = 0.1;
		}

	}

	/****************************************************************************************
	* END of setting values of max_range and min_range considering condition (2.2)			*
	*****************************************************************************************/

	/****************************************************************************************
	* Create random values considering min_Range and max_range								*
	*****************************************************************************************/
	
	for (int i = 0; i < number_of_equations; i++)
	{
		for (int j = 0; j < number_of_random_points; j++)
		{			
			random_number = fRand(max_range[i], min_range[i], ystar[i]);
			random_points[i][j] = random_number;
		}
	}

	/****************************************************************************************
	* END of creating random values															*
	*****************************************************************************************/
	return random_points;
}

/********************************************************************************************************
* END of program part which produces an maximum vector near y* with a specific condition (2.2)			*
*********************************************************************************************************/

// END of 1-2 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// 2-1 -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/********************************************************************************************************
* Calculate the forward problem y = f(x)																*
*********************************************************************************************************/

vector<vector<double > > calculate_forward_problem(vector<vector<double > > x, int xRows, int xCols, int yRows)
{
	// yRows is the number of equations (n)
	// xRows is the number of parameters (m)
	// xCols is the number of randomly created points in 1-1 (l)
	// yCols is the same as xCols -> (l)

	// y: n x l

	int yCols = xCols;
	vector<vector<double > > y(yRows,vector<double>(yCols));
	
	// example of values of y created by forward problem (the real forward problem algorithm has to be called here!)
	int memindex = 0;
	for (int i = 0; i < yRows; i++)
	{
		for (int j = 0; j < yCols; j++)
		{
			y[i][j] = memindex;
			++memindex;
		}
	}

	// end of example values for y

	return y;
}

/********************************************************************************************************
* END of Calculating the forward problem y = f(x)														*
*********************************************************************************************************/

// END of 2-1 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// 2-2 -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/********************************************************************************************************
* Construct a linear approximation of f by calculating the regression coefficients						*
*********************************************************************************************************/

vector<vector<double > > func_calculate_regression_parameters(int m
	, int l
	, int n
	, vector<vector<double > > X
	, vector<vector<double > > Y
	)
{
	/****************************************************************************************************
	* This function calculates the regression coefficient b considering the linear regression formulas  *
	* of																								*
	* https://de.wikipedia.org/wiki/Lineare_Regression#Sch.C3.A4tzung_der_Regressionskoeffizienten_2/	*
	* formula: b = (X^T * X)^(-1) * X^T * Y																*
	*****************************************************************************************************/
	vector<vector<double > > X_tilde(m + 1, vector<double>(l));

	vector<vector<double > > A_tilde(n, vector<double>(m+1));

	vector<vector<double > > Xtranspose(l, vector<double>(m+1));
	vector<vector<double > > Xmultiplication(m+1, vector<double>(m+1));
	vector <vector<double > > Xinverse(m+1, vector<double>(m+1));
	double* Xmultiplication_double_array = new double[m+1*m+1];
	vector<vector<double > > Xinverse_mult_X(m+1, vector<double>(l));

	vector<vector<double > > y_vector(1, vector<double>(l));
	vector<vector<double > > y_vector_transpose(l, vector<double>(1));

	vector<vector<double > > Xinverse_mult_X_mult_Y(m+1, vector<double>(1)); // Column vector	
	vector<vector<double > > Xinverse_mult_X_mult_Y_transpose(1, vector<double>(m+1));

	/*****************************************************************************
	* create X~ out of X according to (2.22) on page B23						 *
	******************************************************************************/

	X_tilde = func_create_x_tilde(X, m, l);

	//matrix_output(X_tilde, m + 1, l, "xtilde", "ende xtilde");
	
	/*****************************************************************************
	* calculate the transpose of matrix X and save it i Xtranspose				 *
	******************************************************************************/

	Xtranspose = func_matrix_transpose(X_tilde, m+1, l);

	//matrix_output(Xtranspose, l, m+1, "Xtranspose", "ende Xtranspose");

	/*****************************************************************************
	* multiply matrix Xtranspose with matrix X and save it in Xmultiplication	 *
	******************************************************************************/

	Xmultiplication = func_matrix_multiplication(X_tilde, Xtranspose,m+1, l, l, m+1);
	//matrix_output(Xmultiplication, m+1, m+1, "Xmultiplication", "ende Xmultiplication");
	
	/*****************************************************************************
	* Convert the vector<vector<double > > matrix Xmultiplication into a double*	 *
	* array. This is necessary because the matrix inversion function needs a	 *
	* double* array.															 *
	* It is saved to Xmultiplication_double_array								 *
	******************************************************************************/

	Xmultiplication_double_array = func_vector_matrix_to_double_star(Xmultiplication, m+1, m+1);
	//cout << "mukt\n";
	//printf("%f %f %f\n", Xmultiplication_double_array[0], Xmultiplication_double_array[1], Xmultiplication_double_array[2]);
	//printf("%f %f %f\n", Xmultiplication_double_array[3], Xmultiplication_double_array[4], Xmultiplication_double_array[5]);
	//printf("%f %f %f\n", Xmultiplication_double_array[6], Xmultiplication_double_array[7], Xmultiplication_double_array[8]);
	
	/*****************************************************************************
	* Calculate the inverse of this matrix in double* array form				 *
	* The result of the inverse is saved to Xmultiplication_double_array		 *
	******************************************************************************/

	inverse(Xmultiplication_double_array, m+1);
	
	/*****************************************************************************
	* The inverse double* Xmultiplication_double_array is saved to Xinverse by	 *
	* converting the double* array back to a vector<vector<double > > matrix		 *
	******************************************************************************/
	Xinverse = func_double_star_to_vector_matrix(Xmultiplication_double_array, m+1, m+1);
	//matrix_output(Xinverse, m + 1, m + 1, "Xinverse", "ende Xinverse");
	
	/*****************************************************************************
	* The inverse matrix is now multiplied by Xtranspose, the result is saved to *
	* Xinverse_mult_XT															 *
	******************************************************************************/
	Xinverse_mult_X = func_matrix_multiplication(Xinverse, X_tilde, m+1, m+1, m+1, l);
	//matrix_output(Xinverse_mult_X, m + 1, l, "Xinverse_mult_X", "ende Xinverse_mult_X");
	
	/*****************************************************************************
	* calculate (2.23) on page B23 for each n									 *
	******************************************************************************/

	for (int i = 0; i < n; i++)
	{
		/*****************************************************************************
		* take row vector i of n x l matrix Y										 *
		******************************************************************************/

		y_vector = func_take_vector_of_matrix(Y,n,l,i,"row"); // 1 x l
		//matrix_output(y_vector,  1, l, "y_vector", "ende y_vector");

		/*****************************************************************************
		* create transpose of of vector i											 *
		******************************************************************************/

		y_vector_transpose = func_matrix_transpose(y_vector, 1, l); // l x 1
		//matrix_output(y_vector_transpose, l, 1, "y_vector_transpose", "ende y_vector_transpose");

		/*****************************************************************************
		* Xinverse_mult_XT is now multiplied by Y. The result is saved to			 *
		* Xinverse_mult_XT_mult_Y.													 *
		******************************************************************************/

		Xinverse_mult_X_mult_Y = func_matrix_multiplication(Xinverse_mult_X, y_vector_transpose, m+1, l, l, 1);
		//matrix_output(Xinverse_mult_X_mult_Y, m+1, 1, "Xinverse_mult_X_mult_Y", "ende Xinverse_mult_X_mult_Y");

		/*****************************************************************************
		* m x 1 vector was created, take transpose									 *
		******************************************************************************/

		Xinverse_mult_X_mult_Y_transpose = func_matrix_transpose(Xinverse_mult_X_mult_Y, m+1, 1); // 1 x m+1
		//matrix_output(Xinverse_mult_X_mult_Y_transpose, 1, m+1, "Xinverse_mult_X_mult_Y_transpose", "ende Xinverse_mult_X_mult_Y_transpose");

		/*****************************************************************************
		* The created 1 x m+1 vector will be added to the matrix A to get a complete *
		* matrix A ( n x m+1 ) with all regression coefficients						 *
		******************************************************************************/

		A_tilde = func_add_vector_to_matrix(Xinverse_mult_X_mult_Y_transpose, 1, m+1, i, "row", A_tilde);
		//matrix_output(A_tilde, n, m + 1, "A_tilde", "ende A_tilde");
	}

	/*****************************************************************************
	* A is the result of this calculation of the regression coefficient			 *
	******************************************************************************/
	return A_tilde;

	/****************************************************************************************************
	* END of function calculating regression coefficient												*
	*****************************************************************************************************/
}

/********************************************************************************************************
* END of Constructing a linear approximation of f by calculating the regression coefficients			*
*********************************************************************************************************/

// END of 2-2 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// 2-3 -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/********************************************************************************************************
* In this section an update vector s will be created according to 2.2.2. on page B24 and 2-3 on page	*
* B21.																									*
*********************************************************************************************************/

vector<vector<double > > func_extact_A_from_A_tilde(vector<vector<double > > A_tilde, int n, int m)
{
	/****************************************************************************************************
	* takes A out of A_tilde																			*
	*****************************************************************************************************/

	vector<vector<double > > A(n,vector<double> (m));

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			A[i][j] = A_tilde[i][j];
		}
	}

	return A;

	/****************************************************************************************************
	* END of taking A out of A_tilde																	*
	*****************************************************************************************************/
}

vector<vector<double > > func_extract_y0_from_A_tilde(vector<vector<double > > A_tilde, int n, int m)
{
	/****************************************************************************************************
	* takes out y0 out of A_tilde																		*
	*****************************************************************************************************/

	vector<vector<double > > y0(n,vector<double>(1));

	for (int i = 0; i < n; i++)
	{
		y0[i][0] = A_tilde[i][m];
	}

	return y0;

	/****************************************************************************************************
	* END of taking y0 out of A_tilde																	*
	*****************************************************************************************************/
}

vector<vector<double > > func_find_update_vector(double* xDach, int m, int n, int l, vector<vector<double > > A, vector<vector<double > > ystar_matrix, vector<vector<double > > X, vector<vector<double > > y0)
{
	/****************************************************************************************************
	* This function creates an update vector s by using (2.25) as the formula for this calculation		*
	* There is a difference between the formula (2.25) and the actual use in this function:				*
	* Instead of a single vector of y* and x and subsequent iteration from 1 to l I use a the n x l		*
	* matrix of all y* and the whole m x l matrix of X. For using this I also use the n x l matrix of	*
	* y0. In each row of the y0 matrix every element is the same, but they differ from row to row		*
	* according to the change of n.																		*
	* So the result of this function is no update VECTOR but an update matrix. For using the vector		* 
	* the corresponding column vector of the update matrix has to be used.								*
	*****************************************************************************************************/

	vector<vector<double > > diag_x_dach(m,vector<double>(m));
	vector<vector<double > > vector_xDach(1,vector<double>(m));
	vector<vector<double > > square_diag_x_dach(m, vector<double>(m));
	vector<vector<double > > A_transpose(m, vector<double>(n));
	vector<vector<double > > A_sq_diag_mult_A_transpose(m, vector<double>(n));

	vector<vector<double > > A_mult_sq_diag(n, vector<double>(m));
	vector<vector<double > > A_mult_sq_diag_mult_Atranspose(n, vector<double>(n));
	double* A_mult_sq_diag_mult_Atranspose_double_star = new double[n*n];
	vector<vector<double > > A_mult_sq_diag_mult_Atranspose_inversion(n, vector<double>(n));

	vector<vector<double > > A_mult_X(n, vector<double>(l));
	vector<vector<double > > third_part_difference(n, vector<double>(l));
	vector<vector<double > > y0_matrix(n, vector<double>(l));

	vector<vector<double > > s(m, vector<double>(l));

	/****************************************************************************************************
	* convert double* array vector to vector<vector<double>>(1,vector<double>(m)) matrix vector			*
	*****************************************************************************************************/
	
	vector_xDach = func_double_star_to_vector_matrix(xDach, 1, m);
	
	/****************************************************************************************************
	* create diagonal matrix on basis of xDach vector													*
	*****************************************************************************************************/
	diag_x_dach = func_vector_to_diag_matrix(vector_xDach, 1, m, "row");

	/****************************************************************************************************
	* calculate square of diag xDach matrix																*
	*****************************************************************************************************/
	square_diag_x_dach = func_matrix_multiplication(diag_x_dach, diag_x_dach, m, m, m, m);

	/****************************************************************************************************
	* calculate the transpose of matrix A																*
	*****************************************************************************************************/
	A_transpose = func_matrix_transpose(A, n, m);

	/****************************************************************************************************
	* calculate the prodcut of (diag(xDach))² * A^(k)^T													*
	*****************************************************************************************************/
	A_sq_diag_mult_A_transpose = func_matrix_multiplication(square_diag_x_dach, A_transpose, m, m, m, n);

	/****************************************************************************************************
	* calculate second part of formula: (A^(k) * (diag(xDach))² * A^(k)^T)^(-1)							*
	* At first calculate: A^(k) * (diag(xDach))²														*
	*****************************************************************************************************/
	A_mult_sq_diag = func_matrix_multiplication(A, square_diag_x_dach, n, m, m, m);

	/****************************************************************************************************
	* calculate the rest of second part of the formula except the inversion								*
	* A^(k) * (diag(xDach))² * A^(k)^T																	*
	*****************************************************************************************************/
	A_mult_sq_diag_mult_Atranspose = func_matrix_multiplication(A_mult_sq_diag, A_transpose, n, m, m, n);

	/****************************************************************************************************
	* calculate inversion according to formula: (A^(k) * (diag(xDach))² * A^(k)^T)^(-1)					*
	*****************************************************************************************************/

	// convert A_mult_sq_diag_mult_Atranspose into double* array matrix
	A_mult_sq_diag_mult_Atranspose_double_star = func_vector_matrix_to_double_star(A_mult_sq_diag_mult_Atranspose, n, n);

	// do the inversion, inverse matrix is saved in A_mult_sq_diag_mult_Atranspose_double_star as double* array
	inverse(A_mult_sq_diag_mult_Atranspose_double_star, n);

	// convert double* array matrix into vector<vector<double > > matrix
	A_mult_sq_diag_mult_Atranspose_inversion = func_double_star_to_vector_matrix(A_mult_sq_diag_mult_Atranspose_double_star, n, n);

	/****************************************************************************************************
	* third part of the formula: (y* - A^(k) * X^(k) - y0^(k))											*
	* at first calculate product A * X																	*
	*****************************************************************************************************/
	A_mult_X = func_matrix_multiplication(A, X, n, m, m, l);

	/****************************************************************************************************
	* create matrix out of y0 col vector																*
	*****************************************************************************************************/
	for (int i = 0; i < l; i++)
	{
		y0_matrix = func_add_vector_to_matrix(y0, n, 1, i, "col", y0_matrix);
	}

	/****************************************************************************************************
	* calcuate difference y* - A^(k) * X^(k) - y0^(k)													*
	*****************************************************************************************************/
	third_part_difference = func_matrix_difference(ystar_matrix, A_mult_X, n, l);
	third_part_difference = func_matrix_difference(third_part_difference, y0_matrix, n, l);

	/****************************************************************************************************
	* calculate product of first, second and third part of the formula									*
	*****************************************************************************************************/
	s = func_matrix_multiplication(A_sq_diag_mult_A_transpose, A_mult_sq_diag_mult_Atranspose_inversion, m, n, n, n);
	s = func_matrix_multiplication(s, third_part_difference, m, n, n, l);

	return s;

	/****************************************************************************************************
	* END of creating update vector/matrix																*
	*****************************************************************************************************/

}

/********************************************************************************************************
* end of creating an update vector																		*
*********************************************************************************************************/

// END of 2-3 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// 2-4 -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/********************************************************************************************************
* diminish the update vector by its half until the sum of x + s is in the domain of f.					*
* Subsequent creating the new x values by adding the update vector										*
*********************************************************************************************************/

vector<vector<double > > func_shrink_vector(vector<vector<double > > input_s, vector<vector<double > > X, double* min_Range, double* max_Range, double factor, int n, int m, int l)
{
	vector<vector<double > > output_s(m, vector<double>(l));
	vector<vector<double > > working_vector(m, vector<double>(1));
	bool bool_half;
	output_s = input_s;

	/****************************************************************************************************
	* go trough all columns of input_vector s															*
	*****************************************************************************************************/
	for (int i = 0; i < l; i++)
	{
		/************************************************************************************************
		* go through all rows of input vector s to save the current column in the "working_vector"		*
		*************************************************************************************************/

		for (int o = 0; o < m; o++)
		{
			working_vector[o][0] = output_s[o][i];
		}
		
		/************************************************************************************************
		* END of saving current column to working_vector												*
		*************************************************************************************************/

		/************************************************************************************************
		* analyse if any sum of x + s-element in the current column is out of range, if so, the			*
		* the whole vector will be multiplied by 0.5													*
		*************************************************************************************************/
		bool_half = true;

		while (bool_half == true)
		{
			bool_half = false;

			/********************************************************************************************
			* check for each element in the current column if the sum of x+s is out of range			*
			*********************************************************************************************/

			for (int j = 0; j < m; j++)
			{
				//matrix_output(working_vector, m, 1, "working_vector", "end working vector");
				/*cout << "\n\n";
				cout << X[j][i];
				cout << " ";
				cout << working_vector[j][0];
				cout << " ";
				cout << max_Range[j];
				cout << " ";
				cout << min_Range[j];
				cout << " \n";*/
				if ((X[j][i] + working_vector[j][0]) > max_Range[j] || (X[j][i] + working_vector[j][0]) < min_Range[j])
				{
					bool_half = true;
				}
			}

			/********************************************************************************************
			* END of checking if sum is out of range													*
			*********************************************************************************************/

			/********************************************************************************************
			* multiply working vector by 0.5															*
			*********************************************************************************************/

			if (bool_half == true)
			{
				working_vector = func_matrix_scalar_multiplication(factor, working_vector, m, 1);
			}

			/********************************************************************************************
			* END of multiplying working vector by 0.5													*
			*********************************************************************************************/
		}

		/************************************************************************************************
		* END of analysing if there is need to multiply by 0.5											*
		*************************************************************************************************/

		/************************************************************************************************
		* change values of s by using the shrinked s													*
		*************************************************************************************************/

		for (int o = 0; o < m; o++)
		{
			output_s[o][i] = working_vector[o][0];
		}

		/************************************************************************************************
		* END of changing values of s by using the shrinked s											*
		*************************************************************************************************/
	}

	return output_s;
}

/********************************************************************************************************
* END of updating x																						*
*********************************************************************************************************/

// END of 2-4 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/// Stage 2 Broyden's method

// 3 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/********************************************************************************************************
* Set up the initial Jacobian approximation for each point x											*
*********************************************************************************************************/

/********************************************************************************************************
* END of setting up initial Jacobian																	*
*********************************************************************************************************/

// END Of 3 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// 4-1 -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/********************************************************************************************************
* Solve forward problem																					*
*********************************************************************************************************/

/********************************************************************************************************
* END of solving forward problem																		*
*********************************************************************************************************/

// END of 4-1 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// 4-2 -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/********************************************************************************************************
* update Jacobian with Broyden'S method																	*
*********************************************************************************************************/
vector<vector<double > > func_update_Jacobian(vector<vector<vector<vector<double > > > > Js_K2, int memindexK2, vector<vector<vector<double > > > s_K2, vector<vector<double > > y, vector<vector<double> >random_points_near_ystar, int n, int m, int l, int xl)
{
	// random_points_near_ystar: n x l
	// y: n x l
	// Js_K2: (K2-(K1+1), vector<vector<vector<double > > >(l, vector<vector<double>>(n, vector<double>(m))));
	// s_K2: (K2-(K1+1), vector<vector<double>>(m, vector<double>(l))); 
	// memindexK2 current run index for the for loop of K2 (starts with 0)

	vector<vector<double > > J(n, vector<double>(m));
	vector<vector<double > > vector_s_k_minus_1();

	vector<vector<double > > y_difference(n, vector<double>(1));
	vector<vector<double > > y_vector(n, vector<double>(1));
	vector<vector<double > > y_star_vector(n, vector<double>(1));

	vector<vector<double > > s_vector(m, vector<double>(1));
	vector<vector<double > > s_vector_transpose(1, vector<double>(m));
	double s_vector_2_norm;

	vector<vector<double > > s_vector_division(1, vector<double>(m));

	vector<vector<double > > second_times_third_part(n, vector<double>(m));

	vector<vector<double > > resulting_J(n, vector<double>(m));

	//for (int xl = 0; xl < l; xl++)
	//{
		/****************************************************************************************
		* save Jacobian(k-1) in J.																*
		*****************************************************************************************/
	cout << "\nmemindexK2: ";
	cout << memindexK2;
	cout << "xl: ";
	cout << xl;
	cout << "\n";
		for (int i = 0; i < number_of_equations; i++)
		{
			for (int j = 0; j < number_of_parameters_from_copasi; j++)
			{
				J[i][j] = Js_K2[memindexK2-1][xl][i][j];
			}
		}

		/****************************************************************************************
		* END of saving Jacobian k-1 in J.														*
		*****************************************************************************************/

		/****************************************************************************************
		* calculate y.j^(k) - y*																*
		*****************************************************************************************/	

		for (int i = 0; i < n; i++)
		{
			// save y vector from y_matrix
			y_vector[i][0] = y[i][xl];

			// save y_star_vector from y_star_matrix
			y_star_vector[i][0] = random_points_near_ystar[i][xl];
		}

		// calculate y_difference: n x 1
		y_difference = func_matrix_difference(y_vector,y_star_vector,n,1);

		/****************************************************************************************
		* END of calculating y.j^(k) - y*														*
		*****************************************************************************************/

		/****************************************************************************************
		* calculate (s.j^(k-1))^T / || s.j^(k-1)||²												*
		*****************************************************************************************/

		// define s_vector from matrix s_K2
		for (int i = 0; i < m; i++)
		{
			s_vector[i][0] = s_K2[memindexK2 - 1][i][xl];
		}

		// calculate transpose of s_vector
		s_vector_transpose = func_matrix_transpose(s_vector, m, 1);

		// calculate 2 norm form s_vector
		s_vector_2_norm = func_2_norm(s_vector, m, 1, "col");

		// calculate division 
		s_vector_division = func_matrix_scalar_multiplication(s_vector_2_norm, s_vector_transpose, 1, m);

		/****************************************************************************************
		* END of calculating (s.j^(k-1))^T / || s.j^(k-1)||²									*
		*****************************************************************************************/

		/****************************************************************************************
		* multiply second and third part of the formula:										*
		* (y.j^(k) - y*) *  ((s.j^(k-1))^T / || s.j^(k-1)||²)									*
		*****************************************************************************************/
		second_times_third_part = func_matrix_multiplication(y_difference, s_vector_division, n, 1, 1, m);

		/****************************************************************************************
		* END of multiplying second and third part of the formula:								*
		* (y.j^(k) - y*) *  ((s.j^(k-1))^T / || s.j^(k-1)||²)									*
		*****************************************************************************************/

		/****************************************************************************************
		* Add this product to J																	*
		* J+ (y.j^(k) - y*) *  ((s.j^(k-1))^T / || s.j^(k-1)||²)								*
		* According to (2.13)
		*****************************************************************************************/

		resulting_J = func_matrix_addition(J, second_times_third_part, n, m);

		/****************************************************************************************
		* END of adding J to product															*
		*****************************************************************************************/

		return resulting_J;
	//}

}
/********************************************************************************************************
* END of updating Jacobian with Broyden's method														*
*********************************************************************************************************/

// END of 4-2 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// 4-3 -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/********************************************************************************************************
* find search direction vector s																		*
*********************************************************************************************************/

vector<vector<double > > func_calculate_search_direction_vector(double* xDach, vector<vector<double > > J, vector<vector<double > > y, vector<vector<double> >random_points_near_ystar, int n, int m, int l, int xl)
{
	// vector<vector<double > > J(vector<vector<double>>(number_of_equations, vector<double>(number_of_parameters_from_copasi)));
	// random_points_near_ystar: n x l
	// y: n x l

	vector<vector<double > > vector_xDach(1, vector<double>(m));
	vector<vector<double > > diag_x_dach(m, vector<double>(m));
	vector<vector<double > > square_diag_x_dach(m, vector<double>(m));
	vector<vector<double > > J_transpose(m, vector<double>(n));
	vector<vector<double > > diag_mult_J_transpose(m, vector<double>(n));

	vector<vector<double > > second_part_first(n, vector<double>(m));
	vector<vector<double > > second_part(n, vector<double>(n));

	double* second_part_double_array;

	vector<vector<double > > second_part_inverse(n, vector<double>(n));
	vector<vector<double > > y_vector(n, vector<double>(1));
	vector<vector<double > > y_star_vector(n, vector<double>(1));
	vector<vector<double > > third_part(n, vector<double>(1));
	vector<vector<double > > first_times_second(m, vector<double>(n));
	vector<vector<double > > first_times_second_times_third(m, vector<double>(1));

	/************************************************************************************************
	* convert double xDach to vector_xDach															*
	*************************************************************************************************/
	vector_xDach = func_double_star_to_vector_matrix(xDach, 1, m);

	/****************************************************************************************************
	* create diagonal matrix on basis of xDach vector													*
	*****************************************************************************************************/
	diag_x_dach = func_vector_to_diag_matrix(vector_xDach, 1, m, "row");

	/****************************************************************************************************
	* calculate square of diag xDach matrix																*
	*****************************************************************************************************/
	square_diag_x_dach = func_matrix_multiplication(diag_x_dach, diag_x_dach, m, m, m, m);

	/****************************************************************************************************
	* calculate transpose of J																			*
	*****************************************************************************************************/
	J_transpose = func_matrix_transpose(J, n, m);

	/****************************************************************************************************
	* calculate product of diag(xDach)² * J^T															*
	*****************************************************************************************************/
	diag_mult_J_transpose = func_matrix_multiplication(square_diag_x_dach, J_transpose, m, m, m, n);

	/****************************************************************************************************
	* second part of the formula: (J * (diag(xDach))² * J^T))^(-1)										*
	* first calculate J * (diag(xDach))²																*
	*****************************************************************************************************/
	second_part_first = func_matrix_multiplication(J, square_diag_x_dach, n, m, m, m);

	/****************************************************************************************************
	* next multiplication of second part: J * (diag(xDach))² * J^T)										*
	*****************************************************************************************************/
	second_part = func_matrix_multiplication(second_part_first, J_transpose, n, m, m, n);

	/****************************************************************************************************
	* convert vector<vector<double > > matrix to double array												*
	*****************************************************************************************************/
	second_part_double_array = func_vector_matrix_to_double_star(second_part, n,n);

	/****************************************************************************************************
	* calculate inverse of second part																	*
	*****************************************************************************************************/
	inverse(second_part_double_array,n);

	/****************************************************************************************************
	* convert inverse double array to vector<vector<double > > matrix										*
	*****************************************************************************************************/
	second_part_inverse = func_double_star_to_vector_matrix(second_part_double_array, n, n);

	/****************************************************************************************************
	* calculate third part (y* - y.j)																	*
	* first get vectors of y out of the y matrices..													*
	*****************************************************************************************************/
	for (int i = 0; i < n; i++)
	{
		// take vector from y_star for each l
		y_star_vector[i][0] = random_points_near_ystar[i][xl];
		// take vector from y for each l
		y_vector[i][0] = y[i][xl];
	}

	/****************************************************************************************************
	* calculate third part of formula (y* - y.j)														*
	*****************************************************************************************************/
	third_part = func_matrix_difference(y_star_vector, y_vector, n, 1);

	/****************************************************************************************************
	* first part times second part																		*
	*****************************************************************************************************/
	first_times_second = func_matrix_multiplication(diag_mult_J_transpose, second_part_inverse, m,n,n,n);

	/****************************************************************************************************
	* product times third part																			*
	*****************************************************************************************************/
	first_times_second_times_third = func_matrix_multiplication(first_times_second, third_part, m, n, n, 1);

	return first_times_second_times_third;
}

/********************************************************************************************************
* END of setting up initial Jacobian																	*
*********************************************************************************************************/

// END of 4-3 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// 4-4 -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/********************************************************************************************************
* diminish the search direction vector by its half until the sum of x + s is in the domain of f.		*
* Subsequent creating the new x values by adding the search direction vector							*
*********************************************************************************************************/

vector<vector<double > > func_shrink_vector_broyden(vector<vector<double > > input_s, vector<vector<double > > X, double* min_Range, double* max_Range, double factor, int n, int m, int l)
{
	vector<vector<double > > output_s(m, vector<double>(1));
	vector<vector<double > > working_vector(m, vector<double>(1));
	bool bool_half;
	output_s = input_s;

	/************************************************************************************************
	* go through all rows of input vector s to save the current column in the "working_vector"		*
	*************************************************************************************************/

	for (int o = 0; o < m; o++)
	{
		working_vector[o][0] = output_s[o][0];
	}

	/************************************************************************************************
	* END of saving current column to working_vector												*
	*************************************************************************************************/

	/************************************************************************************************
	* analyse if any sum of x + s-element in the current column is out of range, if so, the			*
	* the whole vector will be multiplied by 0.5													*
	*************************************************************************************************/
	bool_half = true;

	while (bool_half == true)
	{
		bool_half = false;

		/********************************************************************************************
		* check for each element in the current column if the sum of x+s is out of range			*
		*********************************************************************************************/

		for (int j = 0; j < m; j++)
		{
			//matrix_output(working_vector, m, 1, "working_vector", "end working vector");
			/*cout << "\n\n";
			cout << X[j][i];
			cout << " ";
			cout << working_vector[j][0];
			cout << " ";
			cout << max_Range[j];
			cout << " ";
			cout << min_Range[j];
			cout << " \n";*/
			if ((X[j][0] + working_vector[j][0]) > max_Range[j] || (X[j][0] + working_vector[j][0]) < min_Range[j])
			{
				bool_half = true;
			}
		}

		/********************************************************************************************
		* END of checking if sum is out of range													*
		*********************************************************************************************/

		/********************************************************************************************
		* multiply working vector by 0.5															*
		*********************************************************************************************/

		if (bool_half == true)
		{
			working_vector = func_matrix_scalar_multiplication(factor, working_vector, m, 1);
		}

		/********************************************************************************************
		* END of multiplying working vector by 0.5													*
		*********************************************************************************************/
	}

	/************************************************************************************************
	* END of analysing if there is need to multiply by 0.5											*
	*************************************************************************************************/

	/************************************************************************************************
	* change values of s by using the shrinked s													*
	*************************************************************************************************/

	for (int o = 0; o < m; o++)
	{
		output_s[o][0] = working_vector[o][0];
	}

	/************************************************************************************************
	* END of changing values of s by using the shrinked s											*
	*************************************************************************************************/

	return output_s;
}

/********************************************************************************************************
* END of updating x																						*
*********************************************************************************************************/

// END of 4-4 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/*******************************************************************************************************
* This is the main procedure.																		   *
********************************************************************************************************/

int main()
{

	/********************************************************************************************************
	* STAGE 1																								*
	*********************************************************************************************************/

	/********************************************************************************************
	* STarting with the typical value x-Dach which probably should be the initial numbers for	*
	* the parameters given by the user. This point contains all given numbers for the parameters*
	********************************************************************************************/

	/********************************************************************************************
	* Producing an initial cluster of different random points centered about the typical value  *
	* x-Dach																					*
	*********************************************************************************************/

	/********************************************************************************************
	* variable declaration																		*
	*********************************************************************************************/

	// example values (also needed for non-example
	double* min_range = new double[number_of_parameters_from_copasi];
	double* max_range = new double[number_of_parameters_from_copasi];
	double* xDach = new double[number_of_parameters_from_copasi];
	vector<double> ystar(number_of_equations);

	// example matrices for testing matrixaddition
	vector< vector<double> > matrixA_test(3, vector<double>(4));
	vector< vector<double> > matrixB_test(3, vector<double>(4));
	vector< vector<double> > matrixC_test(3, vector<double>(4));

	// example for testing diagonal vector
	vector<double> diag_elements(3);

	// Example of transposed matrix
	vector<vector<double > > transposed_matrix;

	// Example of testing matrix multiplication
	vector<vector<double > > multiplied_matrix;
	vector<vector<double > > matrixD_test(4, vector<double>(5));

	// Example of scalar matrix mulitplication
	vector<vector<double > > scalar_product_matrix(3, vector<double>(4));

	// Example of calculating determinante
	double determinante;
	vector<vector<double > > double_matrix_3x3(3, vector<double>(3));

	// getting random_point
	double** random_points = 0;
	vector<vector<double > > vector_matrix_random_points_x(number_of_parameters_from_copasi, vector<double>(number_of_random_points));
	// all X matrices are saved in this 3 dimensional array 
	// the first dimension of this matrix has to be K1+1 because of K1+1 is the number of iterations 
	vector<vector<vector<double > > > Xs_K1(K1+1, vector<vector<double> >(number_of_parameters_from_copasi, vector<double>(number_of_random_points)));
	// all s matrices are saved in this 3 dimensional array
	// this is the update vector from Stage 1
	vector<vector<vector<double > > > s_K1(K1+1, vector<vector<double> >(number_of_parameters_from_copasi, vector<double>(number_of_random_points)));

	// getting random points near y*
	vector<vector<double > > random_points_near_ystar;

	// getting y out of the forward problem stage 1, 2-1 
	vector<vector<double > > y;
	vector<vector<double > > x(number_of_parameters_from_copasi, vector<double>(number_of_random_points));

	/*// declaration of variables for regression coefficient calculation, stage 1, 2-2
	int x_rows = 4; // Punkte pro Parameter
	int x_cols = 1; // Sollte immer eins sein, da immer nur ein Paramter geschätzt werden soll, sonst werden die Parameter addiert und ergeben zusammen Y, aber das soll so ja nicht sein. Es soll ja eher so sein, dass ein Parameter einen Wert bekommt aus dem forwadrd Problem (Y) und dann soll für alle Punkte des Paramters der Wert geschätzt werden, der benötigt wird um Y zu erzeugen.  Evtl muss hier vllt auch 2 stehen, wenn der 2. Wert eine 1 bekommt, damit man einen y-achsenabschnitt mit schätzen lassen kann.
	int y_rows = 4; // Ergebnisse aus der Parameter
	vector<vector<double > > X(x_rows, vector<double>(x_cols));
	vector<vector<double > > Y(y_rows, vector<double>(1));
	vector<vector<double > > b(x_cols, vector<double>(1)); // Rows: x_cols; cols: 1*/


	// declaration of variables for regression coefficient calculation, stage 1, 2-2
	
	int x_rows = number_of_parameters_from_copasi; // m
	int x_cols = number_of_random_points; 
	int y_rows = number_of_equations; // n
	int y_cols = number_of_random_points; // l
	vector<vector<double > > X(x_rows, vector<double>(x_cols));
	vector<vector<double > > Y(y_rows, vector<double>(y_cols));
	vector<vector<double > > A_tilde(number_of_random_points, vector<double>(number_of_parameters_from_copasi+1));	

	// find update vector

	vector<vector<double > > s(number_of_parameters_from_copasi,vector<double>(number_of_random_points));
	vector<vector<double > > A(number_of_equations,vector<double>(number_of_parameters_from_copasi));
	vector<vector<double > > y0(number_of_equations,vector<double>(1));

	// Shrink update vector

	vector<vector<double > > shrinked_s(number_of_parameters_from_copasi, vector<double>(number_of_random_points));
	vector<vector<double > > X_K1_plus_1(number_of_parameters_from_copasi, vector<double>(number_of_random_points));

	// Stage 2
	// all Jacobians used in Stage 2 are saved here

  vector< vector< vector< vector< double > > > > Js_K2(K2-(K1+1)+1, vector<vector<vector<double > > >(number_of_random_points, vector<vector<double> >(number_of_equations, vector<double>(number_of_parameters_from_copasi))));
	
	// current Jacobian J
	vector<vector<double > > J(vector<vector<double> >(number_of_equations, vector<double>(number_of_parameters_from_copasi)));

	// all s matrices are saved in this 3 dimensional array
	// this is the search direction vector from Stage 2
	vector<vector<vector<double > > > s_K2(K2 - (K1 + 1)+1, vector<vector<double > > (number_of_parameters_from_copasi, vector<double>(number_of_random_points)));

	vector<vector<double > > s_search_direction_vector(number_of_parameters_from_copasi, vector<double>(1));

	// 4-4 divide by 0.5
	vector<vector<double > > working_vector_x(number_of_parameters_from_copasi, vector<double>(1));
	vector<vector<double > > shrinked_s_broyden(number_of_parameters_from_copasi, vector<double>(1));
	vector<vector<double > > X_K2_plus_1(number_of_parameters_from_copasi, vector<double>(1));
	// all X matrices are saved in this 3 dimensional array 
	// the first dimension of this matrix has to be K1+1 because of K1+1 is the number of iterations 
	vector<vector<vector<double > > > Xs_K2(K2 - (K1 + 1) + 1, vector<vector<double> >(number_of_parameters_from_copasi, vector<double>(number_of_random_points)));
	vector<vector<double > > vector_matrix_random_points_x_broyden(number_of_parameters_from_copasi, vector<double>(1));
	/********************************************************************************************
	* END of variable declaration																*
	*********************************************************************************************/

	/*************************************
	* example values					 *
	**************************************/
	
	min_range[0] = 1.28;
	min_range[1] = 4.1122;
	min_range[2] = 99.1;
	min_range[3] = 3.000009;
	min_range[4] = 1.0;

	max_range[0] = 5.01;
	max_range[1] = 8.6;
	max_range[2] = 100;
	max_range[3] = 4.5;
	max_range[4] = 100;

	xDach[0] = 2;
	xDach[1] = 5;
	xDach[2] = 99.9;
	xDach[3] = 4.2;
	xDach[4] = 5;

	ystar[0] = 4;
	ystar[1] = 7;
	ystar[2] = 99.2;
	ystar[3] = 3.2;
	ystar[4] = 89;
	ystar[5] = 50;
	
	/*************************************
	* END of example values				 *
	**************************************/

	/********************************************************************************************
	* produce array/cluster of random points "random_points" (two dimensional array)			*
	* first dimension: number_of_random_points to be produced, second dimension:				*
	* number_of_parameters_from_copasi which come from copasi									*
	* this is the x-vector from 1-1																*
	*********************************************************************************************/

	srand(time(NULL));

	random_points = func_get_random_points(min_range
		, max_range
		, xDach
		);
	
	// convert double** matrix to a vector<vector<double > > matrix to use it in further steps
	vector_matrix_random_points_x = func_double_star_star_to_vector_matrix(random_points, number_of_parameters_from_copasi, number_of_random_points);
	
	// Copy these x values to Xs_K1 to save them 
	for (int i = 0; i < number_of_parameters_from_copasi; i++)
	{
		for (int j = 0; j < number_of_random_points; j++)
		{			
			Xs_K1[0][i][j] = vector_matrix_random_points_x[i][j];
		}
	}
	/********************************************************************************************
	* end of producing random points															*
	*********************************************************************************************/

	/********************************************************************************************
	* producing random points near y*															*
	* this is y-vector in 1-2																	*
	*********************************************************************************************/

	random_points_near_ystar = func_get_random_vector_near_ystar(ystar);

	/********************************************************************************************
	* END of producing random points near y*													*
	*********************************************************************************************/

	/****************************************************************************************************************
	* go through the algorithm K1 times (number_of_iterations_K1)													*
	*****************************************************************************************************************/

	for (int k = 0; k <= K1; k++)
	{

		/********************************************************************************
		* 2-1 Solve forward problem for each randomly created point x.j					*
		* The function evaluations at for each point of the cluster in this step are	*
		* independent and the most time consuming steps. So for this reason her multiply*
		* CPU cores could be used to increase velocity of the algorihtm.				*
		*********************************************************************************/

			/****************************************************************************
			* The forward problem has to be solved by copasi. So for testing reasons	*
			* I just give some values for f(x).											*
			*****************************************************************************/
			
			// y: n x l
			y = calculate_forward_problem(vector_matrix_random_points_x, number_of_parameters_from_copasi, number_of_random_points, number_of_equations);
			
		/********************************************************************************************
		* END of 2-1 (forward problem)																*
		*********************************************************************************************/

		/********************************************************************************
		* 2-2 Construct linear approximation of f										*
		* the regression function uses the matrix of X^(k) and not the X~ ^(k) matrix	*
		* So matrix which will be included is a m x l matrix.							*
		*																				*
		* The whole <vector<vector<double > > matrix Y is also a parameter of the function* 
		*																				*
		* the function creates a vector<vector<double > > matrix A, which contains the	*
		* regression coefficients according to the matrix A^(k) on page B23.			*
		*********************************************************************************/

			// Frobeniusnorm: DLANGE( NORM, M, N, A, LDA, WORK ) 
			// nach http://www.math.utah.edu/software/lapack/lapack-d/dlange.html
			// WIRD NICHT BENUTZT, DA KEIN MATRIX PROBLEM , SONDERN EIN VEKTORENPROBLEM

		
			// example of calculating regression coefficients

			/*	
			int x_rows = number_of_parameters_from_copasi; // 5
			int x_cols = number_of_random_points; // 10
			int y_rows = number_of_random_points;  // 10
			int y_cols = 1;
			*/
			/*X[0][0] = 2;
			X[0][1] = 5;
			X[1][0] = 1;
			X[1][1] = 3;
				
			Y[0][0] = 0.2;
			Y[0][1] = 0.5;
			Y[1][0] = 0.1;
			Y[1][1] = 0.3;*/
			
			X[0][0] = 10;
			X[0][1] = 11;
			X[0][2] = 12;
			X[0][3] = 13;
			X[0][4] = 14;
			X[0][5] = 15;
			X[0][6] = 16;
			X[0][7] = 17;
			X[0][8] = 18;
			X[0][9] = 19;

			X[1][0] = 20;
			X[1][1] = 21;
			X[1][2] = 22;
			X[1][3] = 23;
			X[1][4] = 24;
			X[1][5] = 25;
			X[1][6] = 26;
			X[1][7] = 27;
			X[1][8] = 28;
			X[1][9] = 29;

			X[2][0] = 4;
			X[2][1] = 351;
			X[2][2] = 82;
			X[2][3] = 633;
			X[2][4] = 134;
			X[2][5] = 935;
			X[2][6] = 10036;
			X[2][7] = 7;
			X[2][8] = 9938;
			X[2][9] = 199;

			X[3][0] = 40;
			X[3][1] = 41;
			X[3][2] = 42;
			X[3][3] = 43;
			X[3][4] = 44;
			X[3][5] = 45;
			X[3][6] = 46;
			X[3][7] = 47;
			X[3][8] = 48;
			X[3][9] = 49;

			X[4][0] = 50;
			X[4][1] = 51;
			X[4][2] = 52;
			X[4][3] = 53;
			X[4][4] = 54;
			X[4][5] = 55;
			X[4][6] = 56;
			X[4][7] = 57;
			X[4][8] = 58;
			X[4][9] = 59;

			Y[0][0] = 0.10;
			Y[0][1] = 0.11;
			Y[0][2] = 0.12;
			Y[0][3] = 0.13;
			Y[0][4] = 0.14;
			Y[0][5] = 0.15;
			Y[0][6] = 0.16;
			Y[0][7] = 0.17;
			Y[0][8] = 0.18;
			Y[0][9] = 0.19;

			Y[1][0] = 0.20;
			Y[1][1] = 0.21;
			Y[1][2] = 0.22;
			Y[1][3] = 0.23;
			Y[1][4] = 0.24;
			Y[1][5] = 0.25;
			Y[1][6] = 0.26;
			Y[1][7] = 0.27;
			Y[1][8] = 0.28;
			Y[1][9] = 0.29;

			Y[2][0] = 0.30;
			Y[2][1] = 0.31;
			Y[2][2] = 0.32;
			Y[2][3] = 0.33;
			Y[2][4] = 0.34;
			Y[2][5] = 0.35;
			Y[2][6] = 0.36;
			Y[2][7] = 0.37;
			Y[2][8] = 0.38;
			Y[2][9] = 0.39;

			Y[3][0] = 0.40;
			Y[3][1] = 0.41;
			Y[3][2] = 0.42;
			Y[3][3] = 0.43;
			Y[3][4] = 0.44;
			Y[3][5] = 0.45;
			Y[3][6] = 0.46;
			Y[3][7] = 0.47;
			Y[3][8] = 0.48;
			Y[3][9] = 0.49;

			Y[4][0] = 0.50;
			Y[4][1] = 0.51;
			Y[4][2] = 0.52;
			Y[4][3] = 0.53;
			Y[4][4] = 0.54;
			Y[4][5] = 0.55;
			Y[4][6] = 0.56;
			Y[4][7] = 0.57;
			Y[4][8] = 0.58;
			Y[4][9] = 0.59;

			Y[5][0] = 0.60;
			Y[5][1] = 0.61;
			Y[5][2] = 0.62;
			Y[5][3] = 0.63;
			Y[5][4] = 0.64;
			Y[5][5] = 0.65;
			Y[5][6] = 0.66;
			Y[5][7] = 0.67;
			Y[5][8] = 0.68;
			Y[5][9] = 0.69;
			
			A_tilde = func_calculate_regression_parameters(number_of_parameters_from_copasi
				, number_of_random_points
				, number_of_equations
				, vector_matrix_random_points_x // X
				, Y
				);
			
			/*
			A = func_calculate_regression_parameters(2
				, 2
				, 2
				, X
				, Y
				);*/


			/*
			cout << "regression parameters: \n";

			for (int i = 0; i < number_of_equations; i++)
			{
				for (int j = 0; j < number_of_parameters_from_copasi+1; j++)
				{
					printf("%.20f", A_tilde[i][j]);
					cout << "  ";
				}
				printf("\n");
			}

			cout << "regression parameters ende \n\n";*/


			/*

			cout << "regression parameters: \n";

			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					printf("%.20f", A[i][j]);
					cout << "  ";
				}
				printf("\n");
			}

			cout << "regression parameters ende \n\n";*/

			// end of example calculating regression coefficients

		/********************************************************************************
		* END of 2-2 constructing linear approximation									*
		*********************************************************************************/

			/****************************************************************************
			* ask if k < K1 because no update step is necessary in the last iteration	*
			* step.																		*
			*****************************************************************************/
			matrix_output(vector_matrix_random_points_x, number_of_parameters_from_copasi, number_of_random_points, "SSSS X", "END SSSS X");
			if (k < K1)
			{
				/// COMMENT: Is this right? Is it really no more necessary or is the goal of the whole process to output the x values? 

				/********************************************************************************
				* 2-3 find an update vector s.j													*
				*********************************************************************************/

				/********************************************************************************
				* extract A and y0 from A_tilde													*
				*********************************************************************************/
	
				A = func_extact_A_from_A_tilde(A_tilde, number_of_equations, number_of_parameters_from_copasi);
				y0 = func_extract_y0_from_A_tilde(A_tilde, number_of_equations, number_of_parameters_from_copasi);
			
				/********************************************************************************
				* END of extracting A and y0 from A_tilde										*
				*********************************************************************************/
			
				s = func_find_update_vector(xDach, number_of_parameters_from_copasi, number_of_equations, number_of_random_points, A, random_points_near_ystar, vector_matrix_random_points_x, y0);
		
				matrix_output(s, number_of_parameters_from_copasi, number_of_random_points, "UPDATE VECTOR s", "end UPDATE VECTOR s");

				/********************************************************************************
				* END of finding an update vector												*
				*********************************************************************************/

				/********************************************************************************
				* 2-4 shrink update vector until it is in domain of f							*
				* Subsequently create new X by adding s											*
				*********************************************************************************/

				shrinked_s = func_shrink_vector(s, vector_matrix_random_points_x, min_range, max_range, 0.5, number_of_equations, number_of_parameters_from_copasi, number_of_random_points);

				X_K1_plus_1 = func_matrix_addition(vector_matrix_random_points_x, shrinked_s, number_of_parameters_from_copasi, number_of_random_points);

				matrix_output(shrinked_s, number_of_parameters_from_copasi, number_of_random_points, "shrinked_s", "END OF shrinked_s");
				matrix_output(X_K1_plus_1, number_of_parameters_from_copasi, number_of_random_points, "X_K1_plus_1", "END OF X_K1_plus_1");

				/********************************************************************************
				 * END of shrinking update vector												*
				 ********************************************************************************/

				 /********************************************************************************
				 * Save matrix of new X values to collecting matrix of Xs						*
				 *********************************************************************************/
				
				for (int i = 0; i < number_of_parameters_from_copasi; i++)
				{
					for (int j = 0; j < number_of_random_points; j++)
					{
						Xs_K1[k + 1][i][j] = X_K1_plus_1[i][j];
					}
				}
				
				/********************************************************************************
				* END of saving matrix to collecting matrix of X								*
				*********************************************************************************/

				/********************************************************************************
				* Save matrix of update vectors to collecting matrix of update vectors			*
				*********************************************************************************/
				
				for (int i = 0; i < number_of_parameters_from_copasi; i++)
				{
					for (int j = 0; j < number_of_random_points; j++)
					{
						s_K1[k][i][j] = shrinked_s[i][j];
					}
				}
				
				/********************************************************************************
				* END of Saving matrix of update vectors to collecting matrix of update vectors	*
				*********************************************************************************/

				/********************************************************************************
				* prepare for next iteration step												*
				*********************************************************************************/

				for (int i = 0; i < number_of_parameters_from_copasi; i++)
				{
					for (int j = 0; j < number_of_random_points; j++)
					{
						vector_matrix_random_points_x[i][j] = X_K1_plus_1[i][j];
					}
				}

				/********************************************************************************
				* END of preparing for next iteration step										*
				*********************************************************************************/

			} // if k < K1

			matrix_output(vector_matrix_random_points_x, number_of_parameters_from_copasi, number_of_random_points, "SSSS X", "END SSSS X");

			/****************************************************************************
			* ask if k < K1 because no update step is necessary in the last iteration	*
			* step.																		*
			*****************************************************************************/

	} // for k <= K1

	/****************************************************************************************************************
	* end of going through the algorithm K1 times																	*
	*****************************************************************************************************************/

	/****************************************************************************************************************
	* END OF STAGE 1																			  				    *
	*****************************************************************************************************************/

// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	/****************************************************************************************************************
	* STAGE 2	Broyden's method																					*
	*****************************************************************************************************************/

	/**************************************************************************************************************
	* 3: set up the initial Jacobian approximation																  *
	***************************************************************************************************************/
	
	for (int l = 0; l < number_of_random_points; l++)
	{
		for (int i = 0; i < number_of_equations; i++)
		{
			for (int j = 0; j < number_of_parameters_from_copasi; j++)
			{
				Js_K2[0][l][i][j] = A[i][j];
			}
		}
	}
	
	matrix_output(A, number_of_equations, number_of_parameters_from_copasi, "FFFFFFF A", "FFFFFFF A end"); 
	matrix_output(Js_K2[0][0], number_of_equations, number_of_parameters_from_copasi, "XXXXXXXXXXXXXXX JsK2", "END XXXXXXXXXXX J");
	matrix_output(Js_K2[0][1], number_of_equations, number_of_parameters_from_copasi, "XXXXXXXXXXXXXXX JsK2", "END XXXXXXXXXXX J");
	matrix_output(Js_K2[0][number_of_random_points-1], number_of_equations, number_of_parameters_from_copasi, "XXXXXXXXXXXXXXX JsK2", "END XXXXXXXXXXX J");

	/**************************************************************************************************************
	* END of 3: setting up the initial Jacobian approximation													  *
	***************************************************************************************************************/

	/**************************************************************************************************************
	* Go trough K+1 to K2 as k																					  *
	***************************************************************************************************************/

	int memindexK2 = 0;
	for (int k = K1 + 1; k <= K2; k++)
	{
		/************************************************************************************************************
		* go through all l																							*
		*************************************************************************************************************/
		for (int xl = 0; xl < number_of_random_points; xl++)
		{
			if (k > 2)
			{
				//break;
			}
			/****************************************************************************************************************
			* set up x values which can be used in Stage 2: vector_matrix_random_points_x_broyden							*
			* But only for k = K1+1 so in first iteration step.																*
			* In the rest of the iteration steps the values of vector_matrix_random_points_x_broyden are updated on the		*
			* bottom of this loop.																							*
			*****************************************************************************************************************/
			if (k == K1 + 1)
			{
				for (int i = 0; i < number_of_parameters_from_copasi; i++)
				{
					vector_matrix_random_points_x_broyden[i][0] = vector_matrix_random_points_x[i][xl];
				}
			}

			/**************************************************************************************************************
			* 4-1 solve forward problem for each x.j																	  *
			***************************************************************************************************************/

			/// solve forward problem with vector_matrix_random_points_x_broyden as X values (number_of_parameters_from_copasi x 1)
			/// use newly created y in this forward problem to go further in the algorithm

			/**************************************************************************************************************
			* END of 4-1																								  *
			***************************************************************************************************************/

			/****************************************************************************************************************
			* ask if k < K2 because in last step it isn't necessary anymore to create search direction vectors.				*
			*****************************************************************************************************************/

			/// Is that true? same reason as in Stage 1 "if (k < K1)".

			if (k < K2)
			{

				/**************************************************************************************************************
				* 4-2 for k != K1 + 1 update Jacobian using Broyden's method												  *
				***************************************************************************************************************/

				if (k != (K1 + 1))
				{
					J = func_update_Jacobian(Js_K2, memindexK2, s_K2, y, random_points_near_ystar, number_of_equations, number_of_parameters_from_copasi, number_of_random_points, xl);
				}

				matrix_output(J, number_of_equations, number_of_parameters_from_copasi, "XXXXXXXXXXXXXXX J", "END XXXXXXXXXXX J");

				/**************************************************************************************************************
				* END of 4-2 search direction vector																		  *
				***************************************************************************************************************/

				/**************************************************************************************************************
				* 4-3 find search direction vector s.j																		  *
				***************************************************************************************************************/

				s_search_direction_vector = func_calculate_search_direction_vector(xDach, J, y, random_points_near_ystar, number_of_equations, number_of_parameters_from_copasi, number_of_random_points, xl);

				/**************************************************************************************************************
				* END of 4-3 search direction vector																		  *
				***************************************************************************************************************/

				/**************************************************************************************************************
				* 4-4 find new points approx. the solution set X* by updating X^(k)											  *
				***************************************************************************************************************/

					/********************************************************************************
					* 2-4 shrink update vector until it is in domain of f							*
					* Subsequently create new X by adding s											*
					*********************************************************************************/
					for (int row_m = 0; row_m < number_of_parameters_from_copasi; row_m++)
					{
						working_vector_x[row_m][0] = vector_matrix_random_points_x[row_m][xl];
					}

					shrinked_s_broyden = func_shrink_vector_broyden(s_search_direction_vector, working_vector_x, min_range, max_range, 0.5, number_of_equations, number_of_parameters_from_copasi, number_of_random_points);

					X_K2_plus_1 = func_matrix_addition(working_vector_x, shrinked_s, number_of_parameters_from_copasi, 1);

					matrix_output(shrinked_s_broyden, number_of_parameters_from_copasi, 1, "shrinked_s_broyden", "END OF shrinked_s_broyden");
					matrix_output(X_K2_plus_1, number_of_parameters_from_copasi, 1, "X_K2_plus_1", "END OF X_K2_plus_1");

					/********************************************************************************
					* END of shrinking search direction vector										*
					********************************************************************************/

					/********************************************************************************
					* Save matrix of new X values to collecting matrix of Xs						*
					*********************************************************************************/
					for (int i = 0; i < number_of_parameters_from_copasi; i++)
					{
						Xs_K2[k - (K1 + 1) + 1][i][xl] = X_K2_plus_1[i][0];
						//vector<vector<vector<double > > > Xs_K2(K2 - (K1 + 1) + 1, vector<vector<double>>(number_of_parameters_from_copasi, vector<double>(number_of_random_points)));
					}

					/********************************************************************************
					* END of saving matrix to collecting matrix of X								*
					*********************************************************************************/

					/********************************************************************************
					* Save matrix of search direction vectors to collecting matrix of update vectors*
					*********************************************************************************/

					for (int i = 0; i < number_of_parameters_from_copasi; i++)
					{
						s_K2[k - (K1 + 1)][i][xl] = shrinked_s_broyden[i][0];
					}
					
					/********************************************************************************
					* END of Saving matrix of update vectors to collecting matrix of update vectors	*
					*********************************************************************************/

					/********************************************************************************
					* prepare for next iteration step												*
					*********************************************************************************/

					for (int i = 0; i < number_of_parameters_from_copasi; i++)
					{
						vector_matrix_random_points_x_broyden[i][0] = X_K2_plus_1[i][0];
					}

					/********************************************************************************
					* END of preparing for next iteration step										*
					*********************************************************************************/

				/**************************************************************************************************************
				* END of 4-4 search direction vector																		  *
				***************************************************************************************************************/
			}

		}

		/************************************************************************************************************
		* END of going through all l																				*
		*************************************************************************************************************/

		++memindexK2;
		//MessageBox(NULL, "pjiJ", "ojio",MB_ICONWARNING | MB_CANCELTRYCONTINUE | MB_DEFBUTTON2);

	}

	/**************************************************************************************************************
	* END of Going trough K+1 to K2 as k																		  *
	***************************************************************************************************************/

	/****************************************************************************************************************
	* END OF STAGE 2																								*
	*****************************************************************************************************************/
	
	/*
	// example output of random points
		cout << "random points x: \n";
		for (int i = 0; i < number_of_parameters_from_copasi; i++)
		{
			for (int j = 0; j < number_of_random_points; j++)
			{
				printf("%.20f", random_points[i][j]);
				cout << "  ";
			}
			printf("\n");
		}	
		printf("\n");
		printf("\n");
		cout << "random points x vector matrix: \n";
		for (int i = 0; i < number_of_parameters_from_copasi; i++)
		{
			for (int j = 0; j < number_of_random_points; j++)
			{
				printf("%.20f", vector_matrix_random_points_x[i][j]);
				cout << "  ";
			}
			printf("\n");
		}
	
		printf("\n");
		printf("\n");
		printf("\n");
		printf("\n");
	// end of example random points



	// example of random values of ystar
		cout << "\n\nrandom GGGGGG ystars:\n\n";

		for (int i = 0; i < number_of_equations; i++)
		{
			for (int j = 0; j < number_of_random_points; j++)
			{
				printf("%.20f", random_points_near_ystar[i][j]);
				cout << "  ";
			}
			printf("\n");
		}
		cout << "\n\nend random GGGGGG ystars:\n\n";
	// end of example of ystar
	*/

	return 0;

}



