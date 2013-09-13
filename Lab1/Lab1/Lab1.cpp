#include "stdafx.h"
#include "DynamicAllocation.h"
#include "Plot.h"
#include <iostream>
#include <cmath>
#include <iomanip> 

using namespace std;
using namespace arma;

// Prototypes:
bool forwardSubstitutionMatrix(matr m_A, arr v_f, int sizeVector);
bool backwardSubstitutionMatrix(matr m_A, arr v_f, arr v_Solution, int sizeVector);
bool forwardSubstitutionVector(arr &v_f, arr &v_b, arr v_a, int sizeVector);
bool backwardSubstitutionVector(arr &v_f, arr &v_b, arr v_a, arr &v_Solution, int sizeVector);
int exA(int sizeVector);
TYPE exB(int sizeVector);
TYPE exB(int sizeVector, bool bIsPartOfD);
TYPE maxRelError(arr numericalVector, arr analyticVector, int n); // ExC
arr analyticVector(int n, TYPE h); // Ex C
TYPE u(TYPE x); // ExC
TYPE luCalling(int sizeVector); // ExD
int exD(int sizeVector); 
void exE(int sizeVector);
TYPE elapsedTime(clock_t start, clock_t finish);
void exCtest(int n);

// Main method
int main(int argc, char* argv[])
{
	// We have to get the size of our Matrix 
	// First, we have to get the number of rows
	char leavingKey = ' ';
	do
	{
		int iNbRow;
		cout << "Enter the number of rows you want \t";
		cin >> iNbRow;

		// Runs some or all exercises (depending on flags set in stdafx.h)
#ifdef EXA
		printf("Part A \n");
		exA(iNbRow);
#elif defined EXB 
		printf ("Part B \n");
		exB(iNbRow);
#elif defined EXC
		printf("Part C \n");
		exCtest(iNbRow);
#elif defined EXD
		printf("Part D \n");
		exD(iNbRow);
#elif defined EXE
		printf ("Part E \n");
		exE(iNbRow);
#endif
		printf( "\n Press q to leave \n" );	
		fflush(stdin);
		leavingKey = getchar();

	} while (leavingKey != 'q') ; // if 'q' then quit
	
	return 0;
}

#pragma region Exercise A
// This function will be the calling function for our a) part
// Reminder : We need 1 matrix and 2 vectors to solve the following equation: - (u(i+1) + u(i-1) -2u(i))/h² = 100e(-10x)
int exA(int sizeVector)
{
	// Declaration of our matrix
	matr m_A = matr(sizeVector, sizeVector);

	// Declaration of our vectors:
	arr v_Solution = arr(sizeVector); // This one is the one we are trying to find
	arr v_f = arr(sizeVector); // This one is the part h²*100e(-10x)
	
	TYPE h = (1.0/(((TYPE) sizeVector)-1.0)); // This is our step length

	// Initialization
	for (int i=0; i<sizeVector;i++)
	{
		v_f[i] = pow(h,2)*100*exp(-(TYPE)i*h);
		std::printf("f[%d]: %f \t", i, v_f[i]);
		for (int j=0; j < sizeVector; j++)
		{
			if ((j== i-1) || (j== i+1))
				m_A(i, j) = -1;
			else if (j == i)
				m_A(i, j) = 2;
		}
	}

	std::printf("\n");

	// Forward Substitution
	forwardSubstitutionMatrix(m_A,v_f,sizeVector);
	// Then backward substitution
	backwardSubstitutionMatrix(m_A,v_f,v_Solution,sizeVector);

	// plot to file - this can be imported into Matlab
	matr xy = matr(sizeVector, 3);
	arr v_Analytic = analyticVector(sizeVector, h);

	// column 0: x (0-->1), column 1: analytical, column 2: numerical
	for( int i = 0; i < sizeVector; i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			xy(i, 0) = (i * h);
			xy(i, 1) = v_Analytic[i];
			xy(i, 2) = v_Solution[i];
		}
	}

	xy.save("plot.txt", raw_ascii);

	printf(" This is our max Rel error : %d", maxRelError(v_Solution,v_Analytic,sizeVector));

	std::printf("\n");
	return 0;
}

// Function used to do the first step of the Gaussion elimination: the Forward Substitution
bool forwardSubstitutionMatrix(matr m_A, arr v_f, int sizeVector)
{
	TYPE coef = 0.0;
	// Our matrix is a tridiagonal matrix, thus, we don't need to do a classical Gauss elimination
	for (int i=1;i < sizeVector; i++)
	{
		coef = m_A(i, i-1)/m_A(i-1, i-1);
		for (int j=1; j< sizeVector; j++)
		{
			m_A(i, j) -= coef*m_A(i-1, j); 
			// printf("m_A[%d][%d]: %f \n", i, j, m_A(i, j)); // debug
		}
		v_f[i] -= v_f[i-1]*coef;
		
		// printf("f[%d]: %f \n", i, v_f[i]); // debug
	}	
	return true;
}

// Function used to do the 2nd step of the Gaussian elimination
// We'll solve the equations' system
bool backwardSubstitutionMatrix(matr m_A, arr v_f, arr v_Solution, int sizeVector)
{
	// We process the n case first, and then, we loop
	v_Solution[sizeVector-1] = v_f[sizeVector-1] / m_A(sizeVector-1, sizeVector-1);
	// printf("u[%d]: %f \n", sizeVector-1, v_Solution[sizeVector-1]); // debug
	for (int i= sizeVector-2; i > -1; i--)
	{
		v_Solution[i] = (v_f[i] - v_Solution[i+1]*m_A(i, i+1) / m_A(i, i));
		// printf("u[%d]: %f \n", i, v_Solution[i]); // debug
	}
	return true;
}
#pragma endregion 

#pragma region Exercise B
// This function is an overloading of the exB(int) function. But when calling exB from exD, we don't want to plot, or things like that. Thus,
// we'll use this boolean
TYPE exB(int sizeVector, bool bIsPartOfD)
{
	// Declaration of our vectors:
	arr v_Solution = arr(sizeVector); // This one is the one we are trying to find
	arr v_f = arr(sizeVector); // This one is the part h²*100e(-10x)
	arr v_b = arr(sizeVector); // The diagonal part of the matrix A
	arr v_a = arr(sizeVector); // And this one will be the one to describe the two sub-diagonals.
	
	TYPE h = (1.0/(((double)sizeVector) - 1.0)); // This is our step length
	TYPE x = 0.0f;
	// Initialization
	for (int i=0; i<sizeVector;i++)
	{
		x = i*h;
		v_f[i] = pow(h,2)*100*exp(-(double)10*x);
		v_b(i) = 2;
		v_a[i] = -1;
	}

	// In ex. D, we need timers
	clock_t start, finish;
	if( bIsPartOfD )
	{
		start = clock();
	}

	// Forward Substitution:
	forwardSubstitutionVector(v_f,v_b,v_a, sizeVector);
	backwardSubstitutionVector(v_f,v_b,v_a,v_Solution,sizeVector);

	if( bIsPartOfD ) // (ex. D) in this case, we are only interested in the time, not the results
	{
		finish = clock();
		return elapsedTime(start, finish);
	}
	else // (ex. B) we want to plot the results
	{
		// plot to file
		matr xy = matr(sizeVector, 3); // column 0: 0-1, column 1: analytical, column 2: numerical
	
		arr v_Analytic = analyticVector(sizeVector, h); // analytic solution

		for( int i = 0; i < sizeVector; i++ )
		{
			for( int j = 0; j < 3; j++ )
			{
				xy(i, 0) = (i * h);
				xy(i, 1) = v_Analytic[i];
				xy(i, 2) = v_Solution[i];
			}
		}

		xy.save("plot.txt", raw_ascii);

		// Calculate max relative error
		TYPE error = maxRelError(v_Solution,v_Analytic,sizeVector);
		printf(" This is our max Rel error : %f", error);

		cout << endl;
	}

	return 0.0;
}

// do ex. B with no timers
inline TYPE exB (int sizeVector)
{
	return exB(sizeVector, false);
}

// Function used to do the first step of the Gaussion elimination: the Forward Substitution
// We'll kill every term which prevents the matrix A to be an upper triangular one.
// Note the &'s needed to make Armadillo objects be passed by reference!
bool forwardSubstitutionVector(arr &v_f, arr  &v_b, arr v_c, int sizeVector)
{
	// First, we need to duplicate v_a, since we have to "nearly diagonals"
	arr v_a = arr(sizeVector);
	// Initialization of our 2nd vector
	for (int i= 0; i< sizeVector; i++)
		v_a[i]=v_c[i];
	
	TYPE coef = 0.0;
	// And then, we compute the forward substitution:
	for (int i= 1; i< sizeVector; i++)
	{
		coef = 1/v_b[i-1];
		//v_a[i] += coef*v_b[i-1];  // This operation just allows us to check if  every term of the first diago is null after the substitution
		v_b[i] -= coef;
		v_f[i] += coef*v_f[i-1];
		//printf("a%d : %f \t b%d : %f \t c%d : %f \n",i, v_a[i],i, v_b[i],i, v_c[i]); // debug
	}

	//for (int i=0; i< sizeVector; i++) // debug
		//printf(" After \tf%d : %f \n", i, v_f[i]); // debug

	return true;
}

// Perform backward substitution
// Note the &'s needed to make Armadillo objects be passed by reference!
bool backwardSubstitutionVector(arr &v_f, arr &v_b, arr v_c, arr &v_Solution, int sizeVector)
{
	// We save the last term:
	v_Solution[sizeVector-1] = v_f[sizeVector-1]/v_b[sizeVector-1];
	//printf(" u%d : %f \t", sizeVector-1, v_Solution[sizeVector-1]);
	// and then we compute what's left
	TYPE x = v_Solution[sizeVector-1];
	for (int i = sizeVector-1; i > 0; i--)
	{
		v_Solution[i-1] = (v_f[i-1] + v_Solution[i]) / v_b[i-1]; // ~2n flops
		x= v_Solution[i-1];
		//printf(" u%d : %f | ", i-1, v_Solution[i-1]); // Not displaying it... Pretty clear for little number of Row, but not in other cases
	}
	printf("\n"); 
	return true;
}
#pragma endregion

#pragma region Exercise C

// For debugging purposes only (n = 5 is a good test case)
void exCtest(int n)
{
	arr v_n = arr(n);
	arr v_a = arr(n);

	// initialize vectors
	// v_a is the "analytic" solution
	// v_n is the "numerical" solution, always off by 1
	for( int i = 0; i < n; i++ )
	{
		TYPE x = 100.0 * ((TYPE)i)/(((TYPE)n)-1.0);
		v_a[i] = x * (100.0 - x);
		v_n[i] = x * (100.0 - x) + 1;
	}

	/* debug
	for( int i = 0; i < n; i++ )
	{
		cout << v_a[i] << " " << v_n[i] << endl;
	} */

	// output max relative error and compare to known value
	TYPE epsilon = maxRelError(v_n, v_a, n);
	printf("epsilon = %f", epsilon);
}

// Calculate maximum relative error, given a vector of analytic values and a vector of numerical ones
TYPE maxRelError(arr numericalVector, arr analyticVector, int n) {
	const TYPE ZEROVALUE = -1000; // the logarithm of anything should be greater than this in practice
	TYPE e_max = ZEROVALUE;

	for( int i = 0; i < n; i++ ) {
		TYPE v_i = numericalVector[i];
		TYPE u_i = analyticVector[i];
		TYPE e_i = 0.0f;

		if (u_i != 0) // computing the relative error makes no sense if u_i is zero (we'd get infinity)
		{
			if (v_i == u_i) // perfect match!
			{
				e_i = ZEROVALUE;
			}
			else
			{
				e_i = log10(abs((v_i - u_i) / u_i));
			}
		}
		else
		{
			e_i = ZEROVALUE; // ignore this case
		}

		if (e_i > e_max)	
		{
			e_max = e_i; // this is our new maximum value

			// printf("i: %d | e_max: %f \t",  i,e_max); // debug
		}		
	}

	return e_max;
}

// Fill a vector with analytic values (takes vector length and step size as input)
arr analyticVector(int n, TYPE h) {
	arr a = arr(n);

	for( int i = 0; i < n; i++ ) {
		TYPE x = i * h; // assumes that x0 = 0
		a[i] = u(x);
	}

	return a;
}

// Analytic solution to the differential equation
inline TYPE u(TYPE x) {
	return (1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x));
}

#pragma endregion

#pragma region Exercise D
/*! This function allows us to time lu() and our tridiagonal solver functions */
int exD(int sizeVector)
{
	// time measurement
	clock_t start, finish;
	TYPE timeToFinish = 0.0;
	
	// Timing the Armadillo's lu function
	printf(" Starting LU decomposition ... \n");
	timeToFinish = luCalling(sizeVector);
	printf(" Elapsed time to do the LU decomposition : %f \n", timeToFinish);
	
	// Timing our tridiagonal solver
	printf(" Starting to use our tridiagonal solver ... \n");
	timeToFinish = exB(sizeVector,true);
	printf(" Elapsed time to solve the system with our algo. : %f \n", timeToFinish);
	
	return 0;
}

// Perform LU decomposition using Armadillo on a tridiagonal matrix
TYPE luCalling(int sizeVector)
{
	try
	{
		// We declare and instantiate our matrices
		matr m_A = matr(sizeVector,sizeVector);
		matr P,L,U;
		// Then we initialize them
		for (int i=0; i< sizeVector; i++)
		{
			for (int j=0; j < sizeVector;j++)
			{
				m_A(i, j) = 0; 
				if (j == i)
					m_A(i, j) = 2;
				else if ((j == i-1) || (j == i+1))
					m_A(i, j) = -1;
			}
		}

		// initialize timers
		clock_t start, finish;

		// time LU decomposition
		start = clock();
		bool success = lu(L, U, P, m_A); // returns false if something went wrong on Armadillo's side of things
		finish = clock();

		if (success)
		{
			return elapsedTime(start, finish);
		}
		else
		{
			printf("An error occured during the Lower Upper decomposition ...\n");
			printf("lu returned false!");

			return 0.0;
		}
	}
	catch( exception e ) // typical cause: out of memory (large matrices)
	{
		printf("An error occured during the Lower Upper decomposition ...\n");
		printf(e.what());

		return 0.0;
	}
}

#pragma endregion

#pragma region Exercise E

// returns number of seconds (as double/float) between two points in time
inline TYPE elapsedTime(clock_t start, clock_t finish)
{
	return (((TYPE)(finish - start))/((TYPE)CLOCKS_PER_SEC)); // avoid long division here
}

// performs matrix multiplication (using pointer matrixes)
// row major order
TYPE** rowMult(TYPE** A, TYPE** B, int Arows, int Acols, int Brows, int Bcols)
{
	// check dimensions
	if( Acols != Brows )
	{
		throw new exception("matrix multiplication undefined - dimensions do not match");
	}

	// else (A*B is defined):

		int k_max = Acols; // this dimension is cancelled in the multiplication

		// resulting n x m matrix
		int n = Arows;
		int m = Bcols;
		TYPE** C = dynamicalMatrix(n, m);

		// use definition of matrix multiplication
		for( int i = 0; i < n; i++ )
		{
			for( int j = 0; j < m; j++ )
			{
				TYPE sum = 0;
				for( int k = 0; k < k_max; k++ )
				{
					sum += A[i][k]*B[k][j];
				}
				C[i][j] = sum;
			}
		}

		return C;
}

// performs matrix multiplication (using pointer matrixes)
// column major order
TYPE** colMult(TYPE** A, TYPE** B, int Arows, int Acols, int Brows, int Bcols)
{
	// check dimensions
	if( Acols != Brows )
	{
		throw new exception("matrix multiplication undefined - dimensions do not match");
	}

	// else (A*B is defined):

		int k_max = Acols; // this dimension is cancelled in the multiplication

		// resulting n x m matrix
		int n = Arows;
		int m = Bcols;
		TYPE** C = dynamicalMatrix(n, m);
	
		// use definition of matrix multiplication
		for( int j = 0; j < m; j++ )
		{
			for( int i = 0; i < n; i++ )	
			{
				TYPE sum = 0;
				for( int k = 0; k < k_max; k++ )
				{
					sum += A[i][k]*B[k][j];
				}
				C[i][j] = sum;
			}
		}

		return C;
}

// compare computing time of multiplying matrices row major vs. column major
// NOTE: Will take up to an hour for n ~ 5000
void exE(int sizeVector)
{
	cout << "Initializing...";

	// initialize n x n matrices to multiply
	int n = sizeVector;
	TYPE** A = randomMatrix(n, n);
	TYPE** B = randomMatrix(n, n);
	TYPE** C;

	// time measurement
	clock_t start, finish;
	TYPE tRow, tCol;

	cout << "done!" << endl << "Row major...";

	// row major
	start = clock();
	C = rowMult(A, B, n, n, n, n);
	finish = clock();
	tRow = elapsedTime(start, finish);

	cout << "done! Time elapsed: " << tRow << endl << "Column major...";

	// column major
	start = clock();
	C = colMult(A, B, n, n, n, n);
	finish = clock();
	tCol = elapsedTime(start, finish);

	cout << "done! Time elapsed: " << tCol << endl;
}

#pragma endregion