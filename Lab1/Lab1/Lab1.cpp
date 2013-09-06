// Lab1.cpp : définit le point d'entrée pour l'application console.
//


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
bool forwardSubstitutionVector(arr v_f, arr v_b, arr v_a, int sizeVector);
bool backwardSubstitutionVector(arr v_f, arr v_b, arr v_a, arr v_Solution, int sizeVector);
int exA(int sizeVector);
int exB(int sizeVector);
TYPE maxRelError(arr numericalVector, arr analyticVector, int n);
arr analyticVector(int n, TYPE h);
TYPE u(TYPE x);

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
		// Then, we have to get the number of columns 
		//int iNbColumns;
		//cout << "Enter the number of columns you want \t";
		//cin >> iNbColumns;
#ifdef EXA
		exA(iNbRow);
#else 
		printf ("Part B \n");
		exB(iNbRow);
#endif
		printf( "\n Press q to leave \n" );	
		leavingKey = getchar();
	} while (leavingKey != 'q') ;

	mat A = randu<mat>(5,5);
	//cout << det(A); // Debugging Test

	getchar(); // pause to be able to see output while debugging

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
	
	TYPE h = ((1.0/(((TYPE) sizeVector) + 1.0))); // This is our step length

	// Initialization
	for (int i=0; i<sizeVector;i++)
	{
		v_f[i] = pow(h,2)*100*exp(-(TYPE)i*h);
		std::printf("f[%d]: %f \t", i, v_f[i]); // XXX A enlever ensuite. Juste pour tester.
		for (int j=0; j < sizeVector; j++)
		{
			if ((j== i-1) || (j== i+1))
				m_A(i, j) = -1;
			else if (j == i)
				m_A(i, j) = 2;
		}
	}

	std::printf("\n");// XXX A enlever ensuite. Juste pour tester.

	// Forward Substitution
	forwardSubstitutionMatrix(m_A,v_f,sizeVector);
	backwardSubstitutionMatrix(m_A,v_f,v_Solution,sizeVector);

	// A enlever plus tard: Vérif !
	TYPE truc = 0.0;
	for (int i=1; i<sizeVector+1;i++)
	{
		truc = (TYPE)i*h;
		TYPE total = u(truc);
		std::printf("U[%d] theorique : %f \t", i, total); // XXX A enlever ensuite. Juste pour tester.
	}
	std::printf("\n");// XXX A enlever ensuite. Juste pour tester.
	return 0;
}

// Function used to do the first step of the Gaussion elimination: the Forward Substitution
// We'll kill every term which prevents the matrix A to be an upper triangular one.
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
			printf("m_A[%d][%d]: %f \n", i, j, m_A(i, j));
		}
		v_f[i] -= v_f[i-1]*coef;
		
		printf("f[%d]: %f \n", i, v_f[i]); // XXX A enlever ensuite. Juste pour tester.
	}	
	return true;
}

// Function used to do the 2nd step of the Gaussian elimination
// We'll solve the equations' system
// TODO: prevent the case m_A[0][0]= 0
bool backwardSubstitutionMatrix(matr m_A, arr v_f, arr v_Solution, int sizeVector)
{
	// We process the n case first, and then, we loop
	v_Solution[sizeVector-1] = v_f[sizeVector-1] / m_A(sizeVector-1, sizeVector-1);
	printf("u[%d]: %f \n", sizeVector-1, v_Solution[sizeVector-1]);// XXX A enlever ensuite. Juste pour tester.
	for (int i= sizeVector-2; i > -1; i--)
	{
		v_Solution[i] = (v_f[i] - v_Solution[i+1]*m_A(i, i+1) / m_A(i, i));
		printf("u[%d]: %f \n", i, v_Solution[i]);
	}
	// Problème dans le coin ... >.<
	return true;
}
#pragma endregion 

int exB (int sizeVector)
{
	// Declaration of our vectors:
	arr v_Solution = arr(sizeVector); // This one is the one we are trying to find
	arr v_f = arr(sizeVector); // This one is the part h²*100e(-10x)
	arr v_b = arr(sizeVector); // The diagonal part of the matrix A
	arr v_a = arr(sizeVector); // And this one will be the one to describe the two sub-diagonals.
	
	TYPE h = (1.0/(((double)sizeVector) + 1.0)); // This is our step length

	// Initialization
	for (int i=0; i<sizeVector;i++)
	{
		v_f[i] = pow(h,2)*100*exp(-(double)i*h);
		v_b[i] = 2;
		v_a[i] = -1;
	}

	// Forward Substitution:
	forwardSubstitutionVector(v_f,v_b,v_a, sizeVector);
	backwardSubstitutionVector(v_f,v_b,v_a,v_Solution,sizeVector);
	
	// XXX :A enlever plus tard: Vérif !
	TYPE truc = 0.0;
	/*for (int i=1; i<sizeVector+1;i++)
	{
		truc = (double)i*h;
		TYPE total = (1 - (1-exp((double)-10))*truc - exp(-(double)10*truc));
		printf("U[%d] theorique : %f \t", i, total); // XXX A enlever ensuite. Juste pour tester.
	}*/

	// plot to file
	matr xy = matr(sizeVector, 3); // column 0: 0-1, column 1: analytical, column 2: numerical
	arr v_Analytic = analyticVector(sizeVector, h);

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

	return 0;
}

// Function used to do the first step of the Gaussion elimination: the Forward Substitution
// We'll kill every term which prevents the matrix A to be an upper triangular one.
bool forwardSubstitutionVector(arr v_f, arr v_b, arr v_c, int sizeVector)
{
	// First, we need to duplicate v_a, since we have to "nearly diagonals"
	arr v_a = arr(sizeVector);
	// Initialization of our 2nd vector
	for (int i= 0; i< sizeVector; i++)
		v_a[i]=v_c[i];
	
	for (int i=0; i< sizeVector; i++)
		printf("Before \tf%d : %f \n", i, v_f[i]);

	TYPE coef = 0.0;
	// And then, we compute the forward substitution:
	printf("a%d : %f \t b%d : %f \t c%d : %f \n",0, v_a[0],0, v_b[0],0, v_c[0]); // XXX : a enlever plus tard
	for (int i= 1; i< sizeVector; i++)
	{
		coef = v_a[i]/v_b[i-1];
		v_a[i] -= coef*v_b[i-1]; // XXX : à enlever plus tard .
		v_b[i] -= coef*v_c[i-1];
		v_f[i] -= coef*v_f[i-1];
		printf("a%d : %f \t b%d : %f \t c%d : %f \n",i, v_a[i],i, v_b[i],i, v_c[i]); // XXX : a enlever plus tard
	}

	for (int i=0; i< sizeVector; i++)// XXX : a enlever plus tard
		printf(" After \tf%d : %f \n", i, v_f[i]); // XXX : a enlever plus tard

	return true;
}

bool backwardSubstitutionVector(arr v_f, arr v_b, arr v_c, arr v_Solution, int sizeVector)
{
	// We save the last term:
	v_Solution[sizeVector-1] = v_f[sizeVector-1]/v_b[sizeVector-1];
	printf(" u%d : %f \t", sizeVector-1, v_Solution[sizeVector-1]); // XXX : a enlever plus tard
	// and then we compute what's left
	for (int i = sizeVector-1; i > 0; i--)
	{
		v_Solution[i-1] = (v_f[i-1] - v_c[i-1]*v_Solution[i])/ v_b[i-1]; // ~2n flops
		printf(" u%d : %f \t", i-1, v_Solution[i-1]); // XXX : a enlever plus tard
	}
	printf("\n"); 
	return true;
}

TYPE maxRelError(arr numericalVector, arr analyticVector, int n) {
	TYPE zeroValue = -1000; // the logarithm of anything should be greater than this
	TYPE e_max = zeroValue;

	for( int i = 0; i < n; i++ ) {
		TYPE v_i = numericalVector[i];
		TYPE u_i = analyticVector[i];

		TYPE e_i;
		try {
			TYPE e_i = log10(abs((v_i - u_i) / u_i)); // could be log10( 0 ) in theory...
		} catch( exception e ) {
			e_i = zeroValue; // ...but we can work around that.
		}

		if( e_i > e_max ) {
			e_max = e_i;
		}
	}

	return e_max;
}

arr analyticVector(int n, TYPE h) {
	arr a = arr(n);

	for( int i = 0; i < n; i++ ) {
		TYPE x = i * h;
		a[i] = u(x);
	}

	return a;
}

TYPE u(TYPE x) {
	return (1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x));
}