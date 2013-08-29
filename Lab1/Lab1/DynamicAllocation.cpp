// Lab1.cpp : définit le point d'entrée pour l'application console.
//

#include "stdafx.h"
#include "DynamicAllocation.h"


/*! We will allocate dynamically a vector which address will be given with th argument.
 *  We'll get the number of columns among the argument too
 */
TYPE* dynamicalVector(int iNbColumn)
{
	TYPE* adVector;
	// We allocate the necessary memory space for the whole vector
	adVector = new TYPE[iNbColumn];
	// And then initialize the vector's coordinates
	for (int i=0; i< iNbColumn;i++)
	{
		adVector[i] = 0.0;
	}
	return adVector;
}

/*! We will allocate dynamically a matrix which address will be given with th argument
 *  We'll get the number of columns and of rows among the argument too
 */
TYPE** dynamicalMatrix(int iNbRow,int iNbColumn)
{
	// We allocate the necessary memory space for the whole matrix
	TYPE ** adMatrix = new TYPE*[iNbRow];
	for (int i =0; i<iNbRow;i++)
	{
		adMatrix[i] = new TYPE[iNbColumn];
	}
	
	//And we initialize our matrix
	for (int i=0; i<iNbRow; i++)
	{
		for (int j=0; j< iNbColumn;j++)
			adMatrix[i][j] = 0.0;
	}
	return adMatrix;
}