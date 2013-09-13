#include "stdafx.h"
#include "Plot.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

	// Write the content of a matrix to file for plotting in Matlab
	// (now redundant because armadillo does the same thing)
	void WriteToFile(TYPE** M, string path, int nRow, int nCol)
	{
		ofstream outfile; // file stream
		outfile.open(path); // attempt to open file for writing

		for( int i = 0; i < nRow; i++ ) // rows
		{
			for( int j = 0; j < nCol; j++ ) // columns
			{
				if( j > 0 ) // no leading space
				{
					outfile << ' ';
				}

				outfile << M[i][j];
			}

			outfile << endl; // newline character
		}

		outfile.close();
	}