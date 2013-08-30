#include "stdafx.h"
#include "Plot.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

	// Write the content of 2 arrays to file for plotting in Matlab
	void WriteToFile(arr x, arr y, string path, int n)
	{
		ofstream outfile;
		outfile.open(path);

		for( int i = 0; i < n; i++ )
		{
			outfile << x[i] << ' ' << y[i] << endl;		// write line
		}

		outfile.close();
	}