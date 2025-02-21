/*
MIT License

© 2025 Lao Chaode

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Author: Lao Chaode
Date: 2025-02-21
*/


#include <stdexcept>
#include <stdio.h>


#include "QRColumnPivotDynamic.h"

int main()
{
	double test_data[] = { -86.0615068994953, -0.0205959716556192, 0, 0,
	-255.507089802444, -3.53229983209699, -635.836267038566, 5.72169364259395E-21,
	-253.17974913249, -2.51170386044137, -2290.73769949779, -10.0998528315302,
	4.28903178185595E-24, 4.25499184707931E-26, -1653.90143245922, -11.0998528315302 };
	
	double J[4][4] = { 0 };
	double Ji[4][4] = { 0 };
	double* point_J[4];
	double* point_Ji[4];
	int k = 0;

	double* t;
	t = &test_data[8];
	double tt = t[0];

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			J[i][j] = test_data[k];
			k++;
		}
	}
	for (int i = 0; i < 4; i++) {
		point_J[i] = J[i];
		point_Ji[i] = Ji[i];
	}
	double b[4] = { 2.0492325907922506 ,
		-54188.500682648504,
		-188565.40871117834,
		 -134700.43885481835
	};
	double x[4] = { 0 };

	QRColumnPivotDynamic pd(4);

	pd.leastSquaresProblem(point_J, x, b, 4);   // solve Jx = b
	for(int i = 0; i < 4; i++) {
		printf("%f\n", x[i]);
	}
	pd.pseudoInverseMat(point_J, point_Ji, 4);  // get the pseudo inverse of J
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 4; j++) {
			printf("%f ", Ji[i][j]);
		}
		printf("\n");
	}
	
	return 0;
}

