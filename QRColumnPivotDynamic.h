
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


#pragma once


class QRColumnPivotDynamic {

private:

	int matSize;
	double** R;
	double** Ri;
	double** Q;
	double** P;
	double* PFix;

	double* colNorm;
	double maxNorm;
	double maxNormId;

	double* V;
	double* RmatMul;
	double** RmatMulSquare;
	double* QmatMul;
	double** QmatMulSquare;

	double* upperB;
	double* Z;
	double** AiColumnPivot;

	int upperRInverse();
	int columnProcess(int columnId);
	int exchangeColumn(int ori_columnId, int dest_columnId);
	int getHouseHolderVector(int columnId);

public:

	QRColumnPivotDynamic(int squareMatSize);
	~QRColumnPivotDynamic();

	int leastSquaresProblem(double** A, double* x, double* b, int len);
	int pseudoInverseMat(double** A, double** Ai, int len);

};

