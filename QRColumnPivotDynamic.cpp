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

#include "QRColumnPivotDynamic.h"

#include "math.h"


// square matrix size
#ifndef QRSign
#define QRSign(x) ((x >= 0) ? (1) : (-1))
#endif // !QRSign



QRColumnPivotDynamic::QRColumnPivotDynamic(int squareMatSize) {	
	int i, j;
	matSize = squareMatSize;
	R = new double* [matSize];
	Ri = new double* [matSize];
	Q = new double* [matSize];
	P = new double* [matSize];
	RmatMulSquare = new double* [matSize];
	QmatMulSquare = new double* [matSize];
	AiColumnPivot = new double* [matSize];
	for (i = 0; i < matSize; i++) {
		R[i] = new double[matSize];
		Ri[i] = new double[matSize];
		Q[i] = new double[matSize];
		P[i] = new double[matSize];
		RmatMulSquare[i] = new double[matSize];
		QmatMulSquare[i] = new double[matSize];
		AiColumnPivot[i] = new double[matSize];
	}
	
	PFix = new double[matSize];
	colNorm = new double[matSize];
	V = new double[matSize];
	RmatMul = new double[matSize];
	QmatMul = new double[matSize];
	upperB = new double[matSize];
	Z = new double[matSize];

	for (i = 0; i < matSize; i++) {
		for (j = 0; j < matSize; j++) {
			R[i][j] = 0;
			Ri[i][j] = 0;
			Q[i][j] = 0;
			P[i][j] = 0;
			RmatMulSquare[i][j] = 0;
			QmatMulSquare[i][j] = 0;
			AiColumnPivot[i][j] = 0;
		}
		V[i] = 0;
		colNorm[i] = 0;
		RmatMul[i] = 0;
		QmatMul[i] = 0;
		upperB[i] = 0;
		Z[i] = 0;
		PFix[i] = 0;
	}

	maxNorm = 0;
	maxNormId = 0;
}


QRColumnPivotDynamic::~QRColumnPivotDynamic() {
	int i, j;
	for (i = 0; i < matSize; i++) {
		delete[] R[i];
		delete[] Ri[i];
		delete[] Q[i];
		delete[] P[i];
		delete[] RmatMulSquare[i];
		delete[] QmatMulSquare[i];
		delete[] AiColumnPivot[i];
	}
	delete[] R;
	delete[] Ri;
	delete[] Q;
	delete[] P;
	delete[] RmatMulSquare;
	delete[] QmatMulSquare;
	delete[] AiColumnPivot;

	delete[] PFix;
	delete[] colNorm;
	delete[] V;
	delete[] RmatMul;
	delete[] QmatMul;
	delete[] upperB;
	delete[] Z;
}


int QRColumnPivotDynamic::leastSquaresProblem(double** A, double* x, double* b, int len) {
	int i, j, k;
	// init
	maxNorm = 0;
	maxNormId = 0;
	for (i = 0; i < matSize; i++) {
		colNorm[i] = 0;
		for (j = 0; j < matSize; j++) {
			colNorm[i] += A[j][i] * A[j][i];	// square matrix, allow to exhange col and row to get column norm
			R[i][j] = A[i][j];
			Q[i][j] = 0;
			//P[i][j] = 0;
		}
		Q[i][i] = 1;
		//P[i][i] = 1;
		PFix[i] = i;
		if (colNorm[i] > maxNorm) {
			maxNorm = colNorm[i];
			maxNormId = i;
		}
	}
	// get Q R P
	for (i = 0; i < matSize - 1; i++) {
		columnProcess(i);
	}

	// upperB = Q^T * b
	for (i = 0; i < matSize; i++) {
		upperB[i] = 0;
		for (j = 0; j < matSize; j++) {
			upperB[i] += Q[j][i] * b[j];
		}
	}
	// R * z = upperB	for upper triangular matrix R use the backward substitution
	for (i = matSize - 1; i > 0; i--) {
		Z[i] = upperB[i] / R[i][i];
		for (j = 0; j <= i - 1; j++) {
			upperB[j] -= Z[i] * R[j][i];
		}
	}
	Z[0] = upperB[0] / R[0][0];
	for (i = 0; i < matSize; i++) {
		x[int(PFix[i])] = Z[i];
	}
	return 0;
}


int QRColumnPivotDynamic::pseudoInverseMat(double** A, double** Ai, int len) {
	int i, j, k;
	// init
	maxNorm = 0;
	maxNormId = 0;
	for (i = 0; i < matSize; i++) {
		colNorm[i] = 0;
		for (j = 0; j < matSize; j++) {
			colNorm[i] += A[j][i] * A[j][i];	// square matrix, allow to exhange col and row to get column norm
			R[i][j] = A[i][j];
			Ri[i][j] = 0;
			Q[i][j] = 0;
		}
		Q[i][i] = 1;
		PFix[i] = i;
		if (colNorm[i] > maxNorm) {
			maxNorm = colNorm[i];
			maxNormId = i;
		}
	}
	// get Q R P
	for (i = 0; i < matSize - 1; i++) {
		columnProcess(i);
	}
	// get upper matrix R inverse
	upperRInverse();
	// Ai*P = R^-1 * Q^T
	for (i = 0; i < matSize; i++) {
		for (j = 0; j < matSize; j++) {
			AiColumnPivot[i][j] = 0;
			for (k = 0; k < matSize; k++) {
				AiColumnPivot[i][j] += Ri[i][k] * Q[j][k];
			}
		}
	}
	// exchange row position after column pivot
	for (i = 0; i < matSize; i++) {
		for (j = 0; j < matSize; j++) {
			Ai[int(PFix[i])][j] = AiColumnPivot[i][j];
		}
	}

	return 0;
}


int QRColumnPivotDynamic::upperRInverse() {
	int i, j, k;
	double sum;
	// upper triangular matrix R use the backward substitution get inverse matrix
	for (i = 0; i < matSize; i++) {
		Ri[i][i] = 1 / R[i][i];
		for (j = i - 1; j >= 0; j--) {
			sum = 0;
			for (k = i; k > j; k--) {
				double t = R[j][k];
				double t1 = Ri[k][i];
				sum += R[j][k] * Ri[k][i];
			}
			Ri[j][i] = -sum / R[j][j];
		}
	}

	return 0;
}


int QRColumnPivotDynamic::columnProcess(int columnId) {
	int i, j, k;
	if (columnId < (matSize - 2)) {
		if (columnId != 0) {
			maxNorm = 0;
			maxNormId = 0;
			for (i = columnId; i < matSize; i++) {
				colNorm[i] -= R[columnId - 1][i] * R[columnId - 1][i];
				if (colNorm[i] > maxNorm) {
					maxNorm = colNorm[i];
					maxNormId = i;
				}
			}
		}
		if (maxNormId != columnId) {
			exchangeColumn(columnId, maxNormId);
		}
	}
	else {
		// last process column
		maxNorm = 0;
		maxNormId = columnId;
		maxNorm = colNorm[columnId] - R[columnId - 1][columnId] * R[columnId - 1][columnId];
	}


	getHouseHolderVector(columnId);
	// R = R - V * (V^T * R)
	// Q = Q - (Q * V) * V^T
	for (i = 0; i < matSize; i++) {
		RmatMul[i] = 0;
		QmatMul[i] = 0;
		for (j = 0; j < matSize; j++) {
			RmatMul[i] += V[j] * R[j][i];
			QmatMul[i] += Q[i][j] * V[j];
		}
		for (j = 0; j < matSize; j++) {
			RmatMulSquare[j][i] = 0;
			RmatMulSquare[j][i] = RmatMul[i] * V[j];
			QmatMulSquare[i][j] = 0;
			QmatMulSquare[i][j] = QmatMul[i] * V[j];
		}
	}
	for (i = 0; i < matSize; i++) {
		for (j = 0; j < matSize; j++) {
			R[i][j] -= RmatMulSquare[i][j];
			Q[i][j] -= QmatMulSquare[i][j];
		}
	}

	return 0;
}


int QRColumnPivotDynamic::exchangeColumn(int ori_columnId, int dest_columnId) {
	int i;
	double temp;
	temp = colNorm[ori_columnId];
	colNorm[ori_columnId] = colNorm[dest_columnId];
	colNorm[dest_columnId] = temp;
	temp = PFix[ori_columnId];
	PFix[ori_columnId] = PFix[dest_columnId];
	PFix[dest_columnId] = temp;
	for (i = 0; i < matSize; i++) {
		temp = R[i][ori_columnId];
		R[i][ori_columnId] = R[i][dest_columnId];
		R[i][dest_columnId] = temp;
		//temp = P[i][ori_columnId];
		//P[i][ori_columnId] = P[i][dest_columnId];
		//P[i][dest_columnId] = temp;
	}

	return 0;
}


int QRColumnPivotDynamic::getHouseHolderVector(int columnId) {
	double temp = 0;
	int i, j, k;

	for (i = 0; i < matSize; i++) {
		if (i >= columnId) {
			V[i] = R[i][columnId];
		}
		else {
			V[i] = 0;
		}
	}

	V[columnId] += QRSign(R[columnId][columnId]) * sqrt(maxNorm);
	for (i = 0; i < matSize; i++) {
		temp += V[i] * V[i];
	}

	if (temp > 0) {
		temp = sqrt(2 / temp);
		for (i = 0; i < matSize; i++) {
			V[i] *= temp;
		}
	}

	return 0;
}