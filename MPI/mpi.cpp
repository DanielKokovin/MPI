#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <Windows.h>
#include "mpi.h"
#include <clocale>
#include <iostream>
# define K 1600
using namespace std;


void init_matrix(double* X, double* Y, int N) {
	int i, j;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) {
			X[i * N + j] = Y[i * N + j] = i + j;
		}
}

void fill(double*& X, double*& Y, double*& Z, double*& res, int& N) {

	X = new double[N * N];
	Y = new double[N * N];
	Z = new double[N * N];
	res = new double[N * N];

	init_matrix(X, Y, N);
	for (int i = 0; i < N * N; i++) {
		Z[i] = 0;
		res[i] = 0;
	}
}

void matrix_multiply(double* X, double* Y, double* Z, int N, int nCommRunk, int nCommSize) {
	int i, j, k;
	for (i = 0; i < N; i++) {
		for (j = 0 + nCommRunk; j < N; j += nCommSize) {
			for (k = 0; k < N; k++) {
				Z[i * N + j] += X[i * N + k] * Y[k * N + j];
			}
		}
	}
}


int main(int argc, char* argv[]) {

	double* X;
	double* Y;
	double* Z;
	double* result = 0;
	int N = K;
	int nCommRunk, nCommSize, namelen, nCounter;
	int nIntervals;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	double t1, t2;

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &nCommRunk);
	MPI_Comm_size(MPI_COMM_WORLD, &nCommSize);
	MPI_Get_processor_name(processor_name, &namelen);

	if (nCommRunk == 0) {
		fill(X, Y, Z, result, N);
		t1 = MPI_Wtime();
		cout << "Multiplying matrices"<< endl;
		for (nCounter = 1; nCounter < nCommSize; nCounter++) {
			MPI_Send(&N, 1, MPI_INT, nCounter, 0, MPI_COMM_WORLD);
			MPI_Send(X, N * N, MPI_DOUBLE, nCounter, 1, MPI_COMM_WORLD);
			MPI_Send(Y, N * N, MPI_DOUBLE, nCounter, 2, MPI_COMM_WORLD);
			MPI_Send(Z, N * N, MPI_DOUBLE, nCounter, 3, MPI_COMM_WORLD);
		}
	}
	else {
		MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		X = new double[N * N];
		Y = new double[N * N];
		Z = new double[N * N];
		MPI_Recv(X, N * N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Y, N * N, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Z, N * N, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	matrix_multiply(X, Y, Z, N, nCommRunk, nCommSize);
	MPI_Reduce(Z, result, N * N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (nCommRunk == 0) {
		t2 = MPI_Wtime();
		cout << "Multipication time " << t2 - t1 << " seconds" << endl;
	}
	MPI_Finalize();
}