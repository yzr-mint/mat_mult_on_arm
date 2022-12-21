#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define SIZE 1024
#define BLOCK 64
#define FOR 10

#include <Eigen/Eigen>
#include <iostream>
#include <chrono>

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <arm_neon.h>
using namespace Eigen;
using namespace std;

Matrix<int, SIZE, SIZE, RowMajor> A1, B1, C1;

int N(SIZE), M(SIZE), L(SIZE);
using namespace std;

int main() {
	cout << "size: " << SIZE << endl;
	srand(time_t(NULL));
	//初始化两个二维数组
	for (int k = 0; k < L; ++k) {
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < N; ++j) {
				A1(i, k) = rand();
				B1(k, j) = rand();
			}
		}
	}

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < FOR; ++i) C1 += A1 * B1;
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	cout << "eigen: " << time_span.count()/FOR << "s" << endl;
}