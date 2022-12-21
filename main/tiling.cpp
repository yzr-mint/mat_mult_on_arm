/*
这是第一步。通过tiling，结果如下。
  T1       T2       2       4       6       8      10
 192  60.065s 60.499s 60.629s 59.999s 59.988s 59.976s
 196  60.084s 60.531s 60.686s 60.026s 60.022s 60.016s
 200  59.838s 60.339s 60.477s 59.818s 59.817s 59.812s
 204  59.459s 60.105s 60.342s 59.615s 59.599s 59.567s
 208  58.163s 58.658s 58.791s 58.11s  58.102s 58.095s
 212  58.205s 58.733s 58.89s  58.195s 58.195s 58.185s
 216  58.542s 59.11s  59.248s 58.532s 58.53s  58.53s
 220  58.306s 58.907s 59.063s 58.326s 58.328s 58.505s
 224  58.17s  58.762s 58.891s 58.125s 58.123s 58.118s
 228  58.47s  59.127s 59.268s 58.489s 58.486s 58.484s
 232  58.68s  62.596s 62.693s 61.844s 61.82s  62.036s
 236  62.024s 62.713s 62.753s 61.868s 61.858s 61.85s
 240  61.88s  62.91s  62.922s 62.05s  62.025s 62.005s
 244  62.26s  63.069s 63.149s 62.23s  62.224s 62.201s
 248  61.794s 62.587s 62.616s 61.751s 62.001s 62.046s
 252  61.488s 62.283s 62.215s 61.446s 61.405s 61.357s
*/
#include <iostream>
#include <chrono>
#include <iomanip>

#define MIN(x,y) (x<y?x:y)
#define MAX(x,y) (x>y?x:y)

#define SAME 512
#define M SAME
#define L SAME
#define N SAME
#define FOR 1
#define FROM 16
#define TO 256
#define STEP 4
int T1 = 1;//分块大小为T1的立方体
int T2 = 3;//将块矩阵分块为大小为T2的立方体

// m| l_ * l| n_
typedef int tp;
tp a_[M][L] = { 0 };
tp b_[L][N] = { 0 };
tp c_[M][N] = { 0 };

using namespace std;

void matmul(int a[M][L], int b[L][N], int c[M][N]) {
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++)
		{
			int s = 0;
			for (int k = 0; k < L; k++)
				s += a[i][k] * b[k][j];
			c[i][j] = s;
		}
	}
}

typedef int tp;
void matmul_0(tp a[M][L], tp b[L][N], tp c[M][N]) {
	register int i, j, k, km, im;
	register tp tmp;
	for (k = 0; k != L; ++k) {
		km = k * M;
		for (i = 0; i != M; ++i) {
			im = i * M;
			tmp = *((tp*)a + im + k);
			for (j = 0; j != N; ++j) {
				*((tp*)c + im + j) += tmp * *((tp*)b + km + j);
			}
		}
	}
}

void matmul_1(tp a[M][L], tp b[L][N], tp c[M][N]) {
	int kfrom, ifrom, jfrom, l, m, n, kto, ito, jto, it, jt, kt;
	register tp tmp;
	register int i, j, k, im, km;// kto, ito, jto, it, jt, kt;
	l = 1 + (L - 1) / T1;
	m = 1 + (M - 1) / T1;
	n = 1 + (N - 1) / T1;

	for (kt = 0; kt != l; ++kt) {
		kfrom = kt * T1;
		kto = MIN(kfrom + T1, L);
		for (it = 0; it != m; ++it) {
			ifrom = it * T1;
			ito = MIN(ifrom + T1, M);
			for (jt = 0; jt != n; ++jt) {
				jfrom = jt * T1;
				jto = MIN(jfrom + T1, N);

				for (k = kfrom; k != kto; ++k) {
					km = k * M;
					for (i = ifrom; i != ito; ++i) {
						im = i * M;
						tmp = *((tp*)a + im + k);
						for (j = jfrom; j != jto; ++j) {
							*((tp*)c + im + j) += tmp * *((tp*)b + km + j);
						}
					}
				}
			}
		}
	}
}

void matmul_2(int a[M][L], int b[L][N], int c[M][N]) {
	int ktt, ktfrom, ktto, itt, itfrom, itto, jtt, jtfrom, jtto, kfrom, ifrom, jfrom, l1, m1, n1, l, m, n, kto, ito, jto, it, jt, kt;
	register int tmp, i, j, k, im, km;// , kto, ito, jto, it, jt, kt;
	
	l = 1 + (L - 1) / T1;
	m = 1 + (M - 1) / T1;
	n = 1 + (N - 1) / T1;

	l1 = 1 + (l - 1) / T2;
	m1 = 1 + (m - 1) / T2;
	n1 = 1 + (n - 1) / T2;

	for (ktt = 0; ktt != l1; ++ktt) {
		ktfrom = ktt * T2;
		ktto = MIN(ktfrom + T2, l);
		for (itt = 0; itt != m1; ++itt) {
			itfrom = itt * T2;
			itto = MIN(itfrom + T2, m);
			for (jtt = 0; jtt != n1; ++jtt) {
				jtfrom = jtt * T2;
				jtto = MIN(jtfrom + T2, n);

				for (kt = ktfrom; kt != ktto; ++kt) {
					kfrom = kt * T1;
					kto = MIN(kfrom + T1, L);
					for (it = itfrom; it != itto; ++it) {
						ifrom = it * T1;
						ito = MIN(ifrom + T1, M);
						for (jt = jtfrom; jt != jtto; ++jt) {
							jfrom = jt * T1;
							jto = MIN(jfrom + T1, N);

							for (k = kfrom; k != kto; ++k) {
								km = k * M;
								for (i = ifrom; i != ito; ++i) {
									im = i * M;
									tmp = *((tp*)a + im + k);
									for (j = jfrom; j != jto; ++j) {
										*((tp*)c + im + j) += tmp * *((tp*)b + km + j);
									}
								}
							}
						}
					}
				}

			}
		}
	}
}

template<typename tp>
double countTime(void f(tp*, tp*, tp*), tp* a, tp* b, tp* c) {

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < FOR; ++i) f(a, b, c);
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

	/*for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			c[i * M + j] = 0;
		}
	}*/

	return time_span.count();
}

int main() {
	//初始化两个二维数组
	for (int k = 0; k < L; ++k) {
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < N; ++j) {
				a_[i][k] = rand();
				b_[k][j] = rand();
			}
		}
	}

	cout << "size: " << SAME << endl;
	//benchmark_les_mul
	std::cout << "banchmark:\t" << countTime(matmul, a_, b_, c_) / FOR << "s" << std::endl;

	std::cout << "After simple improvements:\t" << countTime(matmul_0, a_, b_, c_) << "s" << std::endl;
	cout << setw(4) << "T1" << " " << setw(4) << "T2" << setw(8) << "2" << setw(8) << "4" << setw(8) << "6" << setw(8) << "8" << setw(8) << "10" << endl;
	for (T1 = FROM; T1 != TO; T1 += STEP) {
		std::cout << setw(4) << T1 << "  " << setprecision(5) << countTime(matmul_1, a_, b_, c_) / FOR << "s";
		for (T2 = 2; T2 <= 10; T2 += 2) {
			std::cout << " " << setprecision(5) << countTime(matmul_2, a_, b_, c_) / FOR << "s";
		}
		cout << endl;
	}
	std::cout << "banchmark:\t" << countTime(matmul_0, a_, b_, c_) / FOR << std::endl;

}
