#define NCS_IMPLEMENTATION
#include "../../NCS.h"

#include <chrono>
#include <iostream>

int main()
{
	double xt = 0.000000001, yt = 10.0, xc = 0.0000001, yc = 5.0;
	double tp = 0.0, nl = 1.1;
	double Vs = 2, Vt = 2;


	double gamma0 = 0.1;
	double Tc = 0.1;
	double mi = 0.1;

	auto start = std::chrono::high_resolution_clock::now();
	int num_samples = 1000;
	for (int i = 0; i < num_samples; i++) {
		double res = sc_integrate1D_gl_gpu(i_singlet_uifunt, 0.0, pi, Tc, mi, gamma0, tp);
		Tc += 1e-5;
		mi += 1e-6;
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()/num_samples <<"\n";

	return 0;
}