#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <numeric>
#include <vector>
#include "cmath"
#include <string>
#include <cmath>
#include <cuda_runtime.h>
#include <string>
#include <fstream>

using namespace std;

#define BLOCK_SIZE 16
const double h_pi = 3.141592653597932384;
__constant__ double pi_device = 3.141592653597932384;
const double pi = 3.141592653597932384;
const double delta = 0.0000000001;
const int n = 10000, k = 5;

const double xk[] = { 0.9061798459386641,-0.9061798459386641,0.538469310105683,-0.5384693101056829,0.0 };
const double ak[] = { 0.2369268850561876,0.2369268850561876,0.47862867049936647,0.47862867049936586,0.5688888888888889 };

typedef double (*Under_Integral_Func)(double, double, double, double, double);

__device__ double g(double x, double gamma0)
{
	return gamma0 * (abs(sin(x)));
}

__device__ double ek(double x, double tp)
{
	return -2.0 * (cos(x)) + 2.0 * tp * cos(2 * x);
}

__device__ double sc_uifunc(double x, double Tc, double mi, double gamma0, double tp)
{
	double ksip = ek(x, tp) - mi + g(x, gamma0), ksim = ek(x, tp) - mi - g(x, gamma0);

	return tanh(ksip / (2.0 * Tc)) + tanh(ksim / (2.0 * Tc));
}

__device__ double singlet_uifunt(double x, double Tc, double mi, double gamma0, double tp)
{
	double ksip = ek(x, tp) - mi + g(x, gamma0), ksim = ek(x, tp) - mi - g(x, gamma0);
	double res;

	if (abs(ksim) < 0.0000000000001)
		res = 1.0 / (2.0 * Tc);
	else
		res = tanh(ksim / (2.0 * Tc)) / (2.0 * ksim);

	if (abs(ksip) < 0.0000000000001)
		res += 1.0 / (2.0 * Tc);
	else
		res += tanh(ksip / (2.0 * Tc)) / (2.0 * ksip);

	return res;
}

__device__ double gd(double x, double gamma0)
{
	return gamma0 * (sin(x) * sin(x));
}

__device__ double d2(double x)
{
	return pow(sin(x), 2.0);
}

__device__ double tryplet_uifunt(double x, double Tc, double mi, double gamma0, double tp)
{
	double ksi = ek(x, tp) - mi, eps = 0.000000001, vg = g(x, gamma0);
	double uifunt1 = 0.0, uifunt2 = 0.0, pare = 0.0, uifunt;
	if ((abs(ksi) < eps) && (abs(vg) < eps))
	{
		uifunt1 = 1.0 / (2.0 * Tc);
	}
	else
	{
		if (abs(ksi) < eps)
		{
			uifunt1 = tanh(vg / (2.0 * Tc)) / vg;
			pare = d2(x) - pow(gd(x, gamma0) / vg, 2.0);

			if (abs(vg / Tc) > 100)
			{
				uifunt2 = -1.0 / vg;
			}
			else if (abs(vg / Tc) < eps)
			{
				uifunt2 = 0.0;
			}
			else
			{
				uifunt2 = 1.0 / (Tc + Tc * cosh(vg / Tc)) - sinh(vg / Tc) / (vg + vg * cosh(vg / Tc));
			}

		}
		else
		{
			if (abs(vg) < eps)
			{
				uifunt1 = tanh((ksi) / (2.0 * Tc)) / (ksi);
			}
			else if ((abs(ksi + vg) < eps) || (abs(ksi - vg) < eps))
			{
				uifunt1 = (Tc * tanh(vg / Tc) + vg) / (4.0 * Tc * vg);
				uifunt2 = 0.250 * (tanh(vg / Tc) / vg - 1.0 / Tc);
				pare = d2(x) - pow(gd(x, gamma0) / vg, 2.0);
			}
			else
			{
				uifunt1 = tanh((ksi + vg) / (2.0 * Tc)) / (2.0 * (ksi + vg)) + tanh((ksi - vg) / (2.0 * Tc)) / (2.0 * (ksi - vg));
				uifunt2 = vg / ksi * (tanh((ksi + vg) / (2.0 * Tc)) / (2.0 * (ksi + vg)) - tanh((ksi - vg) / (2.0 * Tc)) / (2.0 * (ksi - vg)));
				pare = d2(x) - pow(gd(x, gamma0) / vg, 2.0);
			}

		}

	}
	uifunt = uifunt1 * d2(x) + uifunt2 * pare;
	return uifunt;
}

__device__ Under_Integral_Func uif_array[] = { sc_uifunc, singlet_uifunt, tryplet_uifunt };
enum UI_Func_Index { i_sc_uifunc = 0, i_singlet_uifunt = 1, i_tryplet_uifunt };

__global__ void integrateKernel(double* result, UI_Func_Index i_f, double a, double h, int n, double* xk, double* ak, double Tc, double mi, double gamma0, double tp) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < n) {
		int j = idx % 5;
		int i = idx / 5;
		double z = ((2.0 * a + h * (2.0 * i + 1.0)) - h * xk[j]) * 0.5;
		result[idx] = ak[j] * ((uif_array[i_f])(z, Tc, mi, gamma0, tp));
	}
}

double sc_integrate1D_gl_gpu(UI_Func_Index i_f, double a, double b, double Tc, double mi, double gamma0, double tp) {
	double h = (b - a) / (n);
	// Host variables
	double* h_result = new double[n * k];

	// Device variables
	double* d_result;
	double* d_xk;
	double* d_ak;
	cudaMalloc((void**)&d_result, n * k * sizeof(double));
	cudaMalloc((void**)&d_xk, k * sizeof(double));
	cudaMalloc((void**)&d_ak, k * sizeof(double));

	// Copy data to device
	cudaMemcpy(d_xk, xk, k * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ak, ak, k * sizeof(double), cudaMemcpyHostToDevice);

	// Launch the kernel
	int num_blocks = (n * k + BLOCK_SIZE - 1) / BLOCK_SIZE;
	integrateKernel << <num_blocks, BLOCK_SIZE >> > (d_result, i_f, a, h, n * k, d_xk, d_ak, Tc, mi, gamma0, tp);

	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed: %s\n", cudaGetErrorString(cudaStatus));
	}
	// Copy the result back to host
	cudaMemcpy(h_result, d_result, n * k * sizeof(double), cudaMemcpyDeviceToHost);

	double result = accumulate(h_result, h_result + n * k, 0.0) * h * 0.5;
	// Cleanup
	cudaFree(d_result);
	cudaFree(d_xk);
	cudaFree(d_ak);
	delete[] h_result;
	return result;
}

double singlet_gap(double Tc, double mi, double gamma0, double Vs, double tp)
{
	return 1.0 - Vs / (2.0 * pi) * sc_integrate1D_gl_gpu(i_singlet_uifunt, 0.0, pi, Tc, mi, gamma0, tp);
}

double sc_occ(double mi, double Tc, double gamma0, double nl, double tp)
{
	return 1.0 - nl - sc_integrate1D_gl_gpu(i_sc_uifunc, 0.0, pi, Tc, mi, gamma0, tp) / (2.0 * pi);
}

typedef double (*Gap_Func)(double, double, double, double, double);


double occ_solve1D_zbr(Gap_Func f, double a, double b, double tol, double T, double gamma0, double nl, double tp)
{
	int itmax = 100;
	double d, r, s, e, p, q, xm, tol1, c, fa, fb, fc;

	fa = f(a, T, gamma0, nl, tp);
	fb = f(b, T, gamma0, nl, tp);

	if (fa * fb > 0.0)
	{
		cout << "occ_zbr err:Takie same znaki!  fa  " << fa << " fb  " << fb << " tc,g0,a,b  " << T << "  " << gamma0 << "  " << a << "  " << b << endl;
		return 0.0;
	}
	c = b;
	fc = fb;
	for (int i = 1; i < itmax; i++)
	{
		if (fb * fc > 0.0)
		{
			c = a;
			fc = fa;
			d = b - a;
			e = d;
		}
		if (abs(fc) < abs(fb))
		{
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0 * delta * abs(b) + 0.5 * tol;
		xm = 0.5 * (c - b);
		if ((abs(xm) < tol1) || (fb == 0)) return b;
		if ((abs(e) > tol1) && (abs(fa) > abs(fb)))
		{
			s = fb / fa;
			if (a == c)
			{
				p = 2.0 * xm * s;
				q = 1.0 - s;
			}
			else
			{
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}
			if (p > 0.0) q = -q;
			p = abs(p);
			if (2.0 * p < min(3.0 * xm * q - abs(tol1 * q), abs(e * q)))
			{
				e = d;
				d = p / q;
			}
			else
			{
				d = xm;
				e = d;
			}
		}
		else
		{
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if (abs(d) > tol1) b += d;
		else
		{
			if (xm > 0.0) b = b + abs(tol1);
			else b = b - abs(tol1);
		}
		fb = f(b, T, gamma0, nl, tp);


	}
	cout << "zbr exeding max iteractions!" << endl;
	return b;
}

struct Result_Pair {
	double T;
	double mi;
};

Result_Pair sc_solve1D_zbr(Gap_Func f, double a, double b, double tol, double gamma0, double Vs, double sxc, double syc, double nl, double tp)
{
	int itmax = 100;
	double d, r, s, e, p, q, xm, tol1, c, fa, fb, fc, eps = 3.0e-8, mi;
	mi = occ_solve1D_zbr(sc_occ, sxc, syc, tol, a, gamma0, nl, tp);
	fa = f(a, mi, gamma0, Vs, tp);
	mi = occ_solve1D_zbr(sc_occ, sxc, syc, tol, b, gamma0, nl, tp);
	fb = f(b, mi, gamma0, Vs, tp);
	Result_Pair res = { b, mi };

	if (fa * fb > 0.0)
	{
		cout << "sc_zbr err:Takie same znaki!  fa  " << fa << " fb  " << fb << " mi,g0,a,b  " << mi << "  " << gamma0 << "  " << a << "  " << b << endl;
		return { 0.0,0.0 };
	}
	c = b;
	fc = fb;
	for (int i = 1; i < itmax; i++)
	{
		if (fb * fc > 0.0)
		{
			c = a;
			fc = fa;
			d = b - a;
			e = d;
		}
		if (abs(fc) < abs(fb))
		{
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0 * eps * abs(b) + 0.5 * tol;
		xm = 0.5 * (c - b);
		if ((abs(xm) < tol1) || (fb == 0)) return res;
		if ((abs(e) > tol1) && (abs(fa) > abs(fb)))
		{
			s = fb / fa;
			if (a == c)
			{
				p = 2.0 * xm * s;
				q = 1.0 - s;
			}
			else
			{
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}
			if (p > 0.0) q = -q;
			p = abs(p);
			if (2.0 * p < min(3.0 * xm * q - abs(tol1 * q), abs(e * q)))
			{
				e = d;
				d = p / q;
			}
			else
			{
				d = xm;
				e = d;
			}
		}
		else
		{
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if (abs(d) > tol1) b += d;
		else
		{
			if (xm > 0.0) b = b + abs(tol1);
			else b = b - abs(tol1);
		}
		mi = occ_solve1D_zbr(sc_occ, sxc, syc, tol, b, gamma0, nl, tp);
		fb = f(b, mi, gamma0, Vs, tp);
		res = { b,mi };
	}

	cout << "zbr exeding max iteractions!" << endl;
	return res;
}

Result_Pair singlet_get_res(double Vs, double nl, double tp, double gamma0, double xt, double yt, double xc, double yc)
{
	double ch, tol = 0.000000001;

	Result_Pair res = sc_solve1D_zbr(singlet_gap, xt, yt, tol, gamma0, Vs, xc, yc, nl, tp);
	/* ch = occ_solve1D_zbr(sc_occ, xc, yc, tol, res.T, gamma0, nl, tp);

	 while(abs(ch-res.mi)>tol)
	 {
		 res=sc_solve1D_zbr(singlet_gap,xt,yt,tol,gamma0,Vs,res.mi *.05,res.mi *1.5,nl,tp);
		 ch=occ_solve1D_zbr(sc_occ,xc,yc,tol,res.T,gamma0,nl,tp);
	 }
	 */
	std::cout << "singlet ch " << res.mi << " t " << res.T << std::endl;
	return res;
}

double tryplet_gap(double Tc, double mi, double gamma0, double Vt, double tp)
{
	const double pi = asin(1.0) * 2.0;
	double temp;

	temp = sc_integrate1D_gl_gpu(i_tryplet_uifunt, 0.0, pi, Tc, mi, gamma0, tp);
	temp = 1.0 - Vt / (2.0 * pi) * (temp);
	return temp;

}

Result_Pair tryplet_get_res(double Vt, double nl, double tp, double gamma0, double xt, double yt, double xc, double yc)
{
	double t, ch, chp, tol = 0.000000001;
	Result_Pair	res = sc_solve1D_zbr(tryplet_gap, xt, yt, tol, gamma0, Vt, xc, yc, nl, tp);

	/*cout<<"t tryplet ch "<<chp<<" t "<<t<<endl;
	t=sc_solve1D_zbr(tryplet_gap,xt,yt,tol,0.0,gamma0,n,k,Vt);
	cout<<"ch tryplet ch "<<chp<<" t "<<t<<endl;
	ch=sc_solve1D_zbr(sc_occ,xc,yc,tol,t,gamma0,n,k,nl);

	chp=0.0;

	while(abs(ch-chp)>tol)
	{
		chp=ch;
		cout<<"t tryplet ch "<<chp<<" t "<<t<<endl;
		t=sc_solve1D_zbr(tryplet_gap,xt,yt,tol,chp,gamma0,n,k,Vt);
		cout<<"ch tryplet ch "<<chp<<" t "<<t<<endl;
		ch=sc_solve1D_zbr(sc_occ,xc,yc,tol,t,gamma0,n,k,nl);
	}

		res[0]=t;
		res[1]=ch;*/
	std::cout << "tryplet ch " << res.mi << " t " << res.T << std::endl;

	return res;
}

struct Density_Result {
	double plus;
	double minus;
};

__device__ Density_Result sc_dos(double o, double g0, double tp)
{
	int npi = 100000;
	double ddos = 1.e-02;
	double step = 2.0 * pi_device / npi;
	double densitytot = 2. * npi;
	double omega, gamma, omp, omm;

	double domega = 2. * ddos;

	double densityp = 0.;
	double densitym = 0.;

	double x = -pi_device;
	for (int j = 0; j < npi; j++)
	{
		omega = ek(x, tp);
		gamma = g(x, g0);
		omp = omega + gamma;
		omm = omega - gamma;
		if (abs(omp - o) < ddos) densityp = densityp + 1.;
		if (abs(omm - o) < ddos) densitym = densitym + 1.;
		x = x + step;
	}
	densitym = densitym / (densitytot * domega);
	densityp = densityp / (densitytot * domega);
	return { densityp,densitym };
}

double fs(double gamma0, double mi, double lambda)
{
	double a = 1 + 0.25 * gamma0 * gamma0 - 0.25 * mi * mi;
	if (a > 0)
	{
		a = -0.5 * mi + lambda * 0.5 * gamma0 * sqrt(a);
	}
	else
	{
		a = -0.5 * mi;
	}
	a = a / (1 + 0.25 * gamma0 * gamma0);
	return a;
}


void sc_tabulate1D(string fname, Result_Pair(*f)(double, double, double, double, double, double, double, double), double a, double b, int N, double Vs, double tp, double nl, double xt, double yt, double xc, double yc)
{
	fstream outfile(fname, fstream::out);
	double g0 = a;
	double h = (b - a) / (N - 1);

	if (!outfile.good())
	{
		cout << "nie otwarty plik!" << endl;
	}

	while (h > 0.0001)
	{
		g0 = g0 + h;
		Result_Pair res = f(Vs, nl, tp, g0, xt, yt, xc, yc);
		if (res.T == 0.0)
		{
			g0 = g0 - h;
			h = h / 2.0;
		}
		else
		{
			outfile << setprecision(10) << g0 << " " << setprecision(10) << res.T << " " << setprecision(10) << res.mi << endl;
		}

	}
}

int main()
{
	double xt = 0.000000001, yt = 10.0, xc = 0.0000001, yc = 5.0;
	double tp = 0.0, nl = 1.2;
	double Vs = 1.2;

	sc_tabulate1D("tryplet.txt", tryplet_get_res, 0.01, 0.5, 100, Vs, tp, nl, xt, yt, xc, yc); //fname,f,a,b,N,Vs,tp,nl,xt,yt,xc,yc 
	return 0;
}
