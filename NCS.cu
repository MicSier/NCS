#define NCS_IMPLEMENTATION
#include "NCS.h"

int main()
{
	double xt = 0.000000001, yt = 10.0, xc = 0.0000001, yc = 5.0;
	double tp = 0.0, nl = 1.1;
	double Vs = 2, Vt = 2;

	//sc_tabulate1D("tryplet.txt", tryplet_get_res, 0.01, 0.5, 100, Vs, tp, nl, xt, yt, xc, yc); //fname,f,a,b,N,Vs,tp,nl,xt,yt,xc,yc 
	coupled_tabulate1D("mixed.txt", coupled_get_res, 0.01, 5.0, 10, Vs, Vt, nl, tp, xt, yt, xc, yc, true);
	return 0;
}
