#define NCS_IMPLEMENTATION
#include "../../NCS.h"

int main()
{
	double xt = 0.000000001, yt = 0.5, xc = 0.0000001, yc = 0.5;
	double tp = 0.0, nl = 1.1;
	double Vs = 2, Vt = 2;

	sc_tabulate1D("tryplet.txt", tryplet_get_res, 0.01, 0.5, 100, Vt, tp, nl, xt, yt, xc, yc);
	return 0;
}
