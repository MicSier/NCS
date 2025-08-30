#define NCS_IMPLEMENTATION
#include "../../NCS.h"
/*
#define LatexrateC_IMPLEMENTATION
#define _CRT_SECURE_NO_WARNINGS
#include "../../LatexrateC.h"

double band(double x)
{
	return ek(x, 0.0);
}

int main()
{
	FILE* doc = create_latex_doc("singlet.tex");

	write_header(doc);
	double xt = 0.000000001, yt = 10.0, xc = 0.0000001, yc = 5.0;
	double tp = 0.0, nl = 1.1;
	double Vs = 2, Vt = 2;

	Named_Function bands = { band, "Energy band"};
	Plot_Config plot_config = default_plot_config();

	calc_and_plot(doc, &bands, plot_config,"");
	write_line(doc, "\\end{document}");
	//sc_tabulate1D("SINGLET.txt", singlet_get_res, 0.01, 10.0, 100, Vs, tp, nl, xt, yt, xc, yc); //fname,f,a,b,N,Vs,tp,nl,xt,yt,xc,yc 
	return 0;
}*/

int main()
{
	double xt = 0.000000001, yt = 1.0, xc = 0.0000001, yc = 5.0;
	double tp = 0.0, nl = 1.1;
	double Vs = 2, Vt = 2;

	sc_tabulate1D("singlet.txt", singlet_get_res, 0.01, 0.5, 10, Vs, tp, nl, xt, yt, xc, yc);
	return 0;
}

