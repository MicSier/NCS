#include <string>
#include <fstream>
#include "gnuplot.h"

using namespace std;

void tabulate1D(string fname,double (*f)(double),double a, double b, int N);

void plot_phase(string nameout,string namein,bool isfromsinglet);

void plot_tc(string nameout,string namein1,string namein2,string namein3,string namein4);

void plot_univ(string nameout,string namein1,string namein2,string namein3,int n);

string path(string name,int n);

void sc_tabulate1D(string fname,double* (*f)(double,double,double,double,double,double,double,double,double*,double*,int,int),double a, double b, int N,double Vs,double tp,double nl,double xt,double yt,double xc,double yc,double xk[],double ak[],int n,int k);

void coupled_tabulate1D(string fname,double* (*f)(double,double,double,double,double,double,double,double,double*,double*,int,int,bool),double a, double b, int N,double Vs,double Vt,double nl,double xt,double yt,double xc,double yc,double xk[],double ak[],int n,int k,bool isfromsinglet);

void coupled_gamma0_fs_tabulate1D(string fname,double* (*f)(double,double,double,double,double,double,double,double*,double*,int,int,bool,double,double),double a, double b, int N,double Vs,double Vt,double nl,double xt,double yt,double xc,double yc,double xk[],double ak[],int n,int k,bool isfromsinglet);

void coupled_universality_tabulate1D(string fname_singlet1,string fname_singlet2,string fname_tryplet1,string fname_tryplet2, string fname);

void sc_dos_tabulate1D(string fname,double* (*f)(double,double,double),double g0,double tp,int N);

void plot_dos(string nameout,string namein,double mi,double V);
