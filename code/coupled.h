
double coupled_fp(double x,double Tc,double mi,double gamma0,double tp);

double coupled_fmdk(double x,double Tc,double mi,double gamma0,double tp);

double coupled_fpdk2 (double x,double Tc,double mi,double gamma0,double tp);

double coupled_gap(double Tc,double mi,double gamma0,int n,int k,double Vs,double Vt,double tp);

double coupled_deltas(double Tc,double mi,double gamma0,int n,int k,double Vs,double tp);

double coupled_simplified_deltas(double Tc,double mi,double gamma0,int n,int k,double Vs,double Vt,double tp);

double coupled_deltat(double Tc,double mi,double gamma0,int n,int k,double Vt,double tp);

double coupled_simplified_deltat(double Tc,double mi,double gamma0,int n,int k,double Vs,double Vt,double tp);

double coupled_simplified_ratio(double Tcs,double mis,double Tct,double mit,double gamma0,int n,int k,double Vs,double Vt,double tp);

double* coupled_get_res(double Vs,double Vt,double nl,double tp,double gamma0,double sxt,double syt,double sxc,double syc,int n,int k,bool isfromsinglet);

double coupled_get_gamma0(double* (*f)(double,double,double,double,double,double,double,double,double,int,int,bool),double Vs,double Vt,double nl,double tp,double sxt,double syt,double sxc,double syc,int n,int k,bool isfromsinglet,double a,double b);

double* eq_gamma0_fs(double gamma0,double Vs,double Vt,double nl,double tp,double sxc,double syc,int n,int k,double pts,double ptt,double pchs,double pcht,bool isfromsinglet);

double* coupled_get_res_both(double Vs,double Vt,double nl,double tp,double gamma0,double sxc,double syc,int n,int k,double pts,double ptt,double pchs,double pcht);

double* eq_gamma0_bz(double gamma0,double Vs,double Vt,double nl,double tp,double sxc,double syc,int n,int k,double pts,double ptt,double pchs,double pcht,bool isfromsinglet);

double* coupled_simplified_gap(double Tc,double mi,double gamma0,int n,int k,double Vs,double Vt,double tp);

double* coupled_get_res_both_simplified(double Vs,double Vt,double nl,double tp,double gamma0,double sxc,double syc,int n,int k,double pts,double ptt,double pchs,double pcht);
